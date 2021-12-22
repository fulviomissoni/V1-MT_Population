function [varargout] = motion_popV1MT(param,stim)
% % % % % % % DEBUG % % % % % % %
debug = false;

file_name='mono_motion_pop';
samples=param.samples;
k0=param.spatFreq;
dur = stim.dur;
if ~isfield(stim,'type')
    error('Define type of the stimulus!! Options are: grat, plaid, or RDS')
end
fieldNames{1} = 'dur'; 
% fieldNames{2} = 'truetheta'; fieldNames{3} = 'vpld';
fieldNames{2} = 'vgrat'; fieldNames{3} = 'theta_g';
argCheck(stim,fieldNames);
%% Stimulus definition 
stimuli = ["plaid","grat","RDS","shift_grat"];
if ~sum(matches(stimuli,stim.type))
    error('Select from the possible stimuli:\n %s, %s, %s',stimuli{1},stimuli{2},stimuli{3})
end
    
switch stim.type
    case 'grat'
        II = sinGrating(samples,samples,dur,[0,0],stim.vgrat,k0,stim.theta_g);
    case 'plaid'
        if ~isfield(stim,'pl_type')
            pl_type = 1;
        else 
            pl_type = stim.pl_type;
        end
        %plaid object
        arg.dur = dur;                      %aperture size      [pixs]
        arg.apert_rad = ceil(samples/2)+2;  %aperture size      [pixs]
        arg.truetheta = stim.truetheta;     %true orientation   [rad]
        arg.vpld = stim.vpld;               %velocity amplitude [pixs/frame]
        arg.k = [k0,k0];                    %spatial freq       [cycle/pix]
        arg.vgrat = stim.vgrat;             %gratings vel       [pixs/frame]
        arg.theta_g = stim.theta_g;         %gratings orient    [rad]
        arg.alpha = 0.5;                    %Alpha channel for transparency
        arg.contrast = [0.5,0.5];         %Contrast of two gratings
        arg.mode = stim.mode;               %stimulus implementation algorithm
        arg.pl_type = pl_type;         %plaid type
        %define plaid object
        pld = plaid(arg);                          
        %generate plaid stimulus
        II = generate_plaid(pld,[],arg.mode);
    case 'RDS'
        %make moving dots
        for i=1:dur
            if i==1
                II{1}(:,:,1) = myRDS(stim.vgrat(1),1,1,1,samples,samples);
            else
                II{1}(:,:,i) = circshift(II{1}(:,:,i-1),stim.vgrat(1),2);
            end
        end
    case 'shift_grat'
        II = sinGrating(samples,samples,dur,[0,0],stim.vgrat(1),k0,stim.theta_g);
        for i=1:floor(1/(2*k0))+1
           II(i+1) = sinGrating(samples,samples,dur,[i,0],stim.vgrat(1),k0,stim.theta_g);
        end
end
if stim.disp == 1
    prompt = 'Press any number to start visualization of visual stimulus\n';
    start = input(prompt);
    if ~isempty(start)
        for i=1:dur
            imagesc(squeeze(II{1}(:,:,i)))
            colorbar
            pause(0.1)
        end
    end
end

%% motion-in-depth descriptors analysis - tuning curves
n_vel = length(param.prefVel);
n_orient = param.nOrient;
% [IIv,IIo] = size(II);
% ph_shift = 0;
% parameters{3,1} = ph_shift;
e=zeros(4,size(II{1},1),size(II{1},2),n_orient,n_vel,size(param.phShift,2),size(II,1),size(II,2)); 
TUN = 0;
% %ocular dominance is 0 or 1 'cause is monocular test
param.ocDom = 1; %force this value
for i=1:size(II,1)
    for j=1:size(II,2)
        %select input
        I = cat(4,zeros(size(II{i,j})),II{i,j});
        [MT,EC21,EC22,EC1,ES] = pop_flow_V1MT(I,param,TUN);
        if debug
            e(8,:,:,:,:,:,i,j) = ES(1,:,:,:,:,:);
            e(7,:,:,:,:,:,i,j) = ES(2,:,:,:,:,:);
            e(6,:,:,:,:,:,i,j) = ES(3,:,:,:,:,:);
            e(5,:,:,:,:,:,i,j) = ES(4,:,:,:,:,:);
        end
        e(4,:,:,:,:,:,i,j) = EC1{1};
        e(3,:,:,:,:,:,i,j) = EC1{2};
        e(2,:,:,:,:,:,i,j) = EC21; 
        e(1,:,:,:,:,:,i,j) = MT; 
        fprintf('Pop-activity, stimulus %d \n',size(II,2)*(i-1)+j);
        clear I
    end
end
varargout{1} = e;
varargout{2} = param;
%Save data in SIMULATIONS Directory
path = 'SIMULATIONS/mono';
OldFolder = cd;
cd(path);
save(file_name,'e','param','-v7.3')
cd(OldFolder)
end

function argCheck(stimulus,fieldName)
j=1;
err = [];
for i=1:length(fieldName)
    if ~isfield(stimulus,fieldName{i})
        err(j) = i;
        j=j+1;
    end
end
if ~isempty(err)
    error('Following stimulus properties are not defined:\n %s %s',fieldName{err})
else
end
end
