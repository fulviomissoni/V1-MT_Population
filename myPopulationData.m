clear variables
close all

addpath FUNCTIONS
%% INITIALIZE
% % % POP PARAMETERS % % % 
% % % %%%%%%%%%%%%%% % % %

% -) SPATIAL FILTER 

% FILTER 11x11
ph_shift_file='FILTERS/myph_shift_PH7_F0.063.mat';
filter_file='FILTERS/Gt11B0.0833f0.25.mat';
k0=0.25;

% FILTER 43x43
% filter_file='FILTERS/Gt43B0.0208f0.063.mat';
% ph_shift_file='FILTERS/myph_shift_PH7_F0.063.mat';
% k0=0.063;           % Spatial Frequency
n_orient = 8;       % Orientation of the spatial filters

% -) TEMPORAL FILTER: 1 -> Gabor; 2 -> exp decay; 3 -> Bio-Filter;
Ft_choice=2;
a=0.2;              % Ocular Dominance

% PARAMETERS NOT USED
energy_th = 2e-9;   %energy thresholding -Not used-> to review!!! 
ori_thr = 2e-9;     %energy thresholding along orient. parameter -Not used-> to review!!!

% % % INPUT DATA % % % 
% % % %%%%%%%%%% % % %
root='IMAGES/Input/images';
img_type='.mat';

slant=4:4:48;
tilt=30:30:360;
OldFolder=cd;
gaze=7;
cont=1;
sy=43; sx=43; n_orient=8; phase_num=7;
n_vel=3;
TUN = 0;
% E_ALL = zeros(sy*sx*n_orient*phase_num,2,length(slant),length(tilt),gaze);
% MT_ALL = zeros(sy*sx*n_orient*phase_num*n_vel,length(slant),length(tilt),gaze);

%%
for tt=1:length(tilt)
    for ss=1:length(slant)
        for kk=1:7
            image=['/video_',num2str(cont),'_gaze_',num2str(kk),'_slant_',num2str(slant(ss)),'_tilt_',num2str(tilt(tt))];
            % image='/video_1_gaze_1_slant_4_tilt_30';
            load([root image img_type])
            [sx,sy,n_frames]=size(IIL);
            % clear disparity_map
            %% Population Response
            I = cat(4,IIL,IIR);
            tic
            [MT, EmtL, EmtR] = MODmyctf_pop_flow_V1MT(I,energy_th,ori_thr,n_orient,ph_shift_file,filter_file,k0,Ft_choice,a,TUN);
%             E = energy_ctf_pop_disparity(II,n_scales,n_filters,energy_th,ori_thr,ph_shift_file,filter_file);         
%             E_ALL(:,:,ss,tt,kk) = cat(2,EmtL(:),EmtR(:),MT(:));
%             if cont==1
%                 [sy,sx,n_orient,n_frames,v,phase_num] = size(MT);
%             end
%             MT = MT(:,:,:,end-n_frames+1:end,:,:);
            MT_ALL(:,ss,tt,kk) = MT(:);
            EmtL = squeeze(EmtL(:,:,:,:,3,:));  %STATIC disparity detectors-> I select only one velocity
            EmtR = squeeze(EmtR(:,:,:,:,3,:));  %STATIC disparity detectors-> I select only one velocity
            E_ALL(:,:,ss,tt,kk) = cat(2,EmtL(:),EmtR(:));
            toc
            fprintf('%d\n',cont)
            cont=cont+1;        
            clear MT EmtL EmtR
        end
    end
    save('E_ALL','E_ALL','-v7.3')
    save('MT_ALL','MT_ALL','-v7.3')
end