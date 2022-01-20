clear
close all
clc

addpath FUNCTIONS
%%  POPULATION INIT

%SPATIAL FILTERS

% % FILTER 11x11
% filter_file='FILTERS/Gt11B0.0833f0.25.mat';     %Spatial domain - components for 8 orientations of 11x11 Gabor 
% % filter_file='FILTERS/Gt15B0.0833f0.250.mat';     %Spatial domain - components for 8 orientations of 11x11 Gabor 
% k0=0.25;        %SPATIAL FREQUENCY
% samples = 11;
% % samples = 15;
% % RELATIVE BANDWIDTH => B=0.0833;

%FILTER 43x43
filter_file='FILTERS/Gt43B0.0208f0.063.mat';
% filter_file='FILTERS/Gt47B0.0800f0.063.mat';
k0=0.063;               %SPATIAL FREQUENCY  [cycle/pix]
samples = 420;          %STIMULUS DIMENSION [pix]
% samples = 47;
% RELATIVE BANDWIDTH => B=0.0208;

%PHASE SHIFT
% ph_ = pi/4:-pi/4:-5/4*pi; %this choice is related to the offset value of phase_shift (pi/2)
ph_ = pi/2;
n_orient = 8;
ph = repmat(ph_,n_orient,1);  %phase_shift values for each orientation channel 

%OCULAR DOMINANCE
a = 1; 
%TEMPORAL FILTER
Ft_choice = 'gabor'; % 'gabor'; 'exp_decay'; 'adelson_bergen'
%PREFERRED VELOCITY
% v = 0;   
v  = linspace(-1,1,11)*2;
% kk = [-3 -1.5  0.5  1.5 3]; %Preferred velocity with Adelson_Bergen

%NORMALIZATION VALUES
alpha = [1;0];          %no normalization

%Organize the input parameters for the functions
param.spatFreq    = k0;             %"k0";
param.ocDom       = a;              %"a";
param.phShift     = ph;             %"ph";
param.nOrient     = n_orient;       %"n_orient";
param.prefVel     = v;              %"v";
param.tempFilt    = Ft_choice;      %"Ft_Choice";
param.spatialFilt = filter_file;    %"filter_file";
param.samples     = samples;        %"samples";
param.normParam   = alpha;          %"normalization factor";

%% EXAMPLE: POP ACTIVITY ON SET OF MOVING RDS (WITH DIFFERENTS VEL VECTOR)

%STIMULUS DEFINITION
stim.type = 'RDS_tuning';
theta_cell = 0:pi/param.nOrient:pi-pi/param.nOrient;        
[vv,tt] = meshgrid(param.prefVel,theta_cell);
stim.stim_size = size(vv);
stim.truetheta = tt(:);
stim.theta_g = stim.truetheta;
stim.vel_stim = vv(:); 
vx = vv.*cos(tt);
vy = vv.*sin(tt);
stim.vgrat = [vx(:),vy(:)];
stim.dur = 72; %duration in frame
stim.mode = 1;
stim.disp = 0; %dispaly stimulus
theta = [3/2*pi-pi/4,3/2*pi+pi/4];
%SIMULATION
[e,param] = motion_popV1MT(param,stim);
th = 2e-2;
%SAVE DATA
path = 'SIMULATIONS';
OldFolder = cd;
cd(path);
save('myvel_tuning_polarRDS_dur72','e','param','stim','-v7.3')
cd(OldFolder)
%BIOGAUTAMA COMPUTING FOR MT PATTERN RESPONSE
theta_cell_OUT = 0:pi/param.nOrient:pi-pi/param.nOrient;

[xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);

load 'SIMULATIONS/BioGautama/GautamaWieghts88_Plaid.mat'
%Explicit intersection of constraints method to compute weigths
%     W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2/(2*0.25^2));
W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2/(2*0.25^2));
%     W2 = (0.5+0.5*cos(2*pi/4*(xx(:).*cos(tt(:)'-tt(:)) - xx(:)')));
%     W2 = exp(-abs(xx(:).*cos(tt(:)'-tt(:)) - xx(:)')/(0.25)).^2;
W2 = W2 - eye(size(W2));

%DESIRED POPULATION ACTIVITY
[dx] = reshape(stim.vel_stim.*cos(stim.truetheta),8,11);
[dy] = reshape(stim.vel_stim.*sin(stim.truetheta),8,11);

sgmx = .5;
sgmy = .5;
G = exp(-((xx.*cos(tt)-dx).^2/(2*sgmx^2)+...
    (xx.*sin(tt)-dy).^2/(2*sgmy^2)));

pop_resp = squeeze(e(3,:,:,:,:));
sze = size(pop_resp);
%     %NORMALIZATION
%     pop_resp = pop_resp./max(pop_resp,[],4);
pop_resp_BioGautama = reshape(reshape(pop_resp,sze(1)*sze(2),[])*W,sze);
pop_resp_BioGautama2 = reshape((reshape(pop_resp,sze(1)*sze(2),[])*W2'),sze);
%     pop_resp_BioGautama = squeeze(mean(mean(pop_resp_BioGautama(61:end-60,61:end-60,:,:),1),2));
%     pop_resp_BioGautama2 = squeeze(mean(mean(pop_resp_BioGautama2(:,:,:,:),1),2));
%     pop_resp = squeeze(mean(mean(pop_resp(:,:,:,:),1),2));
%POP RESPONSE TO MOVING RDS (HORIZONTALLY AT 2 pix/frame)
figure,plot_pop_response(pop_resp,0,0,param.prefVel)
title('POP RESPONSE')
figure,plot_pop_response(pop_resp_BioGautama,0,0,param.prefVel)
title('POP RESP BIO GAUTAMA')
figure,plot_pop_response(pop_resp_BioGautama2,0,0,param.prefVel)
title('POP RESP BIO GAUTAMA2')
%     figure,plot_pop_response(pop_resp_BioGautama2-pop_resp,0,0,param.prefVel)
figure,plot_pop_response(G,0,0,param.prefVel)
%POP RESPONSE TO MOVING RDS (HORIZONTALLY AT 2 pix/frame)
figure,plot_pop_response(pop_resp,vx,vy,param.prefVel)
title('POP RESP')
figure,plot_pop_response(pop_resp_BioGautama,vx,vy,param.prefVel)
title('POP RESP BIO GAUTAMA')
figure,plot_pop_response(pop_resp_BioGautama2,vx,vy,param.prefVel)
title('POP RESP BIO GAUTAMA2')
%     figure,plot_pop_response(pop_resp_BioGautama2-pop_resp,vx,vy,param.prefVel)

%SELECT data
%     etmp = squeeze(e(:,ceil(samples/2),ceil(samples/2),:,:,:,1,1));
%THRESHOLDING
%     M = max(max(etmp,[],3),[],2);
%     M = repmat(M,1,size(etmp,2),size(etmp,3));
%     etmp_norm = etmp./M;
%     mask = abs(etmp_norm)>th;
%     etmp = etmp.*mask;
%     %NORMALIZATION
%     M = max(etmp,[],3);
%     enormtmp = (etmp+M/2)./M;
%     enormtmp(isnan(enormtmp)) = 0;
%     enormtmp(isinf(enormtmp)) = 1;
%BIOGAUTAMA

%     if strcmp(stim.type,'grat')&&length(stim.theta_g)==2
%         %If we use two separate grat Threshold also that response
%         etmp2 = squeeze(e(:,ceil(samples/2),ceil(samples/2),:,:,1,2,2));
%         surf_motion_pop(etmp2,param)
%         %THRESHOLDING
%         M = max(max(etmp2,[],3),[],2);
%         M = repmat(M,1,size(etmp2,2),size(etmp2,3));
%         etmp_norm = etmp2./M;
%         mask = abs(etmp_norm)>th;
%         etmp2 = etmp2.*mask;
%         %TOTAL RESPONSE
%         etmp = etmp + etmp2;
%     end
%     surf_motion_pop(etmp,param)
%     surf_motion_pop(enormtmp,param)
%     %Display Data
%     polarplot_motion_tun()
% end

% %% BIOGAUTAMA TEST
% if bio_gautama == 1
%     
%     lambda = 0.01; %regularization term
% %     stim = init_stimulus([theta(o),theta(o)+pi/6]);
%     stim.type = 'RDS';
%     stim.dur = 27;
%     stim.disp = 0;
% %     stim.disp = 1;
%     %SIMULATION
%     [e,param] = motion_popV1MT(param,stim);
%     th = 2e-2;
%     etmp = squeeze(e);
%     %THRESHOLDING
%     M = max(max(etmp,[],3),[],2);
%     M = repmat(M,1,size(etmp,2),size(etmp,3),1,1);
%     etmp_norm = etmp./M;
%     mask = abs(etmp_norm)>th;
%     etmp = etmp.*mask;
% end

%%  FUNCTIONS
function stim = init_stimulus(varargin)
    stim.type = 'grat';
% %     TYPE II PLAID
    stim.theta_g = 0:pi/8:pi-pi/8; %true orientation
    stim.vgrat = linspace(-1,1,11)*2;

    stim.dur = 27; %duration in frame
    stim.disp = 0; %set to 1 to show visual stimulus in a figure
end