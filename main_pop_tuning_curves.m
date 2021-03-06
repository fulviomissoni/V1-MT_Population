%TUNING CURVES

clear
close all
clc

addpath FUNCTIONS
%%  POPULATION INIT

% % FILTER 11x11
% filter_file='FILTERS/Gt11B0.0833f0.25.mat';     %Spatial domain - components for 8 orientations of 11x11 Gabor 
% % filter_file='FILTERS/Gt15B0.0833f0.250.mat';     %Spatial domain - components for 8 orientations of 11x11 Gabor 
% k0=0.25;        %SPATIAL FREQUENCY
% samples = 11;
% % samples = 15;
% % RELATIVE BANDWIDTH => B=0.0833;

%FILTER 43x43
filter_file='FILTERS/Gt43B0.0208f0.063.mat';
% filter_file='FILTERS/Gt47B0.0210f0.063.mat';
% filter_file='FILTERS/Gt47B0.0270f0.063.mat';
% filter_file='FILTERS/Gt47B0.0300f0.063.mat';
% filter_file='FILTERS/Gt47B0.0400f0.063.mat';
% filter_file='FILTERS/Gt47B0.0800f0.063.mat';
k0=0.063;     %SPATIAL FREQUENCY
samples = 43;
% samples = 47;
% RELATIVE BANDWIDTH => B=0.0208;

% % FILTER 21x21
% filter_file = 'FILTERS/Gt21B0.0417f0.125.mat';
% k0 = 0.125;
% samples = 21;

% ph_ = pi/4:-pi/4:-5/4*pi; %this choice is related to the offset value of phase_shift (pi/2)
ph_ = pi/2;
n_orient = 8;
ph = repmat(ph_,n_orient,1);  %phase_shift values for each orientation channel 
%d_pref = []; %11X11 GABOR
%d_pref = []; %43X43 GABOR

a=0; %OCULAR DOMINANCE
Ft_choice = 'gabor'; % 'gabor'; 'exp_decay'; 'adelson_bergen'
% v = 0;   %Preferred velocity
v  = linspace(-1,1,11)*2;
% v=1;
% kk = [-3 -1.5  0.5  1.5 3]; %Preferred velocity with Adelson_Bergen

%NORMALIZATION VALUES
% alpha = [0,0;1,1];    %my values -> for MID detectors identifies one pixel peak!
alpha = [0,0;1,1];          %no normalization

%Organize the input parameters for the function
param.spatFreq    = k0;             %"k0";
param.ocDom       = a;              %"a";
param.phShift     = ph;             %"ph";
param.nOrient     = n_orient;       %"n_orient";
param.prefVel     = v;              %"v";
param.tempFilt    = Ft_choice;      %"Ft_Choice";
param.spatialFilt = filter_file;    %"filter_file";
param.samples     = samples;        %"samples";
param.normParam   = alpha;          %"normalization factor";
%%
% CHOICES
RDS_Test = 0;
MID_Test = 0;
motion_pop = 0;
bio_gautama = 1;
%% RDS-Test

if RDS_Test==1
    choice = '1D_disp';
%     choice = '2D_disp';
    %For the total simulation one day (at least) is necessary
    RDS_simulation(param,choice)
    %Display Data
    switch(choice)
        case '1D_disp'
            %Disparity Tuning Curves
            plot_disp_tun()
         case '2D_disp'
            %Disparity Tuning Surfaces 
            surf_disp_tun()
    end
end
%% MID-Test - Velocity Tuning Surface
if MID_Test==1
    MID_simulation(param)
    %Display Data
    surf_MID_tun()
end
%% monocular motion direction selectivity, population activity

%Stimulus properties
theta = [3/2*pi-pi/4,3/2*pi+pi/4];
stim = init_stimulus(theta);

if motion_pop==1
    %SIMULATION
    [e,param]= motion_popV1MT(param,stim);
    th = 2e-2;
    %SELECT data
    etmp = squeeze(e(:,ceil(samples/2),ceil(samples/2),:,:,:,1,1));
    %THRESHOLDING
    M = max(max(etmp,[],3),[],2);
    M = repmat(M,1,size(etmp,2),size(etmp,3));
    etmp_norm = etmp./M;
    mask = abs(etmp_norm)>th;
    etmp = etmp.*mask;
%     %NORMALIZATION
%     M = max(etmp,[],3);
%     enormtmp = (etmp+M/2)./M;
%     enormtmp(isnan(enormtmp)) = 0;
%     enormtmp(isinf(enormtmp)) = 1;
    %BIOGAUTAMA

    if strcmp(stim.type,'grat')&&length(stim.theta_g)==2
        %If we use two separate grat Threshold also that response
        etmp2 = squeeze(e(:,ceil(samples/2),ceil(samples/2),:,:,1,2,2));
        surf_motion_pop(etmp2,param)
        %THRESHOLDING
        M = max(max(etmp2,[],3),[],2);
        M = repmat(M,1,size(etmp2,2),size(etmp2,3));
        etmp_norm = etmp2./M;
        mask = abs(etmp_norm)>th;
        etmp2 = etmp2.*mask;
        %TOTAL RESPONSE
        etmp = etmp + etmp2;
    end
    surf_motion_pop(etmp,param)
%     surf_motion_pop(enormtmp,param)
%     %Display Data
%     polarplot_motion_tun()
end

%% BIOGAUTAMA TEST
if bio_gautama == 1

    lambda = 0.01; %regularization term
    for o = 1:length(theta)
%         stim = init_stimulus([theta(o),theta(o)+pi/6]);
        
        %SIMULATION
        [e,param]= motion_popV1MT(param,stim);
        th = 2e-2;
        %SELECT data
        etmp = squeeze(e(:,ceil(samples/2),ceil(samples/2),:,:,:,1,1));
        %THRESHOLDING
        M = max(max(etmp,[],3),[],2);
        M = repmat(M,1,size(etmp,2),size(etmp,3));
        etmp_norm = etmp./M;
        mask = abs(etmp_norm)>th;
        etmp = etmp.*mask;
        EC1 = squeeze(etmp(4,:,:));
        EC1 = EC1(:);
        G = zeros(8,11);
        G(4:6,7:10) = fspecial('gaussian',[3,4]);
        G = G(:);
        W = G*EC1'*inv(lambda*eye(88,88)+EC1*EC1');
        figure, imagesc(W)
        pause
        figure, imagesc(W)
        pause
    end
end

%%  FUNCTIONS
function stim = init_stimulus(theta)
    stim.type = 'plaid';
    % TYPE II PLAID
%     stim.theta_g = [3/2*pi-pi/6,3/2*pi-pi/3]; %true orientation
%     % stim.theta_g = 0;
%     stim.truetheta =  3/2*pi; %true orientation
    % TYPE I PLAID
    stim.theta_g = [theta(1),theta(2)]; %true orientation
    % stim.theta_g = 0;
    stim.truetheta =  3/2*pi; %true orientation
    stim.vpld = [1];
    % stim.vgrat = [2,1];
    if size(stim.theta_g,2)==2
        stim.vgrat = stim.vpld.*[cos(stim.truetheta-stim.theta_g(1)), cos(stim.truetheta-stim.theta_g(2))];
        stim.vgrat = round(stim.vgrat,5,"decimals");
    else
        stim.vgrat = stim.vpld;
    end
    % stim.vgrat = 0;
    stim.dur = 43; %duration in frame
    stim.mode = 1;
    stim.disp = 0; %set to 1 to show visual stimulus in a figure
end