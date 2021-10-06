%TUNING CURVES

clear
% close all
% clc

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
v = 0;   %Preferred velocity
% v  = linspace(-1,1,11)*4;
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
motion_pop = 1;
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
stim.type = 'shift_grat';
% stim.theta_g = [pi/4,pi/4+pi/6]; %true orientation
stim.theta_g = 0;
% stim.truetheta =  pi/2; %true orientation
% stim.vpld = [2];
% stim.vgrat = [2,1];
% stim.vgrat = stim.vpld.*[cos(stim.truetheta-stim.theta_g(1)), cos(stim.truetheta-stim.theta_g(2))];
% stim.vgrat = round(stim.vgrat,2,'significant');
stim.vgrat = 0;
stim.dur = 27; %duration in frame
stim.mode = 2;
if motion_pop==1
    %SIMULATION
    [e,param]= motion_popV1MT(param,stim);
    th=2e-2;
    %DISPLAy data
    etmp = squeeze(e(:,ceil(samples/2),ceil(samples/2),:,:,:,1,1));
    %THRESHOLDING
    M = max(max(etmp,[],3),[],2);
    M = repmat(M,1,size(etmp,2),size(etmp,3));
    etmp_norm = etmp./M;
    mask=abs(etmp_norm)>th;
    etmp=etmp.*mask;
    %NORMALIZATION
    M = max(etmp,[],3);
    enormtmp = (etmp+M/2)./M;
    enormtmp(isnan(enormtmp)) = 0;
    enormtmp(isinf(enormtmp)) = 1;
    if strcmp(stim.type,'grat')&&length(stim.theta_g)==2
        %If we use two separate grat Threshold also that response
        etmp2 = squeeze(e(:,ceil(samples/2),ceil(samples/2),:,:,1,2,2));
%         surf_motion_pop(etmp2,param)
        M = max(max(etmp2,[],3),[],2);
        M = repmat(M,1,size(etmp2,2),size(etmp2,3));
        %NORMALIZATION
        etmp_norm = etmp2./M;
        mask=abs(etmp_norm)>th;
        etmp2=etmp2.*mask;
        %TOTAL RESPONSE
        etmp = etmp + etmp2;
    end
%     surf_motion_pop(etmp,param)
%     surf_motion_pop(enormtmp,param)
%     %Display Data
%     polarplot_motion_tun()
end


