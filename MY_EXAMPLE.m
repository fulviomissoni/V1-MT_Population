clear 
close all
clc

addpath FUNCTIONS 

%%  POPULATION INIT

% FILTER 11x11
filter_file='FILTERS/Gt11B0.0833f0.25.mat';     %Spatial domain - components for 8 orientations of 11x11 Gabor 
k0=0.25;        %SPATIAL FREQUENCY
samples = 11;
% RELATIVE BANDWIDTH => B=0.0833;

% FILTER 43x43
% filter_file='FILTERS/Gt43B0.0208f0.063.mat';
% k0=0.063;     %SPATIAL FREQUENCY
% samples = 43;
% RELATIVE BANDWIDTH => B=0.0208;

ph_ = -pi/4:pi/4:5/4*pi; %this choice is related to the offset value of phase_shift (pi/2)
n_orient = 8;
ph = repmat(ph_,n_orient,1);  %phase_shift values for each orientation channel 
%d_pref = []; %11X11 GABOR
%d_pref = []; %43X43 GABOR

a=0.2; %OCULAR DOMINANCE
Ft_choice='exp_decay'; % 'gabor'; 'exp_decay'; 'adelson_bergen'
% v = [0.4 0.2 -1];   %Preferred velocity
v  = 1; 
% v  = [+1 +0.6 +0.4 0.3 0.2 0 -0.2 -0.3 -0.4 -0.6 -1];
% v  = [-1 -0.4 0 0.4 1];;
% kk = [-3 -1.5  0.5  1.5 3]; %Preferred velocity with Adelson_Bergen
% filters

%NORMALIZATION PARAMETERS
% alpha = [1,1;1e-3,1e-3];   %best parameters for MID analysis but not for disp analysis
alpha = [0,0;1,1];

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

TUN = 0 ;
%% IMAGE LIST
root='IMAGES/';

img_type='.mat';

image='video_slant_12_tilt_180';    %Binocular retinal images from Retinal Project Simulator (Slanted Plane)      
%% TEST
load([root image img_type])
I=cat(4,IIL,IIR);

tic
    [MT, EC21, EC22] = pop_flow_V1MT(I,param,TUN);
toc

[dumb,dumb2,dd,n_frames,vv,ph_n] = size(MT);
size(MT)

%% VISUALIZATION OF MAP RESPONSE

%Phase-Shift Channel
for ksi=1:ph_n
    %Velocity channel
    for v=1:vv
        %Direction/Orientation Channel
%         for d=1:dd
        d=1;
            figure, 
            tmp = squeeze(MT(:,:,d,1,v,ksi));
            imagesc(tmp)
            title(['MT direction=',num2str(d),' velocity=',num2str(v),'phase-shit=',num2str(ksi)])
            figure,
            imagesc(squeeze(EC21(:,:,d,1,v,ksi)));  
            title(['EmtL direction=',num2str(d),' velocity=',num2str(v),'phase-shit=',num2str(ksi)])
            figure, 
            imagesc(squeeze(EC22(:,:,d,1,v,ksi)));
            title(['EmtR direction=',num2str(d),' velocity=',num2str(v),'phase-shit=',num2str(ksi)])
            pause
%         end
    end
    pause
    close all
end