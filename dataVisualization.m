clear 
close all
clc

addpath FUNCTIONS
%% load data and organize them in cell array
lambda = ["0","1e-06","1e-05","00001","00001","0001","001"];
block = ['vel_tuning_PlaidII_lambda'];
mypath = 'SIMULATIONS\';
filename = [mypath,block,convertStringsToChars(lambda(1)),'_difContrasts.mat'];
load(filename)
%size(e) is known
nTheta_g = [6,7];
nVel = [2]; %vel stim tested -2pix/frames and -1.2 pix/frames
nContr_g = [11]; %contrast stim tested are grat1 = [0,0.2,0.4]; grat2 = [1,0.8,0.6];
diff_contrast = [0:0.1:1];
tmp = squeeze(e(3,:,:,:));
theta_cell_OUT = 0:pi/param.nOrient:pi-pi/param.nOrient;
[xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);
W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2/(2*0.25^2)); % my weights
sze = [param.nOrient,numel(param.prefVel),nTheta_g(1),nTheta_g(2),nVel,nContr_g];

pop_resp{1} = reshape(W2*reshape(tmp,sze(1)*sze(2),[]),sze);
for i=2:numel(lambda)
    filename = [mypath,block,convertStringsToChars(lambda(i)),'_difContrasts.mat'];
    load(filename)
    tmp = squeeze(e(3,:,:,:));
    pop_resp{i} = reshape(W2*reshape(tmp,sze(1)*sze(2),[]),sze);
end
%% visualize error percetange of orientation estimate with contrast difference increasing
truetheta = pi;
% theta_pop = 0:pi/param.nOrient:pi-pi/param.nOrient;
[theta1,theta2] = meshgrid([linspace(-3*pi/8,3*pi/8,7)]);
theta_g = [theta1(:),theta2(:)];
theta_g(1:8:end,:) = [];
theta_g1 = reshape(theta_g(:,1),6,7);
theta_g2 = reshape(theta_g(:,2),6,7);
for sel = 1:numel(lambda)
    % Iout = interp2(rho,theta,[I; fliplr(I(1,:))],rhoout,thetaout);
    % theta_pop(1) = -pi;
    %select pop_response for one one plaid velocity
    mypop_resp = squeeze(pop_resp{sel}(:,:,:,:,1,:));
    % mypop_resp = permute(mypop_resp,[5,6,1,2,3,4]);
    % mypop_resp = reshape(mypop_resp,3,3,[]);
    % %select only diag of all 3x3 matrices
    % mypop_resp = mypop_resp(logical(repmat(eye(size(mypop_resp(:,:,1))),[1 1 size(mypop_resp,3)])));
    % mypop_resp = reshape(mypop_resp,3,8,11,6,7);
    % %rimetto le dimensioni in ordine corretto
    % mypop_resp = permute(mypop_resp,[2,3,4,5,1]);
    mypop_resp = reshape(mypop_resp,8,11,6*7,11);
    figure, 
    ind = 1:49;
    ind(1:8:49) = [];
    for i=1:6*7
        %select one combination of theta_g
        tmpmypop_resp = squeeze(mypop_resp(:,:,i,:));
        %salvami, per ogni valore di diff di contrasto, per quale orientamento c'è
        %la risposta massima
        [rho, theta] = meshgrid(param.prefVel,linspace(0,pi,size(tmp,1)+1));
        X = rho.*cos(theta); Y = rho.*sin(theta);
        x = linspace(min(param.prefVel),max(param.prefVel),101);
        [Xout, Yout] = meshgrid(x);
        % [thetaout rhoout] = cart2pol(Xout,Yout);
        thetaout = atan(Yout./Xout);
        thetaout = thetaout + (thetaout<0)*pi;
        thetaout(isnan(thetaout)) = 0;
        rhoout = Xout.*cos(thetaout) + Yout.*sin(thetaout);
        for c=1:11
            tmp = tmpmypop_resp(:,:,c);
            %interpolo per avere una popolazione più fitta
            tmpOut = interp2(rho,theta,[tmp; fliplr(tmp(1,:))],rhoout,thetaout);
            [m,tmpindOr] = max(tmpOut); 
            [m,tmpindVel] = max(m);
            indOr(c) = tmpindOr(tmpindVel);
            indVel(c) = tmpindVel;
%             plot_pop_response(tmpmypop_resp(:,:,c),0,0,param.prefVel);
        end
        % thetaout2 = linspace(0,2*pi,  );
        thetapop = mod(atan2(x(indOr),x(indVel)),2*pi);
        subplot(7,7,ind(i)), plot(diff_contrast, abs(rad2deg( repmat(truetheta,[1 11]) - (thetapop) )))
        ylim([0 70])
        xlim([0 1])
%         figure,plot(diff_contrast, rad2deg( repmat(truetheta,[1 6]) - (thetapop) ))
        
    end
end
