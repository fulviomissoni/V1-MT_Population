clear 
close all
clc

addpath FUNCTIONS
%% load data and organize them in cell array
% lambda = ["0","1e-06","1e-05","00001","00001","0001","001"];
lambda = "1e-03";
% block = 'vel_tuning_All_PlaidII_lambda';
block = 'vel_tuning_PlaidII_lambda';
mypath = 'SIMULATIONS\';
filename = [mypath,block,convertStringsToChars(lambda(1)),'_difContrasts.mat'];
% filename = [mypath,block,convertStringsToChars(lambda(1)),'.mat'];
load(filename)
%size(e) is known
nTheta_g = [6,7];
nVel = [5]; %vel stim tested -2pix/frames and -1.2 pix/frames
nContr_g = [1]; %contrast stim tested are grat1 = [0,0.2,0.4]; grat2 = [1,0.8,0.6];
% diff_contrast = [0:0.1:1];
tmp = squeeze(e(3,:,:,:));
theta_cell_OUT = 0:pi/param.nOrient:pi-pi/param.nOrient;
[xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);
load([mypath,'GautamaWeights88_Plaid.mat'])
W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2/(2*0.25^2)); % my weights
W2 = W2-W2.*eye(size(W2));
sze = [param.nOrient,numel(param.prefVel),nTheta_g(1),nTheta_g(2),nVel,nContr_g];
% sze = [param.nOrient,numel(param.prefVel),nTheta_g(1),nTheta_g(2),param.nOrient,numel(param.prefVel)];
pop_respW{1} = reshape(W*reshape(tmp,sze(1)*sze(2),[]),sze);
% pop_respW{1} = pop_respW{1}/max(pop_respW{1},[],"all");
pop_respW{2} = reshape(W2*reshape(tmp,sze(1)*sze(2),[]),sze);
% pop_respW{2} = pop_respW{2}/max(pop_respW{2},[],"all");

% for i=1:numel(lambda)
%     filename = [mypath,block,convertStringsToChars(lambda(i)),'_difContrasts.mat'];
%     load(filename)
%     tmp = squeeze(e(3,:,:,:));
%     pop_respW1{i} = reshape(W*reshape(tmp,sze(1)*sze(2),[]),sze);
%     pop_respW2{i} = reshape(W2*reshape(tmp,sze(1)*sze(2),[]),sze);
% end
%% Mexican hat weights
sigma_r = 0.3;
sigma_t = 0.3;
K=3;
[xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);
tt = tt + pi*(xx<0);
xx = abs(xx);

% [dx] = stim.vel_stim.*cos(stim.truetheta);
% [dy] = stim.vel_stim.*sin(stim.truetheta);
dx = xx(:).*cos(tt(:));
dy = xx(:).*sin(tt(:));
MX = [];
% for l=1:size(dx,1)
%     X = ( (xx.*cos(tt) - dx(l)).^2 + (xx.*sin(tt) - dy(l)).^2 ) / sigma^2;
    X = ((xx(:) - xx(:)').^2) / sigma_r^2 + ((tt(:) - tt(:)').^2) / sigma_t^2;
%     mx = 1/(pi*sigma^4).*( 1 - 1/2.*(X)).*exp(-X/2);
    mx = 1/(2*pi*sigma_r*sigma_t)*exp(-X/2) - 1/(2*pi*K^2*sigma_r*sigma_t)*exp(-X/(2*K^2));
%     mx = exp(-X/2);

%     MX = [MX mx(:)]; 
    MX = mx;
% end
th = 1e-2;
MX = MX./max(MX,[],'all');
CT = reshape(pop_respW{2},sze(1)*sze(2),[]);
%normalization data for each stimulus is different
% CT = (CT./max(CT,[],1))*10-5;
CT = CT./max(CT,[],1);
CT(CT<th) = 0;
k = 9;
for iteration=3:10
%     text = ['pop_respW',num2str(iteration),'{1} = reshape(MX*CT,sze);'];
%     eval(text)
    CT = MX*CT;
%     CT = exp(CT);
%     CT = (1 + exp(CT));
    CT(CT<0) = 0;
    CT = (CT./max(CT,[],1));
        CT = 1 ./(1+exp(-k*(CT-0.7)));
%     CT = (CT./max(CT,[],1))*10-5;
%     CT(CT<=th) = 0;
    pop_respW{iteration} = reshape(CT,sze);
%     pop_respW{iteration} = reshape((CT./max(CT,[],1))*.5+.5,sze);
%     text = ['CT = reshape(pop_respW',num2str(iteration),'{1},sze(1)*sze(2),[]);'];
%     eval(text)
end
% CT = MX*((CT./max(CT,[],1))*.5+.5);
% CT = MX*(CT);
% CT(CT<th) = 0;
% CT = (CT./max(CT,[],1))*.5+.5;
% CT = CT./max(CT,[],1);
% pop_respW{iteration+1} = reshape(CT,sze);
%% visualize pop activity
for i=1:iteration
    tmp = pop_respW{i}(:,:,1,1,1,1);
    if i<=2
        tmp = tmp./max(tmp,[],'all');
    end
    figure, plot_pop_response(tmp,0,0,param.prefVel)
    colorbar
end
%% visualize error percetange of orientation estimate with contrast difference increasing
truetheta = 0;
% mytheta = 3;
% theta_pop = 0:pi/param.nOrient:pi-pi/param.nOrient;
[theta1,theta2] = meshgrid(linspace(-3*pi/8,3*pi/8,7));
theta_g = [theta1(:),theta2(:)];
theta_g(1:8:end,:) = [];
theta_g1 = reshape(theta_g(:,1),6,7);
theta_g2 = reshape(theta_g(:,2),6,7);
for mytheta=1:numel(truetheta)
    for indW=1:iteration
%         for sel = 1:numel(lambda)
        %select pop_response for one one plaid velocity
        mypop_resp = reshape(pop_respW{indW},sze(1),sze(2),sze(3)*sze(4),sze(5),sze(6));
        figure, 
        ind = 1:49;
        ind(1:8:49) = [];
        for i=1:6*7
            %select one combination of theta_g
%             tmpmypop_resp = squeeze(mypop_resp(:,:,i,mytheta,:));
            tmpmypop_resp = squeeze(mypop_resp(:,:,i,:));
            %salvami, per ogni valore di diff di contrasto, per quale orientamento c'è
            %la risposta massima
            [rho, theta] = meshgrid(param.prefVel,linspace(0,pi,size(tmpmypop_resp,1)+1));
            x = linspace(min(param.prefVel),max(param.prefVel),101);
            [Xout, Yout] = meshgrid(x);
            % [thetaout rhoout] = cart2pol(Xout,Yout);
            thetaout = atan(Yout./Xout);
            thetaout = thetaout + (thetaout<0)*pi;
            thetaout(isnan(thetaout)) = 0;
            rhoout = Xout.*cos(thetaout) + Yout.*sin(thetaout);
            for c=1:5
                tmp = tmpmypop_resp(:,:,c);
                %interpolo per avere una popolazione più fitta
                tmpOut = interp2(rho,theta,[tmp; fliplr(tmp(1,:))],rhoout,thetaout);
%                 figure, plot_pop_response(tmp,0,0,param.prefVel);
%                 [i,indW]
%                 pause(0.15)
%                 close
                [m,tmpindOr] = max(tmpOut); 
                [m,tmpindVel] = max(m);
                indOr(c) = tmpindOr(tmpindVel);
                indVel(c) = tmpindVel;
    %             plot_pop_response(tmpmypop_resp(:,:,c),0,0,param.prefVel);
            end
            % thetaout2 = linspace(0,2*pi,  );
            thetapop = mod(atan2(x(indOr),x(indVel)),2*pi);
%             subplot(7,7,ind(i)), plot(diff_contrast, abs(rad2deg( repmat(truetheta,[1 11]) - (thetapop) )))
            newtruetheta = repmat(truetheta(mytheta),[1 5]) + pi*(param.prefVel(1:nVel)<0);
            subplot(7,7,ind(i)), plot(param.prefVel(1:nVel), abs(rad2deg( newtruetheta - thetapop )))
            ylim([0 50])
            xlim([-2 0])
    %         figure,plot(diff_contrast, rad2deg( repmat(truetheta,[1 6]) - (thetapop) ))
            
        end
    end
    rad2deg(truetheta(mytheta))

%     pause
%     close all
end