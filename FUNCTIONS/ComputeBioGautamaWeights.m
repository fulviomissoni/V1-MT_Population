close all
clear
clc

load 'SIMULATIONS\vel_tuning.mat'
% e = etmp;
theta_cell_OUT = 0:pi/param.nOrient:pi-pi/param.nOrient;

[xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);
[dx,dy] = meshgrid(param.prefVel,param.prefVel);

% sgmx = 1;
% sgmy = 1;
sgmx = 1;
sgmy = 1;

CT = squeeze(e(3,:,:,:,:));
sze = size(CT);
CT = reshape(CT,size(CT,1)*size(CT,2),[]);

for k=1:size(dx,2)
    k
    for l=1:size(dx,1)
        
        g = exp(-((xx.*cos(tt)-dx(l,k)).^2/(2*sgmx^2)+...
            (xx.*sin(tt)-dy(l,k)).^2/(2*sgmy^2)));

        G(:,l+(k-1)*size(dx,2)) = g(:);
        
    end
end

F = CT*CT';
R = CT*G';
% t = 10000;
t = 1000;

% K = (t*eye(size(F))+F);

W = pinv(t*eye(size(F))+F)*R;

G = reshape(G,sze);

CT2 = W*CT;
CT = reshape(CT,sze);

CT2 = reshape(CT2,sze);

figure,
plot_pop_response(CT,param.prefVel,param.prefVel,param.prefVel)
figure,
plot_pop_response(G,param.prefVel,param.prefVel,param.prefVel)
figure,
plot_pop_response(CT2,param.prefVel,param.prefVel,param.prefVel)

