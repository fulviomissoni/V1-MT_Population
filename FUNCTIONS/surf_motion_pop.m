function surf_motion_pop(e,param)

[neuron,n_or,n_vel]= size(e);
if n_or==8
    theta_pref = linspace(0,pi,n_or+1);
    e = cat(2,e,flip(e(:,1,:),3));
end
if n_or==16
    theta_pref = linspace(0,2*pi,n_or+1);
    e = cat(2,e,e(:,1,:));
end
[s,t] = meshgrid(param.prefVel,theta_pref);
[xq,yq] = pol2cart(t,s);

for j=1:neuron
    tmp = squeeze(e(j,:,:));
    figure,
    surf(xq,yq,tmp,'EdgeColor','none','FaceColor','interp')
    axis off
    view(0,90)
    colorbar
end
