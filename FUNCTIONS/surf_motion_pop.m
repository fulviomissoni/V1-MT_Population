function surf_motion_pop(e,param)
[neuron,n_or,n_vel]= size(e);
mytitle = cell(neuron);
mytitle{1} = "MT";
mytitle{2} = "C22";
mytitle{3} = "C14";
mytitle{4} = "C13";
mytitle{5} = "S8";
mytitle{6} = "S7";
mytitle{7} = "S6";
mytitle{8} = "S5";



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
    title(mytitle{j})
end
