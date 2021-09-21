function p = plaid(apert_rad,dur,truetheta,vpld,k,vgrat,theta_g,alpha,contrast)
% MOD = PLAID(THETA1, THETA2): creates a plaid object
% PLAID specifies geometry of plaid movements 
% NB: theta1 and theta2 are expressed in rad

% [p.C1,p.R1,p.C1S]=grating_geometry(grat.theta(1));
% [p.C2,p.R2,p.C2S]=grating_geometry(grat.theta(2));
if ~length(k)==2
    error('Spatial frequency must be defined for both gratings!!')
end
if ~length(vgrat)==2
    error('Velocity must be defined for both gratings!!')
end
if ~length(theta_g)==2
    error('Velocity must be defined for both gratings!!')
end
p.apert_rad = apert_rad;
p.dur = dur;
p.truetheta = truetheta;
p.vpld = vpld;
p.k = k;
p.vgrat = vgrat;
p.theta_g = theta_g;
p.alpha = alpha;
p.c = contrast;
% from structure to class...
p = class(p,'plaid');