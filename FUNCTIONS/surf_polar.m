function surf_polar(I,az,el,normal,disp)

if nargin<2
    az=30;
    el=30;
    normal=1;
end

% I = matrix containig the value of each point in polar coordinates: row =
% orientations, cols = modulus
if normal==1
    I=I./abs(repmat(I(:,floor(size(I,2)/2)+1),[1 size(I,2)]));
end

% [rho theta]=meshgrid(linspace(-1,1,size(I,2)),linspace(0,pi,size(I,1)+1));
[rho theta]=meshgrid(disp,linspace(0,pi,size(I,1)+1));

X=rho.*cos(theta); Y=rho.*sin(theta);

surf(X,Y,[I; fliplr(I(1,:))]); %axis normal
shading interp
fvc=surf2patch(0.9*X,0.9*Y,-1+zeros(size(X)));
% patch(fvc);
% % xlabel('dx')
% % ylabel('dy')
% 
view(az,el);

flag=0;

end
