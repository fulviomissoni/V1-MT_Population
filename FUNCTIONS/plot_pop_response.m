
function plot_pop_response(response,rx,ry,s)

[rx,ry] = meshgrid(rx,ry);

[theta,ro] = cart2pol(rx,ry);

ro = ro/(2*max(ro(:)));

for j=1:size(rx,2)
    for i = 1:size(rx,1)
       posx = ro(i,j)*cos(theta(i,j))+0.5;
       posy = ro(i,j)*sin(theta(i,j))+0.5;
        
        %         axes('Position',[0.075+rem((i-1),8)*0.23/2 0.1+floor((i-1)/8)*0.45/2 0.2/2 0.35/2]);
        %         axes('Position',[0.075+0.5*(j-1)*0.9/9 0.1+(i-1)*0.9/8 0.9/9-0.01 0.9/8-0.01]);
        axes('Position',[posx posy 0.09 0.09]);
        
        surf_polar(response(:,:,i,j),0,90,0,s),colormap(hot(256))
%         colorbar
        set(gca,'Color','none')
        xlim([-3 3]);
        ylim([-3 3]);
        grid off
        if i~=2||j~=1
            set(gca,'xtick',[],'ytick',[]);
            axis off
        else
            xlabel('dx')
            ylabel('dy')
        end

    end
end