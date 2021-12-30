theta_cell_OUT = 0:pi/param.nOrient:pi-pi/param.nOrient;

[xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);
[xxm,ttm] = meshgrid(xx(:),tt(:));

W2 = 1-abs(xx(:).*cos(tt(:)-ttm') - xxm)/4;
pop_resp = squeeze(e(3,:,:,:,:));
sze = size(pop_resp);
%     pop_resp = pop_resp./max(pop_resp,[],4);
% pop_resp_BioGautama = reshape(reshape(pop_resp,sze(1)*sze(2),[])*W,sze);
% pop_resp_BioGautama2 = reshape((reshape(pop_resp,sze(1)*sze(2),[])*W2),sze);
% % pop_resp_BioGautama = squeeze(mean(mean(pop_resp_BioGautama(61:end-60,61:end-60,:,:),1),2));
% pop_resp_BioGautama2 = squeeze(mean(mean(pop_resp_BioGautama2(61:end-60,61:end-60,:,:),1),2));
% pop_resp = squeeze(mean(mean(pop_resp(61:end-60,61:end-60,:,:),1),2));
% figure,plot_pop_response(pop_resp,0,0,param.prefVel)
% % figure,plot_pop_response(pop_resp_BioGautama,0,0,param.prefVel)
% figure,plot_pop_response(pop_resp_BioGautama2,0,0,param.prefVel)