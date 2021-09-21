figure,
for t=2:21
    Atmp = imtranslate(A(:,:,t-1),[0.5,0],'bicubic','OutputView','full');
    dim = size(Atmp,2)-size(A,2);
    Atmp(:,dim) = Atmp(:,end-dim);
    Atmp = Atmp(:,1:end-dim);
    A(:,:,t) = Atmp;
    Atmp(Atmp>0.5)=1;
    Atmp(Atmp<=0.5)=0;
    subplot(2,1,1),imagesc(Atmp)
    colorbar
   	subplot(2,1,2),imagesc(A(:,:,t))
    colorbar
    pause
end
