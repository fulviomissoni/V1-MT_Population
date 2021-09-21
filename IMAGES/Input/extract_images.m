clear variables
close all
%%
slant=4:4:48;
tilt=30:30:360;
OldFolder=cd;
gaze=7;
cont=1;
for tt=1:length(tilt)
    for ss=1:length(slant)
        for kk=1:7
%             if cont>786
                file=['video_',num2str(cont),'_gaze_',num2str(kk),'_slant_',num2str(slant(ss)),'_tilt_',num2str(tilt(tt)),'.mat'];
                load(file)
                cd('images')
                if tmp>30
                    IIL = II_L(:,:,1:5:end);
                    IIR = II_R(:,:,1:5:end);
                else
                    if tmp<3
                        II_L = cat(3,II_L,II_L,II_L);
                        II_R = cat(3,II_R,II_R,II_R);
                    end
                    IIL = II_L(:,:,end-2:end);
                    IIR = II_R(:,:,end-2:end);
                end
                save(file,'IIL','IIR','tmp')
                cd(OldFolder)
                cont=cont+1;
                fprintf('%d\n',cont)
                tmp_ALL(cont)=length(IIL(1,1,:));
%             else
%                 cont=cont+1;
%                 fprintf('%d\n',cont)
%             end
        end
    end
end
            