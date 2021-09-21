function RDS_simulation(param,choice)
%Disparity tuning analysis with Random Dot Stereogram (RDS)
%Parameters contains alla the parameters of the population and choice
%permits to analyze: all the population or a single cell
%%%%
TUN = 1; %select cell in the center of the frame
switch choice
    case '2D_disp'
        tuning_surfaces(param)
    case  '1D_disp'    
        tuning_curves(param)
end

function tuning_surfaces(param)

    %File_name of the simulation data
    file_name='RDS_SurfTuning.mat';

    %INPUT DATA
    % input_file='RDS_correlatedN_seed5000.mat';
    % load(input_file);
    
    %%% disparity descriptors analysis - tuning surfaces
    d=-40:40;
    N_seed=500;
    n_orient=1;
    ph_shift=size(param.phShift,2);
    e=zeros(2,n_orient,ph_shift,length(d),length(d),N_seed); 
    taps = param.samp;
    for g=1:N_seed
            fprintf('finito giro %d\n',g)
        for iHD=1:length(d)
            for iVD=1:length(d)
                %select input
                [I(:,:,1), I(:,:,2)] = myRDS(d(iVD),d(iHD),1,g,taps,taps);
                I = repmat(I,[1 1 1 11]);
                I = permute(I,[1 2 4 3]);
                [MT, EC21, EC22] = pop_flow_V1MT(I,param,TUN);
                e(1,:,:,iHD,iVD,g) = squeeze(EC21);
                e(2,:,:,iHD,iVD,g) = squeeze(EC22); 
                fprintf('%d %d\n',iHD,iVD);
                clear I
            end
        end
        %Save data in SIMULATIONS Directory
        path = ['SIMULATIONS/disparity-tuning'];
        OldFolder = cd;
        cd(path);
        save(file_name,'e','param')
        fprintf('finito giro %d\n',g);
        cd(OldFolder)
    end
    end

function tuning_curves(param)
    
    % disparity descriptors analysis - tuning curves
    d = -20:20;
    N_seed = 100;
    j=1;
    n_orient = 1;
    filename = 'disp_tuning_curve_norm';
    ph_shift = size(param.phShift,2);
    n_cell=7;
    taps = param.samp;
    e = zeros(n_cell,n_orient,ph_shift,length(d),N_seed);
    for g=1:N_seed
        for id=1:length(d)
            [I(:,:,1), I(:,:,2)] = myRDS(d(id),0,1,g,taps,taps);
            I = repmat(I,[1 1 1 11]);
            I = permute(I,[1 2 4 3]);
            [MT,EC21,EC22,C1] = pop_flow_V1MT(I,param,TUN);
            e(1,:,:,id,g) = [squeeze(MT)];
            e(2,:,:,id,g) = [squeeze(EC21)];
            e(3,:,:,id,g) = [squeeze(EC22)];
            e(4,:,:,id,g) = [squeeze(C1{1})];
            e(5,:,:,id,g) = [squeeze(C1{2})];
            e(6,:,:,id,g) = [squeeze(C1{3})];
            e(7,:,:,id,g) = [squeeze(C1{4})];
            fprintf('id = %d\n',id);
            clear I
        end
        j=j+1;
        if j==100
            %Save data in SIMULATIONS Directory
            path = ['SIMULATIONS/disparity-tuning'];
            OldFolder = cd;
            cd(path);
            save(filename,'e','param','d')
            j=1;
            cd(OldFolder);
        end
        fprintf('g = %d\n',g);
    end
    end
end