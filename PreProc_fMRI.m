% Preprocessing steps: 
% 1st level GLM
% Marzieh Fereidouni- 2024

close all; clear all; clc;

% set iCopy to 1 to copy the structural images to the task folder
% set steps 1 based on the needed processes

%% steps to be done
steps=zeros(1, 10);  %% steps that should be run
iFS=1;
if ~iFS, steps(1:5)=1; iCopy=1; else, steps(6:9)=1; iCopy=0; end

iCheckReg=0; tSNR_plot=0;
iSplit=1; iDummy=1; iEndScans=1;  iRest=1;  iMP_Plot=1; 

sn_1=1; sn_end=29;

%% -------------------- subjects ID/ directories/ params --------------------

nSub=30; for iSub=1:nSub, subj_name {iSub} = ['sub-',sprintf('%02d', iSub)]; end

home_dir = 'F:\Projects\PhD\SFB\Project_GPS\Data\Final';
% base_dir = fullfile (home_dir, 'analyzed_data\'); addpath (base_dir); cd(base_dir)
base_dir = fullfile (home_dir, 'analyzed_data'); addpath (base_dir); cd(base_dir)

noise_PCA_dir = fullfile (base_dir, 'PCA_signals', 'ses-2', 'DelcodeGridCell'); addpath (base_dir);
if ~exist(noise_PCA_dir), mkdir(noise_PCA_dir); end

code_dir = 'D:\Projects\PhD\SFB\Project_GPS\Codes\PreProc'; addpath(code_dir)
% addpath(fileparts(mfilename('fullpath')))

spm_dir = fileparts(which('spm')); addpath(spm_dir)

addpath('C:\Toolboxes\MatlabToolboxes\spm12\toolbox')
addpath ('C:\Toolboxes\MatlabToolboxes\spm12\toolbox\jubrain-anatomy-toolbox')

%% params
TR=1.5; voxel_size=2.5; fwhm=2*voxel_size;

seq='ep3d'; ses='ses-3';
Runs_GD={'run-1', 'run-2'}; nRuns=2;
Runs_GD_realign={'run-1', 'run-2'}; 
ref_run = 1; 
% Runs_DG={'run-1', 'run-2'};

task_set = {'task-GridLocalizer', 'task-DelcodeGridCell', 'rest'};%,'task-GridLocalizer', 'task-DelcodeGrid'};
task = task_set{2};

echoes={'echo-1', 'echo-2'}; echo_types={'echo-1', 'echo-2', 'multi_echo'};   
echo_idx=1;   echo_ref=2; csf_filter = ['_CSF_PCA-' echo_types{echo_idx} '_run-'];
filter_fnImages_physnoise = ['^rsub.*', task, '.*nii$'];

nDummy=3; nEchos=2; nRuns_GD=2; pca_comp=5;

noise_type={'MP-ref_acompcorr5', 'mp6_acompcorr5'}; iNoise=1;
noise_reg=[1+pca_comp, 6+pca_comp];
if iNoise==1, rp_file_filter = 'aCompCorr_MP_ref.*'; elseif iNoise==2, rp_file_filter = 'aCompCorr_mp6_.*'; end



%% --------------------  Preproc --------------------
for sn=1:29


    disp(['----------- Preproc for ', subj_name{sn}, '-----------']);

    %% directories and params
    data_dir = fullfile(base_dir, subj_name{sn}, ses);
    func_dir = fullfile(data_dir, 'func');
    anat_dir = fullfile(data_dir, 'anat');
    events_dir = fullfile(data_dir, 'events');

    
    for i_task=1:numel(task_set)
        task_dir_set{i_task} = fullfile(func_dir, task_set{i_task});
    end

    task_dir = task_dir_set{2};
    if ~exist(task_dir), mkdir(task_dir); end

    task_weight = task_set{2}; task_dir_weight = task_dir_set{2};

    ROI_dir = fullfile(task_dir, 'ROIs');
    if ~exist(ROI_dir), mkdir(ROI_dir); end

    echo_dir = fullfile(task_dir, echo_types{echo_idx});
    if echo_idx == 3, echo_dir = fullfile(task_dir, echo_types{echo_idx}, 'w_tSNR'); end
    if ~exist(echo_dir), mkdir(echo_dir); end

    if iCopy
        [structural_T1_fn] = spm_select('FPList',func_dir ,['^' subj_name{sn},'_', ses, '_acq-mprageND_T1w.nii$']);
        [structural_T2_fn] = spm_select('FPList',func_dir ,['^' subj_name{sn},'_', ses, '_acq-tse_T2w.nii$']);
        copyfile(structural_T1_fn, task_dir)
        copyfile(structural_T2_fn, task_dir)
    end

    [structural_T1_fn] = spm_select('FPList',task_dir ,['^' subj_name{sn},'_', ses, '_acq-mprageND_T1w.nii$']);
    [structural_T2_fn] = spm_select('FPList',task_dir ,['^' subj_name{sn},'_', ses, '_acq-tse_T2w.nii$']);


    %% -------------------- STEP 1 -- 4D to 3D, remove dummy and extra scans --------------------

    if steps(1), disp('Step 1 -- Split 4D func images to 3D, and remove dummy and extra scans');

        % computing the valid number of scans of each run of the Gridlocalizer task
        for iRun=1:numel(Runs_GD)
            events_mat_file =readtable(spm_select('FPList', fullfile(events_dir, task), ...
                ['^sub.*', 'task-gc.*', Runs_GD{iRun} '.*.txt$']));
            run_duration = round(table2array(events_mat_file(end, 2)) + table2array(events_mat_file(end, 3)));
            disp(['valid scans of run-', num2str(iRun), ' = ', num2str(round(run_duration/TR)+2)])
            nScans(iRun)=round(run_duration/TR)+2;
        end

        % split 4D -> 3D and remove the first 3 dummy scans for all the tasks
        % and also remove the last extra scans for the linear-PI task
        disp ('DelcodeGridCell task')
        for iEcho=1:nEchos
            for iRun=1:numel(Runs_GD)
                disp ([Runs_GD{iRun}, ' --- ' , echoes{iEcho}])

                func3D_file_path = fullfile(task_dir, echoes{iEcho}, Runs_GD{iRun});
                if ~exist(func3D_file_path), mkdir(func3D_file_path); end

                func4D_file=spm_select('FPlist', func_dir, ...
                    ['^sub.*',task, '.*', Runs_GD{iRun}, '.*', echoes{iEcho},'.*.nii$']);

                % 4D to 3D
                if iSplit, spm_file_split(func4D_file, func3D_file_path); disp (['4D -->> 3D']); end
                func3D_file=spm_select('FPlist', func3D_file_path, ['^sub.*', task, '.*nii$']);
                func3D_file_cell = mat2cell(func3D_file, ones(size(func3D_file,1),1));
                nScan_total=size(func3D_file,1);

                % remove dummy scans
                if iDummy, disp('remove dummy scans'); delete(func3D_file_cell{1:nDummy}); end

                % remove extra scans
                if iEndScans
                    if nScan_total-nDummy > nScans(iRun)
                        delete(func3D_file_cell{nScans(iRun)+nDummy+1:end}); disp('remove extra scans'); end
                end
                clear func3D_file func3D_file_cell func4D_file
            end
        end

        disp('Step 1 -- Done!');
    end

    %% --------------------  STEP 2 -- Realign (estimate and reslice) --------------------
    if steps(2), disp('Step 2 -- Realignment');

        % select all the func images
        func_Realign_cell = cell(nEchos, nRuns);

        for iEcho=1:nEchos
            % GL task: run-2 is put first as the reference
            for iRun=1:nRuns
                func_Realign = spm_select('FPlist', fullfile(task_dir, echoes{iEcho}, Runs_GD_realign{iRun}), ...
                    ['^sub.*', task, '.*nii$']);
                func_Realign_cell(iEcho, iRun) = {mat2cell(func_Realign, ones(size(func_Realign,1),1))};
            end
        end

        PreProc_realignment_estimation(func_Realign_cell, nEchos)


        % plotting FD and MP
        if iMP_Plot
            for iEcho=1%:nEchos
                filter_rp = ['rp_sub.*' echoes{iEcho}, '.*'];
                PreProc_MP_plot(task_dir, iEcho, filter_rp, ref_run)
            end
        end

        % % checking the realignment
        % % rfunc_echo=spm_select('FPListRec', echo_dir , ['^rsub.*.nii$']);
        % rfunc_echo=spm_select('FPListRec', func_dir, ['^rsub.*echo-1.*.nii$']);
        % % rfunc_echo=spm_select('FPListRec', func_dir, ['^rsub.*echo-2.*.nii$']);
        %
        % spm_check_registration(rfunc_echo)
        % %         spm_check_registration(rfunc_echo2)
        % %         spm_check_registration
        % set(gcf, 'Resize', 'on')
        disp('Step 2 - Done!');
    end

    
    %% --------------------  coregisteration: estiamtion / segmentation/ reslice --------------------

    % Ref
    % two options: 1) ref image: first volume of first run, 2) ref image: mean image

    % 1) ref image: first volume of the reference run
    func4D_file=spm_select('FPlist', fullfile(task_dir, echoes{echo_ref}, Runs_GD_realign{1}), ['^sub.*', task, '.*nii$']);
    func4D_file = mat2cell(func4D_file, ones(size(func4D_file,1),1));
    ref_img={func4D_file{1}};

    %     % 2) ref image: mean image
    %     ref_img={spm_select('FPlist', taskDir, ['^meansub.*bold.nii$'])};

    %% STEP 3 -- Coregister structural image to first dynamic image (estimate only)
    if steps(3),  disp('Step 3 -- Coregister structural image to first functional image');

        source_img = {structural_T1_fn};

        % ROI_Files=spm_select('FPList', ROI_dir, '^sub.*ROI.*nii$');
        % ROI_Files=deblank(mat2cell(ROI_Files, ones(size(ROI_Files,1),1)));
        % other_img = ROI_Files;
        other_img= {''};

        i_est = 1; i_reslice = 0;
        PreProc_coregisteration(ref_img, source_img, other_img, i_est, i_reslice)

        disp('Step 3 - Done!');
    end

    %% STEP 4 -- Segmentation of coregistered structural image into GM, WM, CSF, etc
    % (with implicit warping to MNI space, saving forward and inverse transformations)
    if steps(4), disp('Step 4 -- Segmentation');

        PreProc_segmentation(structural_T1_fn)

        disp('Step 4 - done!');
    end

    %% STEP 5 -- Reslice all to functional-resolution image grid
    if steps(5),  disp('Step 5 -- Reslice all to functional-resolution image grid');

        % Source
        % ROI_Files=spm_select('FPList', ROI_dir, '.*ROI.*nii$');
        % ROI_Files=deblank(mat2cell(ROI_Files, ones(size(ROI_Files,1),1)));

        source_img = {};
        segmented_structural_fn=spm_select('FPList', task_dir, '^c.*T1w.nii$');
        source_img=[mat2cell(segmented_structural_fn, ones(size(segmented_structural_fn,1),1)); structural_T1_fn];
        % source_img=[mat2cell(segmented_structural_fn, ones(size(segmented_structural_fn,1),1));...
        %     structural_T1_fn; structural_T2_fn; ...
        %     deblank(ROI_Files)];

        other_img = {''};

        i_est = 0; i_reslice = 1;
        PreProc_coregisteration(ref_img, source_img, other_img, i_est, i_reslice)

        disp('Step 5 - done!');
    end

    %% -------------------- Step 6: Generating masks --------------------

    if steps(6), disp('Step 6 -- Generating brain mask');

        % first the the Brain (w/o skulp) mean image should be generated with FS

        ref_mask=[ref_img{:} ',1'];
        EPI_img_Brain = spm_select('FPList', task_dir, '^rsub.*mprageND_T1w.*Brain.*nii'); % generated by ANTs

        PreProc_MaskGeneration(structural_T1_fn, EPI_img_Brain, ref_mask, task_dir, ses) % ref_mask is used to generate mask images headers

        disp('Step 6 - done!');
    end

    %% -------------------- Step 7 -- Echo combination  --------------------
    if steps(7)
        disp('Step 7 -- Combination of echoes');

        TE =[12.6, 33.5];

        brain_mask = spm_select('FPList', task_dir,  '.*brain_mask_dil.nii');

        ref_scans = cell(nEchos,1);

        % GridLocalizer task
        func_scans = cell(nEchos, nRuns);

        for iEcho=1:nEchos
            for iRun=1:nRuns
                func_fn = spm_select('FPlist', fullfile(task_dir, echoes{iEcho}, Runs_GD{iRun}), ['^rsub.*', task, '.*nii$']);
                func_scans{iEcho, iRun} = mat2cell(func_fn, ones(size(func_fn,1),1));
                if iRun==1
                    ref_scans{iEcho}=func_scans{iEcho, iRun}(1:80);
                end
            end
        end

        iRest=1; iComb=2; 
        PreProc_CombineEchoes(iRest, iComb, subj_name{sn}, ses, ref_scans, task_dir, func_scans, brain_mask, TE)

        clear func_scans

        disp('Step 7 - done!');
    end

    if iCheckReg

        iEcho=3;
        echo_dir_checkreg = fullfile(task_dir, echo_types{iEcho});
        if iEcho == 3, echo_dir_checkreg = fullfile(task_dir, echo_types{iEcho}, 'w_tSNR'); end

        rfunc_combine = spm_select('FPListRec',echo_dir_checkreg, ['^rsub.*.nii$']);
        spm_check_registration(rfunc_combine)
        set(gcf, 'Resize', 'on')
    end

    %% STEP 8: Physiological noise
    if steps(8)
        disp('Step 8 -- computing physiological noise');

        noise_mask = double(spm_read_vols(spm_vol(spm_select('FPList', task_dir,  '.*CSF_mask_erode.nii'))));
        noise_mask_Nan = noise_mask; noise_mask_Nan((noise_mask_Nan==0))=NaN;  % removing non-mask area from the pca computation

        GM_mask = double(spm_read_vols(spm_vol(spm_select('FPList', task_dir,  '.*GM_mask.nii'))));
        GM_mask_Nan = GM_mask;  GM_mask_Nan((GM_mask_Nan==0))=NaN;  % removing non-mask area from the pca computation

        rp_filter='^MP_ref.*';

        func_scans = cell(1,nRuns);

        % generating rp files for the multi-echo scans from echo-1 and echo-2
        rp_val={}; e_rp={}; d_rp={};

        for iRun=1:nRuns
            rp_file = spm_select('FPListRec', fullfile(task_dir, echo_types{1}), [rp_filter 'run-' num2str(iRun)]);
            [d_rp{iRun}, e_rp{iRun}, f_rp] = fileparts(rp_file);
            rp_val{iRun} = load(rp_file);
        end

        for iRun=1:nRuns
            func_fn = spm_select('FPlist', fullfile(echo_dir, Runs_GD{iRun}), filter_fnImages_physnoise);
            func_scans{1,iRun} = mat2cell(func_fn, ones(size(func_fn,1),1));
        end


        % figure;
        corr_gm_noise=0; corr_rp_noise=0; corr_rp_gm=0;

        for iRun=1:nRuns
            disp(['---Run-', num2str(iRun)])

            func_scan_cell = struct;
            func_scan_cell.info = spm_vol(cell2mat(func_scans{1,iRun}));

            temp_info = func_scan_cell.info(1); [d_fn, e_fn, f_fn] = fileparts(temp_info.fname);

            func_scan_cell.vol =  spm_read_vols(func_scan_cell.info);
            func_scan_cell.vol_noise_masked = noise_mask_Nan.*func_scan_cell.vol;
            func_scan_cell.vol_GM_masked = GM_mask_Nan.*func_scan_cell.vol;

            t_signal=size(func_scan_cell.info,1);
            func_scan_cell.noise_siganl = reshape(func_scan_cell.vol_noise_masked(~isnan(func_scan_cell.vol_noise_masked)), [], t_signal);
            func_scan_cell.GM_siganl = reshape(func_scan_cell.vol_GM_masked(~isnan(func_scan_cell.vol_GM_masked)), [], t_signal);

            [noise_coeff,noise_score,noise_latent, noise_TSQUARED, noise_EXPLAINED] = pca(func_scan_cell.noise_siganl);
            [GM_coeff,GM_score, GM_latent, GM_TSQUARED, GM_EXPLAINED] = pca(func_scan_cell.GM_siganl);

            noise_coeff_cumsum = cumsum(noise_EXPLAINED(1:5))'
            % GM_coeff_cumsum = cumsum(GM_EXPLAINED);

            nuisance_reg  = [rp_val{iRun}, noise_coeff(:, 1:5)];
            nuisance_reg_fname=fullfile(echo_dir, Runs_GD{iRun}, ['aCompCorr_' e_rp{iRun} '.txt']);
            save(nuisance_reg_fname,'nuisance_reg','-ascii');

            % noise plots
            figure;
            for iPlot=1:15
                subplot(15,1,iPlot)
                plot(noise_coeff(:,iPlot));
            end
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
            saveas (gcf, fullfile (noise_PCA_dir, [subj_name{sn} csf_filter num2str(iRun) '.png'])) ;

            % GM plots
            % figure;
            % for iPlot=1:15
            %     subplot(15,1,iPlot)
            %     plot(GM_coeff(:,iPlot));
            % end
            % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
            % saveas (gcf, fullfile (noise_PCA_dir, [subj_name{sn} '_GM_PCA_multi-echo_run-' num2str(iRun) '.png'])) ;
        end
        close all

        disp('Step 8 - done!');
    end

    %% STEP 9 -- Gaussian kernel smoothing of realigned data
    if steps(9)

        disp('STEP 9 -- Gaussian kernel smoothing of realigned data');

        for iEcho=[1, 3]

            disp(['smoothing - ', echo_types{iEcho}])
            echo_dir_smooth = fullfile(task_dir, echo_types{iEcho});
            if iEcho == 3, echo_dir_smooth = fullfile(task_dir, echo_types{iEcho}, 'w_tSNR'); end

            rfunc_fn_cell={};
    
            rfunc_fn=spm_select('FPlistRec', echo_dir, '^rsub.*.nii$');
            rfunc_fn_cell = deblank(mat2cell(rfunc_fn, ones(size(rfunc_fn,1),1)));
    
            smooth = struct;
            smooth.matlabbatch{1}.spm.spatial.smooth.data =rfunc_fn_cell;% vertcat(rfunc_fn_cell{:});
            smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
            smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            smooth.matlabbatch{1}.spm.spatial.smooth.im = 0;
            smooth.matlabbatch{1}.spm.spatial.smooth.prefix = 's';
            spm fmri
            cfg_util('run',smooth.matlabbatch);
        end

        disp('Step 9 - done!');

    end

    %% STEP 10 -- computing sSNR and tSNR
    if steps(10)

        % ROIs should be co-registered to the mean_EPI space first in ANTs

        disp('STEP 10 -- Calculating sSNR and tSNR');

        % metric_filename = fullfile(echo_dir, [subj_name{sn}, '_', ses, '_acq-', seq, ['task-', task],  '_metrics.txt']);
        % disp(' ');
        % disp(['Metric file: ' metric_filename]);
        % fp_Metrics = fopen(metric_filename, 'w+');
        % fprintf(fp_Metrics, 'Participant = %s; ses = %s; task = %s' , subj_name{sn}, ses, task);


        % for iEcho=1:nEchos
        % Calculate sSNR
        cd (task_dir)

        tSNR_temp = load(fullfile(base_dir, 'tSNR_maps', ses, task, 'tSNR.mat'));
        nVoxel_temp = load(fullfile(base_dir, 'tSNR_maps', ses, task, 'Num_voxel.mat'));

        brain_mask=spm_read_vols(spm_vol(spm_select('FPList', task_dir,  '.*brain_mask.nii')));
        brain_mask_Nan = brain_mask;  brain_mask_Nan((brain_mask_Nan==0))=NaN;  % removing non-mask area from the pca computation

        GM_mask = double(spm_read_vols(spm_vol(spm_select('FPList', task_dir,  '.*GM_mask.nii'))));
        GM_mask_Nan = GM_mask;  GM_mask_Nan((GM_mask_Nan==0))=NaN;  % removing non-mask area from the pca computation

        mean_img_info = spm_vol(spm_select('FPlist', fullfile(echo_dir, Runs_GD{2}), '^rsub.*00001.*.*nii$'));
        mean_img = spm_read_vols(mean_img_info);
        mean_Img_masked = brain_mask_Nan.* mean_img;

        sSNR(sn) = mean(mean_Img_masked, [1,2,3], 'omitnan')./std(mean_Img_masked,[],[1,2,3], 'omitnan');
        display([subj_name{sn}, ' - sSNR_brain_task-', task,' = ',num2str(sSNR)]);
        sSNR_temp = load(fullfile(base_dir, 'tSNR_maps', ses, task, 'sSNR.mat'));
        % sSNR_temp.sSNR(sn) = sSNR(sn);
        sSNR_temp.sSNR(end+1) = sSNR(sn);
        save(fullfile(base_dir, 'tSNR_maps', ses, task, 'sSNR.mat') , 'sSNR' );
        % fprintf(fp_Metrics, '\n\nsSNR_WholeBrain = %2.2f\n', sSNR);
        % end

        % Calculate tSNR
        rfunctional4D_fn=spm_select('FPListRec', echo_dir, ['^rsub.*', task, '.*nii$']);

        rfunctioanl4D_Run_info = spm_vol(rfunctional4D_fn);
        rfunctioanl4D_Run_img = spm_read_vols(rfunctioanl4D_Run_info);
        rfunctioanl4D_img_brain_masked = rfunctioanl4D_Run_img.*brain_mask_Nan;
        tSNR_brain = mean(rfunctioanl4D_img_brain_masked,4, 'omitnan')./std(rfunctioanl4D_img_brain_masked,[],4, 'omitnan');
        display(['tSNR_brain_task' task ,'=',num2str(mean(tSNR_brain,[1, 2, 3], 'omitnan'))]);

        rfunctioanl4D_img_GM_masked=rfunctioanl4D_Run_img.*GM_mask_Nan;
        tSNR_GM=mean(rfunctioanl4D_img_GM_masked,4)./std(rfunctioanl4D_img_GM_masked,[],4);
        display(['tSNR_GM_task' task ,'=',num2str(mean(tSNR_GM,[1, 2, 3], 'omitnan'))]);

        % save tSNR_brain map
        tSNR_brain_info=mean_img_info;
        tSNR_brain_info.fname=[echo_dir filesep [subj_name{sn}, '_', ses, '_tSNR_brain_task-' task, '.nii']];
        spm_write_vol(tSNR_brain_info, tSNR_brain);

        tSNR_brain_info.fname=[echo_dir filesep [subj_name{sn}, '_', ses, '_tSNR_GM_task-' task, '.nii']];
        spm_write_vol(tSNR_brain_info, tSNR_GM);
        %
        % fprintf(fp_Metrics, '\ntSNR_WholeBrain = %2.2f\n', mean(tSNR_brain,[1, 2, 3], 'omitnan'));
        % fprintf(fp_Metrics, '\ntSNR_GM = %2.2f\n', mean(tSNR_GM,[1, 2, 3], 'omitnan'));

        % tSNR for ROIs
        ROI_File = spm_select('FPList', ROI_dir, '^coreg.*nii$');
        ROI_File = deblank(mat2cell(ROI_File, ones(size(ROI_File,1),1)));

        for iROI=1 : numel(ROI_File)
            fROI=ROI_File{iROI};
            fROI_info = spm_vol(fROI);
            [d_fROI, e_fROI, r_fROI]=fileparts(fROI_info.fname);

            fROI_img = (spm_read_vols(fROI_info));
            fROI_img=round(fROI_img);

            fROI_img_Nan = fROI_img;
            fROI_img_Nan((fROI_img_Nan==0))=NaN;

            e_fROI2=e_fROI(end-11:end);
            display(['Num of voxels in ', e_fROI2,' = ',num2str(sum((fROI_img>0),[1 2 3]))]);

            fROI_img_4D = fROI_img_Nan.*(rfunctioanl4D_img_brain_masked);
            tSNR_fROI = mean(fROI_img_4D,4, 'omitnan')./std(fROI_img_4D,[],4, 'omitnan');
            tSNR_temp.tSNR_all{end+1, iROI} = mean(tSNR_fROI,[1, 2, 3], 'omitnan');
            nVoxel_temp.nVoxel{end+1, iROI} =  sum((fROI_img>0),[1 2 3]);
            % tSNR_temp.tSNR_all{sn, iROI} = mean(tSNR_fROI,[1, 2, 3], 'omitnan');
            % nVoxel_temp.nVoxel{sn, iROI} =  sum((fROI_img>0),[1 2 3]);

            fROI_Img_thr = fROI_img(tSNR_fROI>10);
            % nVoxel_temp.nVoxel_thr{sn, iROI} =  sum((fROI_Img_thr>0),[1 2 3]);
            nVoxel_temp.nVoxel_thr{end+1, iROI} =  sum((fROI_Img_thr>0),[1 2 3]);
            
            % fprintf(fp_Metrics, '\n%s\nnumber of voxels = %d\ntSNR = %.2f\n', ...
        end

        % tSNR_temp.tSNR_all{sn, iROI+1} =  mean(tSNR_GM,[1, 2, 3], 'omitnan');
        tSNR_temp.tSNR_all{end+1, iROI+1} =  mean(tSNR_GM,[1, 2, 3], 'omitnan');


        % fclose(fp_Metrics);
        save(fullfile(base_dir, 'tSNR_maps', ses, task, 'sSNR.mat') , 'sSNR_temp' );
        tSNR_all = tSNR_temp.tSNR_all;
        save(fullfile(base_dir, 'tSNR_maps', ses, task, 'tSNR.mat') , 'tSNR_all' );
        save(fullfile(base_dir, 'tSNR_maps', ses, task, 'Num_voxel.mat') , 'nVoxel_temp' );

        disp('Step 10 - done!');
        close all
    end

end

if tSNR_plot
    tSNR_temp = load(fullfile(base_dir, 'tSNR_maps', ses, task, 'tSNR.mat'));
    nVoxel_temp = load(fullfile(base_dir, 'tSNR_maps', ses, task, 'Num_voxel'));

    ROI_names = {'left_ERC', 'left_HPC', 'left_PHC', 'right_ERC', 'right_HPC', 'right_PHC', 'GM'};

    figure();
    h = boxplot(cell2mat(tSNR_temp.tSNR_all));

    % Customize the x-axis labels
    set(gca, 'XTickLabel', ROI_names);
    xtickangle(45); % Rotate labels for better readability if needed

    % Create a "fake" legend
    hold on;
    for i = 1:length(ROI_names)
        plot(nan, nan, 'square', 'DisplayName', ROI_names{i});
    end
    hold off;

    title('tSNR for 29 participants - Grid Localizer task')
end
%
