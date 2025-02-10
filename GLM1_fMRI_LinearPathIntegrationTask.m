% 1st level GLM

close all; clear all; clc;

%% -------------------- subjects ID/ directories/ params --------------------
iDropout = 0; iEstimate=1; iContrast=1;

nSub=30; for iSub=1:nSub, subj_name {iSub} = ['sub-',sprintf('%02d', iSub)]; end

home_dir = 'F:\Projects\PhD\SFB\Project_GPS\Data\Final';
% base_dir = fullfile (home_dir, 'analyzed_data\'); addpath (base_dir); cd(base_dir); 
base_dir = fullfile (home_dir, 'analyzed_data'); addpath (base_dir); cd(base_dir); 

code_dir = 'D:\Projects\PhD\SFB\Project_GPS\Codes\PreProc'; addpath(code_dir)
% addpath(fileparts(mfilename('fullpath')))

spm_dir = fileparts(which('spm')); addpath(spm_dir)

addpath('C:\Toolboxes\MatlabToolboxes\spm12\toolbox')
addpath ('C:\Toolboxes\MatlabToolboxes\spm12\toolbox\jubrain-anatomy-toolbox')

%% params
seq='ep3d'; ses='ses-1';
Runs={'run-1', 'run-2', 'run-3', 'run-4'};
echo_types={'echo-1', 'echo-2', 'multi_echo'}; iEcho=1;
TR=1.5; voxel_size=2.5; fwhm=2*voxel_size;
nDummy=3; nEchos=2; nRuns=4; nTasks=4; nSpeed=4; speed_val = 8.5.*[1, 2 , 3, 4];
pca_comp=5; iFIR=0;         
extra_explanation = '_';


condition_names = {'Mov_RT', 'Mov_RT_pmod', 'Mov_RT_pmod_stick', ...
    'Mov_RT_noCatch', 'Mov_RT_noCatch_pmod', 'Mov_RT_noCatch_pmod_stick', ...
    'Mov_RT_4SpeedReg', 'Mov_RT_4SpeedReg_stick', 'Mov_RT_4SpeedReg_noCatch', 'Mov_RT_4SpeedReg_noCatch_stick', ...
    'Mov_T1_1.5s_T2_7.5s_duration_1.5s_RT_noCatch_pmod', 'CentDist_noCatch_pmod', 'Mov_RT_noCatch_pmod_random', ...
    'Mov_T1T2_RT_noCatch_pmod_9s_1.5s', 'CentDistAbs_noCatch_pmod'};
n_regs = [2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 4, 4, 6, 4];

reg_idx = 9; %if reg_idx==11, extra_explanation = '_9s_1.5s_'; end

condition_name=condition_names{reg_idx}; n_regress=n_regs(reg_idx); 

noise_type={'MP-ref_acompcorr5', 'mp6_acompcorr5' ,'MP_ref'}; iNoise=1;
noise_reg=[1+pca_comp, 1, 1, 6+pca_comp];
n_reg=([1:numel(Runs)]-1)*(n_regress + noise_reg(iNoise));

if iNoise==1, rp_file_filter = 'aCompCorr_MP_ref.*';
elseif iNoise==2, rp_file_filter = 'aCompCorr_mp6_.*';
elseif iNoise==3, rp_file_filter = 'MP_ref.*';
end


for sn=3:29
    % sn=1;
    disp(['----------- GLM for ', subj_name{sn}, '-----------']);
    %% directories and params
    data_dir = fullfile(base_dir, subj_name{sn}, ses);
    func_dir = fullfile(data_dir, 'func');
    anat_dir = fullfile(data_dir, 'anat');
    events_dir = fullfile(data_dir, 'events');
    ROI_dir = fullfile(data_dir, 'ROI_coreg2T1');

    task_set = {'task-LinearPI', 'task-SceneLocalizer', 'task-MotionLocalizer', 'rest'};%,'task-GridLocalizer', 'task-DelcodeGrid'};
    task = task_set{1}; task_dir=fullfile(func_dir,task);
    echo_dir = fullfile(task_dir, echo_types{iEcho});
    if iEcho==3, echo_dir = fullfile(echo_dir, 'w_tSNR'); end

    glm_dir=fullfile(echo_dir, ['GLM', '_', condition_name, extra_explanation , noise_type{iNoise}]);
    % if iFIR, glm_dir=fullfile(task_dir, echo_types{iEcho}, ['FIR', '_', glm_name]); end
    if ~exist(glm_dir), mkdir(glm_dir); end


    %% STEP 1 -- first-level GLM ---  model specification and estimation

    if iEstimate
        first_level_GLM=struct;

        event_filter=['sub.*', task,'.*events.*mat'];
        events_mat_file = spm_select('FPList', fullfile(events_dir, task, condition_name), event_filter);
        % events_mat_file=spm_select('FPList', fullfile(events_dir, task_set{1}, condition_name), filter);
        events_mat_file = mat2cell(events_mat_file, ones(size(events_mat_file,1),1)); %convert to cell array

        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.dir = {glm_dir};
        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1; % reference slice

        % model specification
        for iRun=1:numel(Runs) % for each run
            scan_filter=['^srsub.*nii']; FilePathListThisRun=[];
            [FilePathListThisRun] = spm_select('FPList',fullfile(echo_dir, Runs{iRun}), scan_filter);
            FilePathListThisRun = mat2cell(FilePathListThisRun, ones(size(FilePathListThisRun,1),1));%convert to cell arrays
            %scans
            first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).scans= FilePathListThisRun;

            % set to empty
            if ~iFIR
                first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond = struct('name', {}, 'onset', {}, 'duration', {} , 'pmod', {});%, 'tmod', {}, 'pmod', {});
                % multi cond files
                first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).multi = {events_mat_file{iRun}}; %cond_timing .mat file
            else
                first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(1).name = 'Motion';
                events_mat_FIR = load(events_mat_file_FIR{iRun});

                first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(1).onset = cell2mat(events_mat_FIR.onsets);
                first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(1).duration = 10;
                first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(1).tmod = 0;
                first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
                first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(1).orth = 1;
                first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).multi = {''}; %cond_timing .mat file
            end

            % movement files
            first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).regress = struct('name', {}, 'val', {});
            rp_file = spm_select('FPlistRec', echo_dir, [rp_file_filter, 'run-',num2str(iRun), '.*txt']);
            % rp_file = spm_select('FPlistRec', fullfile(task_dir, 'echo-1'), [rp_file_filter, 'run-',num2str(iRun), '.*txt']);
            rp_file = mat2cell(rp_file, ones(size(rp_file,1),1));
            first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).multi_reg = rp_file;

            if iFIR
                first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).multi_reg = {''};
            end

            % hpf
            first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).hpf = 128;
        end

        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        if ~iFIR
            first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        end
        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.5;

        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.mask ={fullfile(func_dir, 'brain_mask.nii')};
        first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'FAST'; %'AR(1)'
        if iFIR
            first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length = 20;
            first_level_GLM.matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order = 10;
        end

        % model estimation
        estName=fullfile(glm_dir,'SPM.mat'); % SPM.mat for estimation
        first_level_GLM.matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {estName};
        first_level_GLM.matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        first_level_GLM.matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

        % run model specification and estimation
        spm fmri
        cfg_util('run',first_level_GLM.matlabbatch);

        disp('Step 10 - done!');

    end
    %% STEP 2 -- Contrast calcualtion
    if iContrast
        disp('STEP 11 -- Contrasts');
        cd (glm_dir); load SPM;
        %             rm con*
        %             rm spmT*

        try SPM=rmfield(SPM,'xCon'); catch clc; end

        switch (condition_name)
            case {condition_names{1}, condition_names{4}}
                nCotrasts=2;
                con_all=zeros(nCotrasts,size(SPM.xX.X,2));
                cNames = {'motionVbaseline','responseVbaseline'};
                con_all(1,n_reg+[1]) = 1; % 'motion vs. baseline'
                con_all(2,n_reg+[2]) = 1; % 'response vs. baseline'

            case {condition_names{2}, condition_names{3}, condition_names{5}, condition_names{6}, condition_names{13}}
                nCotrasts=4;
                con_all=zeros(nCotrasts,size(SPM.xX.X,2));
                cNames = {'motionVbaseline','speed_pmod', 'responseVbaseline', 'speed_pmodVmotion'};
                con_all(1,n_reg+[1]) = 1; % 'motion/avg vs. baseline'
                con_all(1,n_reg+[2]) = 0; % 'motion/pmod vs. baseline'

                con_all(2,n_reg+[1]) = 0; % 'motion/avg vs. baseline'
                con_all(2,n_reg+[2]) = 1; % 'motion/pmod vs. baseline'

                con_all(3,n_reg+[3]) = 1; % 'response vs. baseline'

                con_all(4,n_reg+[1]) = -1; % 'motion/avg vs. baseline'
                con_all(4,n_reg+[2]) = 1; % 'motion/pmod vs. baseline'

            case {condition_names{7}, condition_names{8}, condition_names{9}, condition_names{10}}
                
                nCotrasts=9; cNames = {'speed_mod', 'speed_exp_mod', 'speed_quad_mod', 'speed_log_mod', ...
                    'responseVbaseline', 'S1', 'S2', 'S3', 'S4'};
                con_all=zeros(nCotrasts,size(SPM.xX.X,2));
                
                z_scored_speed = zscore(speed_val);
                z_scored_exp_speed = zscore(exp(speed_val));
                z_scored_quad_speed = zscore(speed_val.^2);
                z_scored_log_speed = zscore(log(speed_val));

                for i=1:nSpeed
                    con_all(1,n_reg+[i]) = z_scored_speed(i); % 'speed modulation'
                end

                for i=1:nSpeed
                    con_all(2,n_reg+[i]) = z_scored_exp_speed(i); % 'exp speed modulation'
                end

                for i=1:nSpeed
                    con_all(3,n_reg+[i]) = z_scored_quad_speed(i); % 'quad speed modulation''
                end

                for i=1:nSpeed
                    con_all(4,n_reg+[i]) = z_scored_log_speed(i); % 'log speed modulation''
                end

                con_all(5,n_reg+[5]) = 1; % 'response vs. baseline'

                for iCon=1:nSpeed
                    con_all(iCon+5, n_reg+iCon) = 1; % Si
                end

            case {condition_names{11}, condition_names{14}}
                nCotrasts=9;
                con_all=zeros(nCotrasts,size(SPM.xX.X,2));
                cNames = {'motion_T1Vbaseline','speed_pmod_T1Vbaseline', 'motion_T2Vbaseline','speed_pmod_T2Vbaseline', ...
                'motion_T2-T1', 'speed_mod_T2-T1',  'responseVbaseline', 'motion_T1-T2', 'speed_mod_T1-T2'};
                % T1
                con_all(1,n_reg+[1]) = 1; % 'motion/avg vs. baseline'
                con_all(2,n_reg+[2]) = 1; % 'motion/pmod vs. baseline'

                %T2
                con_all(3,n_reg+[3]) = 1; % 'motion/avg vs. baseline'
                con_all(4,n_reg+[4]) = 1; % 'motion/pmod vs. baseline'

                %T2-T1
                con_all(5,n_reg+[1]) = -1; % 
                con_all(5,n_reg+[3]) = 1; % 

                con_all(6,n_reg+[2]) = -1; % 
                con_all(6,n_reg+[4]) = 1; % 

                % 'response vs. baseline'
                con_all(7,n_reg+[5]) = 1; 

                %T1-T2
                con_all(8,n_reg+[1]) = 1; % 
                con_all(8,n_reg+[3]) = -1; % 

                con_all(9,n_reg+[2]) = 1; 
                con_all(9,n_reg+[4]) = -1; % 

            case {condition_names{12}, condition_names{15}}
                nCotrasts=3;
                con_all=zeros(nCotrasts,size(SPM.xX.X,2));
                cNames = {'CentDist','CentDist_pmod', 'responseVbaseline'};
                con_all(1,n_reg+[1]) = 1; % 'motion/avg vs. baseline'

                con_all(2,n_reg+[2]) = 1; % 'motion/pmod vs. baseline'

                con_all(3,n_reg+[3]) = 1; % 'response vs. baseline'
        end


        %_____t-test for phases
        for conI=1:numel(cNames)
            con = con_all(conI,:);
            SPM.xCon(conI)=spm_FcUtil('Set',sprintf('%s',cNames{conI}), 'T', 'c',con',SPM.xX.xKXs);
        end
        SPM=spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save SPM SPM;
        disp('Step 11 - done!');
    end

    if iDropout
        glm_mask=spm_read_vols(spm_vol(fullfile(glm_dir, 'mask.nii')));
        ERC_l_img=spm_select('FPList', ROI_dir, ['^r.*', subj_name{sn}, '.*', ses, '.*ROI_left_ERC.nii$']);
        ERC_l_mask=spm_read_vols(spm_vol(ERC_l_img));

        ERC_r_img=spm_select('FPList', ROI_dir, ['^r.*', subj_name{sn}, '.*', ses, '.*ROI_right_ERC.nii$']);
        ERC_r_mask=spm_read_vols(spm_vol(ERC_r_img));

        noise_type{iNoise}
        % unique(ERC_l_mask)
        ERC_l_vox=sum(ERC_l_mask, 'all')
        ERC_r_vox=sum(ERC_r_mask, 'all')

        ERC_l_drop_vox=sum(ERC_l_mask.*glm_mask, 'all')
        ERC_r_drop_vox=sum(ERC_r_mask.*glm_mask, 'all')

    end
end

