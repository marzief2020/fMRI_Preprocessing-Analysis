%-----------------------------------------------------------------------
% Second-Level Analysis: Model Specifacation and Estimation

% Hint: the format of folders should be like BIDS otherwise you need to
% make some changes in the directory names.

% Marzieh Fereidouni_ 2024
%-----------------------------------------------------------------------

close all; clear all; clc;

%% -------------------- subjects ID/ directories/ params --------------------

iSub=1; for sub_idx=[2:6, 8:19, 21:29], subj_name{iSub} = ['sub-',sprintf('%02d', sub_idx)]; iSub=iSub+1; end

home_dir = 'F:\Projects\PhD\SFB\Project_GPS\Data\Final';
% base_dir = fullfile (home_dir, 'analyzed_data\'); addpath (base_dir); cd(base_dir)
base_dir = fullfile (home_dir, 'analyzed_data'); addpath (base_dir); cd(base_dir)

code_dir = 'D:\Projects\PhD\SFB\Project_GPS\Codes\PreProc'; addpath(code_dir)
% addpath(fileparts(mfilename('fullpath')))

spm_dir = fileparts(which('spm')); addpath(spm_dir);     fs=filesep;

addpath('C:\Toolboxes\MatlabToolboxes\spm12\toolbox')
addpath ('C:\Toolboxes\MatlabToolboxes\spm12\toolbox\jubrain-anatomy-toolbox')

echo_dir = 'ses-1/func/task-LinearPI/echo-1'; echo_type = 'GLM2_LinearPI_Echo1';

% Set the task that we want to preprocess its related images to 1
cNames = {'motionVbaseline','speed_pmod', 'responseVbaseline'};
% cNames = {'speedVbaseline', 'responseVbaseline'};

% Set the names of subjects and runs for each task, here
subjects = [subj_name];
task_name='LinearPI';
extra_explanation = '_';

condition_names = {'Mov_RT', 'Mov_RT_pmod', 'Mov_RT_pmod_stick', ...
    'Mov_RT_noCatch', 'Mov_RT_noCatch_pmod', 'Mov_RT_noCatch_pmod_stick', ...
    'Mov_RT_4SpeedReg', 'Mov_RT_4SpeedReg_stick', 'Mov_RT_4SpeedReg_noCatch', 'Mov_RT_4SpeedReg_noCatch_stick', ...
    'Mov_T1_1.5s_T2_7.5s_duration_1.5s_RT_noCatch_pmod', 'CentDist_noCatch_pmod', 'Mov_RT_noCatch_pmod_random', ...
    'Mov_T1T2_RT_noCatch_pmod_9s_1.5s', 'CentDistAbs_noCatch_pmod'};

reg_idx=9;
n_regs = [2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 4, 4, 6, 4];

noise_type={'MP-ref_acompcorr5', 'mp6_acompcorr5'}; iNoise=1;

% con_case = 'pmod_reg';
con_case = 'separate_reg';
% con_case = 'speed_mod_T1-T2';
% con_case = 'Cent_reg';

wcon_prefix='wcon_000';
% wcon_prefix='w_ants_con_000';

switch con_case
    case 'pmod_reg'
        % con_names = {'motionVbaseline', 'speed_pmod', 'responseVbaseline', 'speed_pmodVmotion'};
        con_names = {'motionVbaseline', 'speed_pmod', 'responseVbaseline'};
        for iCon=1:numel(con_names), con_image{iCon} = [wcon_prefix, num2str(iCon), '.nii']; end
        % condition_name=condition_names{5}; n_regress=n_regs(5); % 'Mov_RT_noCatch_pmod'
        % condition_name=condition_names{13}; n_regress=n_regs(13); % 'Mov_RT_noCatch_pmod_random'

    case 'separate_reg'
        con_names = {'speed_mod', 'speed_exp_mod','speed_quad_mod' ,'speed_log_mod' , 'responseVbaseline', 'S1', 'S2', 'S3', 'S4'};
        for iCon=1:numel(con_names), con_image{iCon} = [wcon_prefix, num2str(iCon), '.nii']; end
        % condition_name=condition_names{9}; n_regress=n_regs(9); % 'Mov_RT_4SpeedReg_noCatch'

    case 'speed_mod_T1-T2'
        con_names = {'motion_T1Vbaseline','speed_pmod_T1Vbaseline', 'motion_T2Vbaseline','speed_pmod_T2Vbaseline', ...
            'motion_T2-T1', 'speed_mod_T2-T1',  'responseVbaseline', 'motion_T1-T2', 'speed_mod_T1-T2'};
        for iCon=1:numel(con_names), con_image{iCon} = [wcon_prefix, num2str(iCon), '.nii']; end
        % condition_name = condition_names{11}; n_regress=n_regs(11); % 'Mov_RT_4SpeedReg_noCatch'
        % extra_explanation = '_9s_1.5s_';

    case 'Cent_reg'
        con_names = {'CentDist', 'CentDist_pmod', 'responseVbaseline'};
        % for iCon=1:numel(con_names), con_image{iCon} = ['wcon_000', num2str(iCon), '.nii']; end
        for iCon=1:numel(con_names), con_image{iCon} = [wcon_prefix, num2str(iCon), '.nii']; end
        % condition_name=condition_names{12}; n_regress=n_regs(12); % 'Mov_RT_noCatch_pmod'
        
end
condition_name=condition_names{reg_idx}; n_regress=n_regs(reg_idx); 

glm_dir = ['GLM_', condition_name, extra_explanation, noise_type{iNoise}];
% glm2_dir_all = fullfile(base_dir, 'GLM_2nd\ANTs\ses-1\task-LinearPI\sub_26\');

glm2_dir_all = fullfile(base_dir, echo_type);
if ~exist(glm2_dir_all), mkdir(glm2_dir_all); end


nCon = numel(con_names); if strcmp(con_case, 'separate_reg'), nCon=5; end

for iCon=1:nCon

    glm2_dir= fullfile(glm2_dir_all, ['GLM2_', condition_name, extra_explanation, noise_type{iNoise}], con_names{iCon});
    if ~exist(glm2_dir), mkdir(glm2_dir); end

    %% Second Level Analysis

    % con_scans = cellfun(@(c) [base_dir fs c fs ['ses-1/func/task-LinearPI/multi_echo/w_tSNR' fs glm_dir fs con_image{iCon}]], subjects','uni',false);
    con_scans = cellfun(@(c) [base_dir fs c fs [echo_dir fs glm_dir fs con_image{iCon}]], subjects','uni',false);
    glm_2nd = struct;
    glm_2nd.matlabbatch{1}.spm.stats.factorial_design.dir = {glm2_dir};
    glm_2nd.matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = con_scans;
    glm_2nd.matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    glm_2nd.matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    glm_2nd.matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    glm_2nd.matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    glm_2nd.matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    glm_2nd.matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    glm_2nd.matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    glm_2nd.matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    glm_2nd.matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
    glm_2nd.matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
    glm_2nd.matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
    glm_2nd.matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
    glm_2nd.matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
    glm_2nd.matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
    glm_2nd.matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
    glm_2nd.matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
    glm_2nd.matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
    glm_2nd.matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    spm('defaults','fmri');
    spm_figure('GetWin','Graphics');
    spm_jobman('initcfg');
    spm_jobman('run', glm_2nd.matlabbatch);

end
