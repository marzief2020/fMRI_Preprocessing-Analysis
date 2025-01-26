function PreProc_coregisteration(ref_img, source_img, other_img, i_est, i_reslice)

if i_est
    coreg_estimate = struct;
    % Ref

    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.ref = ref_img;

    % Source
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.source = source_img;

    % Other: ROIs
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.other = other_img;

    % Eoptions
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    % Run
    cfg_util('run',coreg_estimate.matlabbatch);
end

if i_reslice
    coreg_reslice = struct;
    % Ref
    coreg_reslice.matlabbatch{1}.spm.spatial.coreg.write.ref = ref_img;

    % Source
    coreg_reslice.matlabbatch{1}.spm.spatial.coreg.write.source = source_img;

    % Roptions
    coreg_reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    coreg_reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    coreg_reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    coreg_reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

    % Run
    cfg_util('run',coreg_reslice.matlabbatch);
end

