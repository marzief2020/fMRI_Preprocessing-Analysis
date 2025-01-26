function PreProc_CombineEchoes(iRest, iComb, subj_name, ses, ref_scans, task_dir, func_scans, brain_mask_dil, TE)

nEchos = size(TE, 2); nRuns = size(func_scans, 2); n_beforeEcho = 17; mFactor=10^4;
brain_mask_vol = spm_read_vols(spm_vol(brain_mask_dil));
multi_echo_dir = fullfile(task_dir, 'multi_echo');  if ~exist(multi_echo_dir), mkdir(multi_echo_dir); end
comb_type = {'tSNR', 'w_tSNR'};

%% generating weight images from the rest echo scans

if iRest
    % sub-directiory for tSNR or w_tSNR
    multi_echo_sub_dir = fullfile(task_dir, 'multi_echo', comb_type{iComb});  if ~exist(multi_echo_sub_dir), mkdir(multi_echo_sub_dir); end

    rest_scan = struct; norm_image = 0;
    for iEcho = 1:nEchos
        % reading the scan volumes and mask them
        rest_scan(iEcho).info = spm_vol(ref_scans{iEcho}); temp_info = rest_scan(iEcho).info{1};
        temp_vol = cellfun(@spm_read_vols, rest_scan(iEcho).info, 'UniformOutput', false);
        rest_scan(iEcho).vol = cellfun(@(x) bsxfun(@times, x, brain_mask_vol), temp_vol, 'UniformOutput', false);

        % Convert the cell array to a 4D matrix
        rest_scan_4D = cat(4, rest_scan(iEcho).vol{:});

        % Calculate the tSNR
        rest_scan(iEcho).tSNR = mean(rest_scan_4D, 4, 'omitnan')./std(rest_scan_4D, 0, 4, 'omitnan');

        % wrting the tSNR images
        tsnr_info = temp_info;
        tsnr_info.fname=fullfile(multi_echo_sub_dir,  [subj_name, '_', ses, '_', comb_type{iComb},'_rest_echo-' num2str(iEcho), '.nii']);
        if iComb==1, output_scan = rest_scan(iEcho).tSNR; end%output_scan(isnan(output_scan)) = 0;
        if iComb==2, output_scan = TE(iEcho).*rest_scan(iEcho).tSNR; end
        spm_write_vol(tsnr_info, output_scan);

        % generating the norm image
        norm_image =  norm_image + output_scan; clear output_scan rest_scan_4D
    end

    % writing the norm image
    norm_image(norm_image<0)=0;
    norm_info = temp_info; norm_info.fname=fullfile(multi_echo_sub_dir,  [subj_name, '_', ses, '_norm_', comb_type{iComb},'_rest.nii']);
    spm_write_vol(norm_info, norm_image); clear norm_info

    for iEcho = 1:nEchos
        if iComb==1, temp_tsnr_mat = rest_scan(iEcho).tSNR; end
        if iComb==2, temp_tsnr_mat = TE(iEcho).*rest_scan(iEcho).tSNR; end

        temp_tsnr_mat(temp_tsnr_mat<0)=0;
        rest_scan(iEcho).wImg = (mFactor).*round((temp_tsnr_mat ./ norm_image), 4);
        output_scan = rest_scan(iEcho).wImg; output_scan(isnan(output_scan)) = 0;
        weight_info = temp_info; weight_info.fname=fullfile(multi_echo_sub_dir,  [subj_name, '_', ses, '_', 'weight_', comb_type{iComb}, '_echo-' num2str(iEcho), '.nii']);
        spm_write_vol(weight_info, output_scan); clear weight_info output_scan temp_tsnr_mat
    end
end
%% combine the echoes of the main runs

% w1=TE(1)/(TE(1)+TE(2)); w2=TE(2)/(TE(1)+TE(2));

    task_scan = struct; 
    multi_echo_dir = fullfile(task_dir, 'multi_echo', comb_type{iComb});

    if ~iRest
        rest_scan = struct;
        for iEcho = 1:nEchos
            weigth_filter = ['.*weight.*', comb_type{iComb}, '.*echo-' num2str(iEcho) '.*'];
            rest_scan(iEcho).wImg = spm_read_vols(spm_vol(spm_select('FPList', multi_echo_dir, weigth_filter)));
        end
    end

    disp(['Starting echo combination for ', comb_type{iComb}])

    for iRun = 1:nRuns
        disp(['--------- Run-', num2str(iRun), ' ---------'])

        multi_echo_run_dir = fullfile(multi_echo_dir, ['run-', num2str(iRun)]);
        if ~exist(multi_echo_run_dir), mkdir(multi_echo_run_dir); end

        % read the scans of each echo in the run and multiply them with the weighted tSNR rest scan
        for iEcho = 1:nEchos
            task_scan(iRun).Echo(iEcho).info = spm_vol(func_scans{iEcho, iRun}); scan_info = task_scan(iRun).Echo(1).info{1};
            task_scan(iRun).Echo(iEcho).vol = cellfun(@spm_read_vols, task_scan(iRun).Echo(iEcho).info, 'UniformOutput', false);
            task_scan(iRun).Echo(iEcho).w_vol = cellfun(@(x) bsxfun(@times, x, rest_scan(iEcho).wImg./(mFactor)), task_scan(iRun).Echo(iEcho).vol, 'UniformOutput', false);
        end
        disp('echos are weighted!')

        % combine the weighted scans from each echo
        nVols = size(task_scan(iRun).Echo(iEcho).info,1); % number of volumes in the run
        w_vol_merge = cat(2, task_scan(iRun).Echo.w_vol);
        task_scan(iRun).combined_vol = arrayfun(@(idx) sum(cat(4, w_vol_merge{idx, :}), 4), 1:nVols, 'UniformOutput', false)';
        disp('weighted echos are summed!')

        % Generate the headerCell using cellfun
        [d_scan, e_scan, f_scan] = fileparts(scan_info.fname);
        task_scan(iRun).headerCell =  task_scan(iRun).Echo(1).info;
        task_scan(iRun).headernamesCell = cellfun(@(x) strcat(multi_echo_run_dir,filesep,[e_scan(1:end-n_beforeEcho),'bold_',sprintf('%05d', x), '.nii']), num2cell(1:nVols), 'UniformOutput', false)';
        task_scan(iRun).headerCell = cellfun(@(hdr, name) setfield(hdr, 'fname', name), task_scan(iRun).headerCell, task_scan(iRun).headernamesCell, 'UniformOutput', false);
        disp('combined-echo names are generated!')

        % write the combined scans 
        cellfun(@(hdr, img) spm_write_vol(hdr, img),  task_scan(iRun).headerCell,  task_scan(iRun).combined_vol);
        disp('combined echos are generated!')
    end
end

           