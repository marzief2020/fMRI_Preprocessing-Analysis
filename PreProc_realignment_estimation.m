%% Realignment

function PreProc_realignment_estimation(func_Realign_cell, nEchos)  % funce_realign is the scans from echo-2

% estimation for echo-1
spm fmri
data_fn = func_Realign_cell(1,:);
realign_estimate = struct;
realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.data = data_fn;
realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 1;  % reference to the mean
realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
cfg_util('run',realign_estimate.matlabbatch);


% apply the estimation to echo-2
for iEcho=2:nEchos
    for iRun=1:size(func_Realign_cell, 2)
        for iScan=1:size(func_Realign_cell{1,iRun}, 1)

            scan_echo1_1 = deblank(func_Realign_cell{1,iRun}(iScan)); % Pick out current image
            % if iScan==1
            %     M = spm_get_space(scan_echo1_1{:}); % Read its voxel-to-world info
            % else
            %     scan_echo1_2 = deblank(func_Realign_cell{1,iRun}(iScan-1)); % Pick out current image
            %     M1 = spm_get_space(scan_echo1_1{:}); % Read its voxel-to-world info
            %     M2 = spm_get_space(scan_echo1_2{:}); % Read its voxel-to-world info
            %     M = (M1+M2)./2;
            % end
            M = spm_get_space(scan_echo1_1{:}); % Read its voxel-to-world info
            scan_echo2 = deblank(func_Realign_cell{2, iRun}(iScan)); % Pick out current image
            spm_get_space(scan_echo2{:}, M); % write this matrix into the header of the 2nd TE time-series
            clear scan_echo1_1 scan_echo2 M M1 M2
        end
    end
end


% reslice all the scans
realign_write = struct;
realign_write.matlabbatch{1}.spm.spatial.realign.write.data=vertcat(func_Realign_cell{:});
realign_write.matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 1];
realign_write.matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
realign_write.matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
realign_write.matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
realign_write.matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
cfg_util('run',realign_write.matlabbatch);




%%
% Unwarping

% spm fmri
% realign_estimate_unwarp = struct; % realign and unwarp
% data_realign = func_Realign_cell(1,:);
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.data.scans = vertcat(data_realign{:});
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan = '';
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
% realign_estimate_unwarp.matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'ru';
% cfg_util('run',realign_estimate_unwarp.matlabbatch);


% reslice all the scans
% realign_write = struct; % realign and unwarp
% realign_write.matlabbatch{1}.spm.spatial.realignunwarp.write.data=vertcat(func_Realign_cell{:});
% matlabbatch{1}.spm.spatial.realignunwarp.write.uwroptions.uwwhich = [2 1];
% matlabbatch{1}.spm.spatial.realignunwarp.write.uwroptions.rinterp = 4;
% matlabbatch{1}.spm.spatial.realignunwarp.write.uwroptions.wrap = [0 0 0];
% matlabbatch{1}.spm.spatial.realignunwarp.write.uwroptions.mask = 1;
% matlabbatch{1}.spm.spatial.realignunwarp.write.uwroptions.prefix = 'ru';
% cfg_util('run',realign_write.matlabbatch);
