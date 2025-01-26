function PreProc_MaskGeneration(structural_img, main_img_Brain, ref_mask, task_dir, ses)

mask_threshold = 0.5;
ref_mask_info = spm_vol(ref_mask);

[d, e, f] = fileparts(structural_img);
output.rstructural_fn = [d filesep 'r' e f];
output.rgm_fn = [d filesep 'rc1' e f];
output.rwm_fn = [d filesep 'rc2' e f];
output.rcsf_fn = [d filesep 'rc3' e f];
%     output.rbone_fn = [d filesep 'rc4' f e];
%     output.rsoft_fn = [d filesep 'rc5' f e];
%     output.rair_fn = [d filesep 'rc6' f e];

[GM_img_bin, WM_img_bin, CSF_img_bin] = createBinarySegments(output.rgm_fn, output.rwm_fn, output.rcsf_fn, mask_threshold);
I_GM = find(GM_img_bin); I_WM = find(WM_img_bin); I_CSF = find(CSF_img_bin);


% brain_mask_raw = GM_img_bin | WM_img_bin | CSF_img_bin;
brain_mask_raw = spm_read_vols(spm_vol(main_img_Brain));
brain_mask_raw(brain_mask_raw>0)=1;
brain_mask_raw(brain_mask_raw<0)=0;
brain_mask = imfill(brain_mask_raw, 'holes');
se = strel('disk', 1); brain_mask_dil = imfill(imdilate(brain_mask_raw, se), 'holes');

% brain_mask_info=spm_vol([ref_img ',1']);
% [d_fn,e_fn,f_fn]=fileparts(ref_img);
brain_mask_info = ref_mask_info;
brain_mask_info.fname=[task_dir filesep [e(1:6) '_' ses '_brain_mask.nii']];
spm_write_vol(brain_mask_info, brain_mask);

brain_mask_dil_info = ref_mask_info;
brain_mask_dil_info.fname=[task_dir filesep [e(1:6) '_' ses '_brain_mask_dil.nii']];
spm_write_vol(brain_mask_dil_info, brain_mask_dil);

GM_mask_info = ref_mask_info;
GM_mask_info.fname=[task_dir filesep [e(1:6) '_' ses '_GM_mask.nii']];
spm_write_vol(GM_mask_info, GM_img_bin);

WM_mask_info = ref_mask_info;
WM_mask_info.fname=[task_dir filesep [e(1:6) '_' ses '_WM_mask.nii']];
spm_write_vol(WM_mask_info, WM_img_bin);

CSF_mask_info = ref_mask_info;
CSF_mask_info.fname=[task_dir filesep [e(1:6) '_' ses '_CSF_mask.nii']];
spm_write_vol(CSF_mask_info, CSF_img_bin);

se = strel('disk',1); 
CSF_mask_erode =  imerode(CSF_img_bin,se);
CSF_mask_erode_info = ref_mask_info;
CSF_mask_erode_info.fname=[task_dir filesep [e(1:6) '_' ses '_CSF_mask_erode.nii']];
spm_write_vol(CSF_mask_erode_info, CSF_mask_erode);

se = strel('disk',1); 
WM_mask_erode =  imerode(WM_img_bin,se);
WM_mask_erode_info = ref_mask_info;
WM_mask_erode_info.fname=[task_dir filesep [e(1:6) '_' ses '_WM_mask_erode.nii']];
spm_write_vol(WM_mask_erode_info, WM_mask_erode);

% combining the white and CSF masks and erode them (to prevent regressing out the GM signal during aCompcorr)
noise_mask_info = ref_mask_info;
noise_mask = WM_img_bin | CSF_img_bin;
se = strel('disk',1); 
noise_mask_erode = imerode(noise_mask,se);
noise_mask_info.fname=[task_dir filesep [e(1:6) '_' ses '_noise_mask.nii']];
spm_write_vol(noise_mask_info, noise_mask_erode);

 