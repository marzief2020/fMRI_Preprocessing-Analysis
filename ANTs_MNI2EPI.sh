#!/bin/bash

# ANTs bash code for MNI --> T1 --> EPI

# subjects folders are in $base_dir/subjects: sub-01, sub-02, ....
export base_dir=/mnt/d/ANTs

# iterate over each subject, transform ROIs in the MNI space to the functional space of each participant via T1
for sub_dir in $base_dir/subjects/*; do

	
	subj_id=$(basename $sub_dir)
	
	# define ouput folder for the registered images 
	export out_dir=$sub_dir/ANTs_outputs
	[ ! -d "$out_dir" ] && mkdir -p "$out_dir"
	
	# define output for the registered ROIs 
	export ROI_dir=$out_dir/ROIs
	[ ! -d "$ROI_dir" ] && mkdir -p "$ROI_dir"

		# registering the mean epi image to T1 and T1 to MNI space 
		antsRegistrationSyN.sh -d 3 -f $sub_dir/*T1w_Brain.nii -m $sub_dir/meansub*Brain.nii -o $out_dir/boldToT1_ -t r 
		antsRegistrationSyN.sh -d 3 -f $base_dir/tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz  -m $sub_dir/*T1w_Brain.nii -o $out_dir/t1ToMNI_res-01_  -t s 
		
		# get the ROIs in the MNI space (in MNI_ROIs folder), and transform them to the functional space 
		for ROI_files in  $base_dir/MNI_ROIs/*.nii; do
			ROI_name=$(basename $ROI_files)
		 
			antsApplyTransforms \
				-d 3 \
				-i $base_dir/MNI_ROIs/"$ROI_name" \
				-o $ROI_dir/"r${subj_id}_${ses}_ROI_${ROI_name}" \
				-n GenericLabel \
				-r $sub_dir/meansub*Brain.nii \
				-t [ $out_dir/boldToT1_0GenericAffine.mat, 1 ] \
				-t [ $out_dir/t1ToMNI_res-01_0GenericAffine.mat, 1 ] \
				-t $out_dir/t1ToMNI_res-01_1InverseWarp.nii.gz	
				
				# -- interpolation NearestNeighbor \
				
		done
		
done