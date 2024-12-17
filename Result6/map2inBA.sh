#!/usr/bin/env bash
export SUBJECTS_DIR=/home/sunyue/recon1/

for subj in `ls ./recon1 -I fsaverage`
do
	mri_surf2surf --hemi lh --srcsubject fsaverage --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.BA_aparc.annot \
		      --trgsubject $subj --trgsurfval $SUBJECTS_DIR/$subj/label/lh.BA.annot

	mris_anatomical_stats -a $SUBJECTS_DIR/$subj/label/lh.BA.annot \
		      -f $SUBJECTS_DIR/$subj/stats/lh.aparc.BAatlas.stats -b $subj lh

	mri_surf2surf --hemi rh --srcsubject fsaverage --sval-annot $SUBJECTS_DIR/fsaverage/label/rh.BA_aparc.annot \
		      --trgsubject $subj --trgsurfval $SUBJECTS_DIR/$subj/label/rh.BA.annot

	mris_anatomical_stats -a $SUBJECTS_DIR/$subj/label/rh.BA.annot \
		      -f $SUBJECTS_DIR/$subj/stats/rh.aparc.BAatlas.stats -b $subj rh
done