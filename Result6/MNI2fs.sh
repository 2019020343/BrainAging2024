#!/usr/bin/env bash
export SUBJECTS_DIR=/home/sunyue/recon1/
cd $SUBJECTS_DIR

mri_vol2surf --mov /home/sunyue/recon1/HOA112_1mm.nii --mni152reg --hemi lh --out_type mgh --o lh.HOA112.mgh
for ((j=1;j<=112;j+=2))
do
    mri_vol2label --i lh.HOA112.mgh --id $j --l /home/sunyue/recon1/fsaverage/label/lh.num$j.label --surf fsaverage lh
done

mris_label2annot --s fsaverage --h lh --ctab /home/sunyue/MyFreeSurferColorLUT.txt --a myaparc \
				 --l $SUBJECTS_DIR/fsaverage/label/lh.num1.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num3.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num5.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num7.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num9.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num11.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num13.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num15.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num17.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num19.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num21.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num23.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num25.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num27.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num29.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num31.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num33.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num35.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num37.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num39.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num41.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num43.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num45.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num47.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num49.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num51.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num53.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num55.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num57.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num59.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num61.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num63.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num65.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num67.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num69.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num71.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num73.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num75.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num77.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num79.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num81.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num83.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num85.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num87.label \
				 --l $SUBJECTS_DIR/fsaverage/label/lh.num89.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num91.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num93.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num95.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num97.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num99.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num103.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num107.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num109.label \
                 --l $SUBJECTS_DIR/fsaverage/label/lh.num111.label 




mri_vol2surf --mov /home/sunyue/recon1/HOA112_1mm.nii --mni152reg --hemi rh --out_type mgh --o rh.HOA112.mgh
for ((j=2;j<=112;j+=2))
do
    mri_vol2label --i rh.HOA112.mgh --id $j --l /home/sunyue/recon1/fsaverage/label/rh.num$j.label --surf fsaverage rh
done

mris_label2annot --s fsaverage --h rh --ctab /home/sunyue/MyFreeSurferColorLUT3.txt --a myaparc \
				 --l $SUBJECTS_DIR/fsaverage/label/rh.num2.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num4.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num6.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num8.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num10.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num12.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num14.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num16.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num18.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num20.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num22.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num24.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num26.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num28.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num30.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num32.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num34.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num36.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num38.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num40.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num42.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num44.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num46.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num48.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num50.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num52.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num54.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num56.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num58.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num60.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num62.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num64.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num66.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num68.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num70.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num72.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num74.label \
				 --l $SUBJECTS_DIR/fsaverage/label/rh.num76.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num78.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num80.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num82.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num84.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num86.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num88.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num90.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num92.label \
				 --l $SUBJECTS_DIR/fsaverage/label/rh.num94.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num96.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num98.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num100.label \
				 --l $SUBJECTS_DIR/fsaverage/label/rh.num108.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num110.label \
                 --l $SUBJECTS_DIR/fsaverage/label/rh.num112.label