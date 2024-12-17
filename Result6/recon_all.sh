#!/usr/bin/env bash
export SUBJECTS_DIR=/home/sunyue/recon-UPENN
for subj in `ls ./test`
do
    recon-all -s $subj -i ./test/$subj/*.nii -sd /home/sunyue/recon-UPENN/ -all -qcache
done

