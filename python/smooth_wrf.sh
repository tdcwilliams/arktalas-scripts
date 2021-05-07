#! /bin/bash -x
indir=/cluster/work/users/timill/wrf_downscaled
outdir=/cluster/work/users/timill/wrf_downscaled_smooothed/
prefix=breakup_march2013_r10_ctrl_
for infile in $indir/$prefix*
do
    ./smooth_wrf.py $infile $outdir/$(basename $infile) -sf 8
done
