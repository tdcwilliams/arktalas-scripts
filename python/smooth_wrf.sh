#! /bin/bash -x
set -e
indir=/cluster/work/users/timill/wrf_downscaled
mkdir -p $indir
if [ 1 -eq 0 ]
then
    rm -f $indir/*
    cp /cluster/projects/nn2993k/sim/data/WRF/start_20130213_nudging/*r10* $indir
fi

for n in 4 8
do
    outdir=/cluster/work/users/timill/wrf_downscaled_smoothed_${n}0km/
    mkdir -p $outdir
    rm -f $outdir/*

    mkdir -p $outdir
    for f in $indir/*r10_nudge*
    do
        bname=$(basename $f)
        s=($(echo $bname | grep "LANDMASK")) || \
            singularity exec --cleanenv $PYNEXTSIM_SIF \
            ./smooth_wrf.py $f $outdir/$bname -sf $n
    done
done
