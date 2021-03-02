#! /bin/bash -x
d0=20130101
d1=20130531
dt=$d0
while [ $dt -le $d1 ]
do
    singularity exec --cleanenv $PYNEXTSIM_SIF ./proc_ecmwf_fc.py $dt
    dt=$(date -d "$dt +1day" "+%Y%m%d")
done
