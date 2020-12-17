#! /bin/bash -x
# 
# In new container do:
# export DEBIAN_FRONTEND=noninteractive
# apt-get update
# apt-get install texlive-pictures texlive-science texlive-latex-extra
# pip install pylatex

[[ $# -ne 2 ]] && exit 1
REMOTE_DIR=$1
LOCAL_DIR=$2

# download the mean maps for observations in $obs
obs=()
obs+=(eval-cs2smos)
obs+=(eval-osisaf-conc)
obs+=(eval-osisaf-drift)
obs+=(eval-osisaf-drift-mu10kpd)

pdf_list=()
for obs_dir in "${obs[@]}"
do
    remote_dir=$REMOTE_DIR/$obs_dir
    local_dir=$LOCAL_DIR/$obs_dir
    for subdir in maps_bias maps_fcst maps_obs
    do
        mkdir -p $local_dir/$subdir
        scp $remote_dir/$subdir/*mean* $local_dir/$subdir
    done

    # download the summary figs
    scp $remote_dir/*.* $local_dir/

    # collect figs into pdf
    python pylatex_nxs.py $local_dir
    pdf=$local_dir/collected_maps.pdf
    [[ -f $pdf ]] && pdf_list+=($pdf)
done

# merge into one pdf
[[ ${#pdf_list[@]} -lt 2 ]] && exit 0
suffix=$(basename $LOCAL_DIR)
pdfunite ${pdf_list[@]} $LOCAL_DIR/collected_maps_$suffix.pdf
