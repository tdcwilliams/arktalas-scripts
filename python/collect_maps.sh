#! /bin/bash -x
# needs to run on a machine or inside a container
# with python (specifically the pylatex library) and pdfunite installed

[[ $# -ne 2 ]] && exit 1
ROOT_DIR=$1

# download the mean maps for observations in $obs
obs=()
obs+=("eval-cs2smos")
obs+=("eval-osisaf-conc")
obs+=("eval-osisaf-drift")
obs+=("eval-osisaf-drift-mu10kpd")

pdf_list=()
for sub_dir in "${obs[@]}"
do
    obs_dir=$ROOT_DIR/$sub_dir
    [[ ! -d $obs_dir ]] && continue

    # collect figs into pdf
    ./collect_maps_pylatex.py $obs_dir
    pdf=$obs_dir/collected_maps.pdf
    [[ -f $pdf ]] && pdf_list+=($pdf)
done

# merge into one pdf
[[ ${#pdf_list[@]} -eq 0 ]] && echo no pdf files created - quitting && exit 0
outfile=$ROOT_DIR/collected_maps_$(basename $ROOT_DIR).pdf
[[ ${#pdf_list[@]} -eq 1 ]] && cp ${pdf_list[0]} $outfile && exit 0
pdfunite ${pdf_list[@]} $outfile
