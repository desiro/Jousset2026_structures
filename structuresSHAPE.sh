#/bin/bash

################################################################
## run with installation:
# bash run.sh install

## run without installation:
# bash run.sh

################################################################
## paths
mdir="."

# scripts
getStSHP="$mdir/getStructuresSHAPE.py"
cmpStSHP="$mdir/compareStructuresSHAPE.py"

# data
data="$mdir/data"

# results
results="$mdir/results"
resGetStSHP="$results/getStructuresSHAPE"
resCmpStSHP="$results/compareStructuresSHAPE"

# conda
if [ "$1" = "install" ]; then
  conda env create -f "$mdir/envStructuresSHAPE.yml"
fi

# mkdir
mkdir -p $resGetStSHP $resCmpStSHP

################################################################
## main
function main() {
    ################################################################
    echo "# Step 0: Get structures for SHAPE."
    files=("SC35M-WT" "SC35M-rS1mut" "SC35M-rS3mut")
    getStructuresSHAPE files

    echo "# Step 1: Compare generated structures."
    files1=("SC35M-WT" "SC35M-WT")
    files2=("SC35M-rS1mut" "SC35M-rS3mut")
    compareStructuresSHAPE files1 files2
}

################################################################
## functions
function getStructuresSHAPE() {
    local -n one=$1
    for i in "${!one[@]}"
    do
        conda run -n strctures_SHAPE python $getStSHP \
            --prefix $resGetStSHP/${one[i]} --overwrite --pickleData --viennaRNAtemperature 32.0 \
            --structureEnergy -6.0 --interactionDistance 200 \
            --genomeFasta $data/${one[i]}.fa --dataSHAPE $data/${one[i]}
    done
}

function compareStructuresSHAPE() {
    local -n one=$1
    local -n two=$2
    for i in "${!one[@]}"
    do
        conda run -n strctures_SHAPE python $cmpStSHP \
            --prefix $resCmpStSHP/${one[i]}_${two[i]} --overwrite --pickleData --viennaRNAtemperature 32.0 \
            --genomeFasta1 $data/${one[i]}.fa --dataSHAPE1 $data/${one[i]} --structureFile1 $resGetStSHP/${one[i]}/${one[i]}_subsets.pcl \
            --genomeFasta2 $data/${two[i]}.fa --dataSHAPE2 $data/${two[i]} --structureFile2 $resGetStSHP/${two[i]}/${two[i]}_subsets.pcl
    done
}

################################################################
## call
main
