#!/bin/bash

###############################################################################################################################################
## This script is used to prepare the protein prior to inserting it into the lipid membrane. In this step, firstly, the protein is cleaned up 
## and then the protonation states of the residues and termini are assigned. At the end, the orientation of the protein with respect to a flat
## DOPC membrane is determined using PPM3 software. The software is locally installed. The output PDB structure is inserted into a POPC lipid 
## bilayer. Alternatively, users can download the pre-oriented protein from OPM database and then inserted into the bilayer. 

## Usage:

##    bash Protein_Prep.sh --input input.pdb --dest ./ --chain AB

## Written by Nidhin Thomas (dr.thomas@nidhin-thomas.com)

###############################################################################################################################################

usage() { echo -e "\nThis script is used to prepare the protein for the MD simulations and then insert into a POPC lipid membrane using PACKMOL.\n";

          echo -e "bash Protein_Prep.sh [--input <input.pdb>] [--dest <destination directory>] [--chain <chain IDs of the protein>] [--help]\n" 1>&2;

          echo -e "Example: bash Protein_Prep.sh --input input.pdb --dest ./  --chain AB \n";

          exit 1; 
        }

for arg in "$@"; do
  shift
  case "$arg" in
    '--help')      set -- "$@" '-h'    ;;
    '--input')     set -- "$@" '-f'    ;;
    '--dest')      set -- "$@" '-d'    ;;
    '--chain')      set -- "$@" '-c'   ;;
    *)             set -- "$@" "$arg"  ;;
    
  esac
done

# Default behavior
dest="./";
input="";
chain="A";

OPTIND=1
while getopts "hf:d:c:" opt
do
  case "$opt" in
    'h') usage; exit 0 ;;
    'f') input=$OPTARG ;;
    'd') dest=$OPTARG  ;;
    'c') chain=$OPTARG ;;
    '?') usage >&2; exit 1 ;;
  esac
done
shift $(expr $OPTIND - 1)

script_dir=$( dirname -- "$0"; )
ppm_dir=${script_dir}/ppm

conda activate MemProtMD

## If the code is run in a linux cluster, ensure that PPM can be run from the command line.
module load gcc
module load amber/22
export AMBERHOME=/amber/22

function pqr2pdb () {
inputPDB=$1
script_dir=$2
destination=$3
pdbfilename=$(basename "$inputPDB")
pdbname=${pdbfilename%.*}
pdb2pqr30 --with-ph 7.4 --ff=AMBER --titration-state-method propka --ffout=AMBER --pdb-output ${destination}/${pdbname}_pka7.pdb ${inputPDB} ${destination}/${pdbname}_pka7.pqr
propka3 ${inputPDB}
mv ${pdbname}.pka ${destination}/${pdbname}.pka
python ${script_dir}/pdbGMXprep.py --pkafile ${destination}/${pdbname}.pka --pH 7.4 --mutatedPDB ${destination}/${pdbname}_protonated.pdb --pdb ${destination}/${pdbname}_pka7.pdb
}

function ppmOrient () {
inputPDB=$1
ppm_dir=$2
script_dir=$3
destination=$4
chain=$5
pdbfilename=$(basename "$inputPDB")  
pdbname=${pdbfilename%.*}
python ${script_dir}/prepPPMinput.py --pdb ${inputPDB} --ppm ${ppm_dir} --dest ${destination} --chain ${chain}
${ppm_dir}/immers<${destination}/${pdbname}.in>${destination}/${pdbname}.out
sed '/HETATM|DUM/d' ${destination}/${pdbname}out1.pdb > ${destination}/${pdbname}_ppm3.pdb
packmol-memgen --pdb ${destination}/${pdbname}_ppm3.pdb --lipids POPC --ratio 1 --preoriented --salt --salt_c Na+ --saltcon 0.15 --dist 10 --dist_wat 15 --notprotonate --nottrim
}

pqr2pdb ${input} ${script_dir} ${dest}
ppmOrient ${dest}/*_protonated.pdb ${ppm_dir} ${script_dir} ${dest} ${chain}
