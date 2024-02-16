# MemProt4MD

A tool to embed transmembrane proteins into lipid membranes for MD simulations using PPM and Packmol.

## MemProt4MD
This tool is designed to generate the input structures for perform equilibrium MD simulations of proteins embedded in lipid bilayers. User is expected to have fundamental knowledge about proteins, lipid membranes, and experience with python and bash scripting. This tool integrates two tools commonly used for generating Protein-Lipid membrane complexes: PPM3 and Packmol. Users can either orient the protein with respect to a dummy lipid membrane using this tool or download the pre-oriented protein from PPM 3.0 web server. The tool includes the following:

1. Initial input protein preparation using PROPKA+PDB2PQR  
2. Generating the correct orientation of proteins within lipid membranes using PPM local version  
3. Insertion of proteins into lipid membranes of choice using PACKMOL-MEMGEN  
4. Protein-membrane complex visualization   

## How is the repository arranged?

**src**: src contains the scripts relevant for this tool. Protein_Prep.sh is the master wrapper bash sript used in this tool.

**example**: example directory contains the input PDB file and output files from PPM and Packmol tools.

**notebook**: Python notebook for visualization and validation of the structure.

## Prerequisites

PPM local version can be downloaded from `https://opm.phar.umich.edu/ppm_server3`. 

Install conda environment using:  

`conda env create --file environment.yml`

Activate the conda environment:  

`conda activate MemProtMD`

Packmol is provided with AmberTools. AmberTools is installed along with Amber software. AmberTools can also installed through Anaconda and is already included in the conda environment.

### 1. Initial input protein preparation using PROPKA+PDB2PQR
The input protein structure may have problems with the protonation state of the residues, terminis and disulfide bond definition. These problems could originate by the way these structures were generated. Therefore, the input structures need to be cleaned well prior to any modeling. I highly recommend using preliminary cleaning of the pdb file using tools such as "pdb-tools" or Parmed or any other simple PDB editing tools. Structures from Rosetta may have residue-level energy values written at the end of the PDB files. These lines should be deleted during this stage to obtain reasonable structures without failure. Just like GROMACS workflow, the input PDB structure should be checked for proper protonation. If the protonation states of the residues are different from the original files, user can essentially see those suggestions from PROPKA and PDB2PQR. If the protonation state, recommended is different than the original information, pdbGMXprep.py will read the output from PROPKA and modify the input PDB files to reflect the correct protonation states suggested by PROPKA. GROMACS typically determines the salt bridge states during "pdb2gmx". However, the disulfide bonds in the proteins or between proteins, in AMBER MD package are determined during the initial phase itself. Correct naming of the chain IDs, residues, and termini are required to successfully run these steps.  

<pre>
```
function pqr2pdb () {
inputPDB=$1
script_dir=$2
destination=$3
pdbfilename=$(basename "$inputPDB")  
pdbname=${pdbfilename%.*}
pdb2pqr30 --with-ph 7.4 --ff=AMBER --titration-state-method propka --ffout=AMBER \ 
    --pdb-output ${destination}/${pdbname}_pka7.pdb ${inputPDB} ${destination}/${pdbname}_pka7.pqr
propka3 ${inputPDB}
mv ${pdbname}.pka ${destination}/${pdbname}.pka
python ${script_dir}/pdbGMXprep.py --pkafile ${destination}/${pdbname}.pka --pH 7.4 \
    --mutatedPDB ${destination}/${pdbname}_protonated.pdb --pdb ${destination}/${pdbname}_pka7.pdb
}  
```
</pre>

### 2. Generating the correct orientation of proteins within lipid membranes using PPM3 local version
Typically, the orientation of transmembrane proteins in lipid membranes are determined from OPM database or by running it in PPM server. Howeve, they also provide a local version of this software and that has been added to this pipeline under "PPM" directory. The correct orientation of the protein with respect to a planar membrane is determined by embedding the protein into a DOPC bilayer. The details of the protocol and the collection of output results are provided in the link below. You can also find a detailed instruction from the orignal authors in a word document in PPM directory under the name "ppm3_instructions.docx".PPM3 requires gfortran. The format of the input file (*.inp) is given in "example" directory. The example bash script is given below. The executable and input files are all expected to be in the same working directory.

`module load gcc`

`./immers<2membrane.inp>output_file` 

The output PDB file contains the re-oriented protein structure along with dummy beads to represent the membrane-water interface (typically phosphorus atoms in the DOPC or Carbonyl group at the interface). These "HETATM" lines need to be removed before going into the next stage where the protein is inserted into a lipid membrane with composition of our choice. These lines can be removed using the "sed" command.

`sed '/HETATM\|DUM/d' input.pdb > output.pdb`

### 3. Insertion of proteins into lipid membrane using PACKMOL-MEMGEN
PACKMOL-MEMGEN is pre-installed in AMBER module. When AMBER is loaded onto the terminal, user can use PACKMOL-MEMGEN to generate either membrane-only systems or proteins embeddged into the lipid membrane. The input arguments required for this stage can be found using the commands below. The orientation the protein with respect to lipid membrane is already determined using PPM3-server local version.The descriptions of the command line arguments are provided below:

**--pdb pdbname_from_ppm3.pdb**: input transmembrane protein PDB  
**--lipids POPC**: build POPC membrane (Composition of the lipids can be modified)  
**--ratio 1**: POPC ratio of 1 (Different composition can be tried here)  
**--preoriented**: the PPM coordinates orientated along the z-axis.  
**--salt --salt_c Na+ --saltcon 0.15**: Add 0.15 M of NaCl salt to the water layer  
**--dist 10**: minimum maxmin value for x, y, z to box boundary  
**--dist_wat 15**: water layer thickness of 15 Angstrom  
**--notprotonate --nottrim**: No need to further process input protein PDB file  

`packmol-memgen --helppackmol-memgen --pdb destination/pdbname_from_ppm3.pdb --lipids POPC --ratio 1 --preoriented --salt --salt_c Na+ --saltcon 0.15 --dist 10 --dist_wat 15 --notprotonate --nottrim`

### 4. Visualization

To ensure that the protein is properly embedded in the bilayer and the system is solvated appropriately, visualize the system using nglview or other tools such as PyMol or UCSF ChimeraX.