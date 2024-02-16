# MembProtMD_AMBER
==============================

A workflow to automate MD simulations of membrane-protein complexes using AMBER

## Membrane-Protein MD simulations in AMBER
This workflow is designed to perform equilibrium MD simulations of proteins embedded in lipid bilayers. User is expected to have fundamental knowledge in molecular dynamics and structural chemistry as a pre-requisite to run this workflow. The workflow is divided into five stages:  

1. Initial input protein preparation using PROPKA+PDB2PQR  
2. Generating the correct orientation of proteins within lipid membranes using PPM local version  
3. Insertion of proteins into lipid membranes of choice using PACKMOL-MEMGEN  
4. Build system topology using tleap  
5. Equilibrium MD simulation of complex system in AMBER PMEMD
6. Post-simulation analyses of the complex (cluster analysis, PCA, Free energy, Structure analysis)    

### 1. Initial input protein preparation using PROPKA+PDB2PQR
The input protein structure may have problems with the protonation state of the residues, terminis and disulfide bond definition. These problems could originate by the way these structures were generated. Therefore, the input structures need to be cleaned well prior to any modeling. I highly recommend using preliminary cleaning of the pdb file using tools such as "pdb-tools" or Parmed or any other simple PDB editing tools. Structures from Rosetta may have residue-level energy values written at the end of the PDB files. These lines should be deleted during this stage to obtain reasonable structures without failure. Just like GROMACS workflow, the input PDB structure should be checked for proper protonation. If the protonation states of the residues are different from the original files, user can essentially see those suggestions from PROPKA and PDB2PQR. If the protonation state, recommended is different than the original information, pdbGMXprep.py will read the output from PROPKA and modify the input PDB files to reflect the correct protonation states suggested by PROPKA. GROMACS typically determines the salt bridge states during "pdb2gmx". However, the disulfide bonds in the proteins or between proteins, in AMBER MD package are determined during the initial phase itself. Correct naming of the chain IDs, residues, and termini are required to successfully run these steps.  

<pre>bash
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
    
### 4. Build system topology using tleap
The system topology files contain the force field parameters required to perform the MD runs. In AMBER, there are multiple ways to build the topology. Here, we use tleap to build the topology. The system topology can be built either interactively or in single step. Users are encouraged to use the interactive version of tleap from command line for new systems. However, if the runs need to be parallelized, topology building in single step is preferred. In order to run tleap, In Lovelace and Marjorie, users need to load AMBER package.  

For interactive session:  
`module load amber/22`
`tleap`

Enter each line in topol_build.leap one by one. Modify the lines if user has additional molecules in the system or  user wants to change the system force field parameters.  

For execution of topology building in single step:

`tleap -f system_build.leap`

Each line in the topol_build.leap is explained below.  

**source leaprc.protein.ff19SB**: AMBERFF19SB force field is used for protein-in-solution  
**source leaprc.water.tip3p**: TIP3P water molecules are added to solvate the system  
**source leaprc.gaff2**: If small molecules or ligands are added to the system, they need be to parametrized. In AMBER, we can use GAFF2 to parametrize these molecules using AM1-BCC theory.  
**source leaprc.lipid21**: Latest lipid force field parametes used by AMBER. The choice of the lipid force field is dependent on the combination of FF being used in the system.  
**loadamberparams frcmod.ionsjc_tip3p**: Monovalent ion parameters for Ewald and TIP3P  
**water.receptor = loadpdb protein.pdb**: Adding the transmembrane protein to the FF  
**membrane = loadpdb POPC_amber.pdb**: Adding the membrane lipid to the FF  
**rism_wat = loadpdb rism_sele.pdb**: Adding the initial water molecules to the FF  
**system = combine{ receptor membrane rism_wat }**: Combining these FF parametes together  
**set system box {XX.XXX XX.XXX XX.XXX }**: System box size is added  
**saveamberparm system protein_membrane.prmtop protein_membrane.inpcrd**: Saving the parameters  
**savepdb system protein_membrane.pdb**: Save the PDB files after building the topology  
**quit**: Quit tleap after saving the files  

**Hydrogen Mass Repartitioning** : Optionally user can also increase the simulation timestep to 4fs from 2fs with the help of hydrogen mass repartitioning. This can be done with the help of the Parmed package.  

`parmed -i hmass_parmed.in`

### 5. Equilibrium MD simulation of complex system in AMBER PMEMD
Equilibrium MD simulation can be performed as a single package. It includes initial energy minimization, NVT equilibration, NPT equilibration and production run. We perform these steps in GPU (minimization can be done in CPU). This step can either be a part of the entire workflow or separate. If user is exploring the structure and setting up the system for the first time, we recommend using the stand alone bash script "run_MD.sh" to perform the simulation. This step gives control over the type of input files and the simulation parameters such as timestep, thermostat and barostat variables. The simulation steps and associated input files are discussed below.  

**01_Min.in**: very short minimization on CPU with pmemd. This is advised for membrane systems, which often have bad clashes of lipid chains to resolve. The CPU code is more robust in dealing with these than the GPU  
**02_Min2.in**: longer minimization on GPU  
**03_Heat.in, 04_Heat2.in**: heating to 303 K, with restraints on protein and membrane lipids  
**05_Back.in**: run 1 ns NPT with restraints on receptor backbone atoms  
**06_Calpha.in**: run 1 ns NPT with restraints on receptor carbon-alpha atoms  
**07_Prod.in**: run 100 ns NPT equilibration, all restraints removed  
**08_Long.in**: run 500 ns NPT production, with Monte Carlo barostat and hydrogen mass repartitioning  

The script can be run from command line:  

`sbatch run_MD.sh`

### 6. Post-simulation analyses of the complex