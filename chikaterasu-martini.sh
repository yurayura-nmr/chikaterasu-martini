#!/bin/sh

# =============================================================================
# Chikaterasu Martini - Automated Coarse-Grained MD Setup
# 
# Author: Erik Walinda
# Affiliation: Kyoto University, Graduate School of Medicine
# 
# Description: Automated setup for Martini 2.2 force field simulations
# Tested with: GROMACS 2022, 2025
# Repository: https://github.com/yurayura-nmr/chikaterasu-martini
# =============================================================================

# =============================================================================
# CONFIGURATION PARAMETERS
# 
# Modify these variables to configure your simulation run.
# Most MD parameters are defined in the mdp files.
# =============================================================================

# Debug level controls execution stopping points:
#   0 - Full production MD run [NPT]; no intermediate stops
#   1 - Stop after topology generation [pdb2gmx]
#   2 - Stop after solvation
#   3 - Stop after counterion addition  
#   4 - Stop before production MD (after equilibration)
debug_level=1

# Number of independent trajectories to simulate (range: 1-10)
nruns=1

# Manual box dimensions configuration
# Set to false for automatic box size calculation
box_manual=true
box_dim="7.500   7.500   7.500"

# File naming conventions
protein_name="Abeta-CG"  # PDB filename (without extension)
top="Abeta"              # Topology filename (without extension)

# =============================================================================
# SIMULATION PREREQUISITES
#
# Required files for system setup:
#   - ./gromacs/coord/${protein_name}.pdb     : Starting conformation
#   - ./${top}.top                            : Main topology file
#   - ./Protein_A+Protein_B.itp               : Protein chain topologies
#   - ./chika_mdp/martini.itp                 : Martini forcefield parameters
#
# Example for K48-linked diubiquitin:
#   - ./k48.top
#   - ./Protein_A+Protein_B.itp  
#   - ./gromacs/coord/k48-cg.pdb
# =============================================================================

# =============================================================================
# COMMAND LINE ARGUMENT PARSING
#
# If debug level is provided as command line argument, it overrides the
# configured value above.
# =============================================================================
if [ -z "$1" ]; then
    read -p "[Chikaterasu-martini] Command line arguments are empty. Using manually set debug level $debug_level." dummy
else
    read -p "[Chikaterasu-martini] Command line arguments provided. Using first argument as debug level $1." dummy
    debug_level=$1
fi

# =============================================================================
# DIRECTORY STRUCTURE SETUP
#
# Creates and cleans necessary directory structure for simulation files.
# Previous run directories are removed to ensure clean state.
# =============================================================================
mkdir -p gromacs
rm -rf gromacs/top
rm -rf gromacs/solvation
rm -rf gromacs/addions
rm -rf gromacs/emin

mkdir -p gromacs/top
mkdir -p gromacs/solvation
mkdir -p gromacs/addions
mkdir -p gromacs/emin
mkdir -p gromacs/coord
mkdir -p runs
mkdir -p runs/nvt
mkdir -p runs/npt
mkdir -p runs/md
mkdir -p custom_analysis

# =============================================================================
# TOPOLOGY GENERATION NOTE
#
# IMPORTANT: This script assumes topology has been pre-generated using martinize.py
# 
# Example commands for topology generation:
#   K48-linked diubiquitin (merged chains):
#     python2.7 martinize.py -f 1aar_modified.pdb -o 1aar_modified.top \
#              -x 1aar_modified-CG.pdb -dssp dssp -p backbone -merge A,B
#
#   Monoubiquitin:
#     python2.7 martinize.py -f 1UBQ.pdb -o single-ubq.top \
#              -x 1UBQ-CG.pdb -dssp dssp -p backbone
# =============================================================================

# =============================================================================
# BOX CREATION AND VACUUM MINIMIZATION
#
# Steps:
#   1. Create simulation box around protein using editconf
#   2. Optionally apply manual box dimensions
#   3. Prepare and run vacuum energy minimization
# =============================================================================

echo "[Chikaterasu-martini] Creating simulation box and running vacuum minimization..."

# Create initial simulation box
cd gromacs/coord
gmx editconf -f $protein_name.pdb -d 1.0 -bt triclinic -o $protein_name.gro

# Apply manual box dimensions if configured
if [ "$box_manual" = true ]; then
    sed -i '$ d' $protein_name.gro
    echo "$box_dim" >>$protein_name.gro
    echo "[Chikaterasu-martini] Applied manual box dimensions: $box_dim"
fi

# Prepare vacuum minimization run
cd ../top/
cp ../../chika_mdp/martini.itp .
cp ../../*.top .
cp ../../*.itp .

gmx grompp -p $top.top -f ../../chika_mdp/minimization.mdp -c ../coord/$protein_name.gro -o ../emin/minimization-vac.tpr

# Execute vacuum minimization
cd ../emin/
gmx mdrun -deffnm minimization-vac -v

# Exit if debug level 1 (stop after topology generation)
if [ "$debug_level" = 1 ]; then
    echo "[Chikaterasu-martini] Debug level 1 reached - exiting after topology generation"
    exit 0
fi

: '
*************************************************************
Solvate the protein and perform minization of solvated system
*************************************************************
'

cd ../solvation
cp ../../chika_mdp/water.gro .
gmx solvate -cp ../emin/minimization-vac.tpr -cs water.gro -radius 0.21 -o system_solvated.gro

# Now we need to count how many waters we added
cp ../top/$top.top ../top/system_solvated.top
echo -n "\nW\t\t" >>../top/system_solvated.top
grep -c W system_solvated.gro >>../top/system_solvated.top

# Solution minimization
cd ../emin
gmx grompp -p ../top/system_solvated.top -c ../solvation/system_solvated.gro -f ../../chika_mdp/minimization.mdp -o minimization-sol.tpr
gmx mdrun -deffnm minimization-sol -v

if [ "$debug_level" = 2 ]; then
    echo "[Chikaterasu] Debug level 2 set. Exiting after solvation."
    exit 1
fi

: '
*************************************************************
Add ions

IONS NOT IMPLEMENTED YET
In our cases of K48-linked ubiquitin the protein itself was
neutral, so no opportunity yet to test ions yet.
*************************************************************
'

cd ../ions/

gmx grompp -f ../../chika_mdp/equilibration.mdp -p ../top/system_solvated.top -c ../solvation/system_solvated.gro -o genion.tpr -maxwarn 2
gmx genion -s genion.tpr -pname NA+ -nname CL- -neutral -o ions_added.gro

# Update names to avoid mismatch downstream
#awk '{if ($2 == "NA") $2 = "NA+"; print}' ions_added.gro > temp && mv temp ions_added.gro

awk '
NR <= 2 { print; next }

NF == 6 {
    resnum   = substr($0,  1, 5)
    resname  = substr($0,  6, 5)
    atomname = substr($0, 11, 5)
    atomnum  = substr($0, 16, 5)
    x        = substr($0, 21, 8)
    y        = substr($0, 29, 8)
    z        = substr($0, 37, 8)

    # If atom name is exactly " NA", replace with "NA+"
    if (atomname ~ /^ *NA$/) {
        atomname = " NA+"
    }

    printf "%s%s%s%s%s%s%s\n", resnum, resname, atomname, atomnum, x, y, z
    next
}

{ print }
' ions_added.gro >temp && mv temp ions_added.gro

#cd ../..

# Solution minimization again with ions
cd ../emin

# Now we need to count (and as of now manually add) how many ions we added
cp ../top/system_solvated.top ../top/system_solvated_ions.top
#exit 1

# Na: 0 --> 3
# W : N --> N-3
echo "\nNA+\t\t3" >>../top/system_solvated_ions.top

sed -i '3i#include "../../chika_mdp/martini_v2.0_ions.itp"' ../top/system_solvated_ions.top
awk '/^\s*W[ \t]+[0-9]+/ {$2 = $2 - 3} {print}' ../top/system_solvated_ions.top >temp && mv temp ../top/system_solvated_ions.top

gmx grompp -p ../top/system_solvated_ions.top -c ../ions/ions_added.gro -f ../../chika_mdp/minimization.mdp -o minimization-sol.tpr

gmx mdrun -deffnm minimization-sol -v

cd ../..

if [ "$debug_level" = 3 ]; then
    echo "[Chikaterasu-martini] Debug level 3 set. Exiting after adding counterions."
    exit 1
fi

: '
*************************************************************
Start the final equilibration and production MD loop
*************************************************************
'

for i in $(seq 1 $nruns); do
    # Equilibration
    cd runs/
    mkdir -p npt
    cd npt
    cp ../../gromacs/top/*.top .
    cp ../../gromacs/top/*.itp .
    cp ../../gromacs/emin/minimization-sol.gro .

    gmx grompp -p system_solvated.top -c minimization-sol -f ../../chika_mdp/equilibration.mdp -o equilibration.tpr -maxwarn 2
    gmx mdrun -deffnm equilibration -v

    if [ "$debug_level" = 4 ]; then
        echo "[Chikaterasu-martini] Debug level 4 set. Exiting after equilibration."
        exit 1
    fi

    # 5. Production
    cd ..
    mkdir -p md
    cd md

    cp ../npt/*.top .
    cp ../npt/*.itp .
    cp ../npt/equilibration.gro .
    cp ../npt/equilibration.tpr .

    gmx grompp -p system_solvated.top -f ../../chika_mdp/dynamic.mdp -o dynamic.tpr -c equilibration.tpr -maxwarn 2
    gmx mdrun -deffnm dynamic -nb gpu -v

    cd ../..

    # Move the data into their run directories
    mkdir -p runs/md_$i
    mkdir -p runs/nvt_$i
    mkdir -p runs/npt_$i

    mv runs/md/* runs/md_$i/
    mv runs/nvt/* runs/nvt_$i/
    mv runs/npt/* runs/npt_$i/

    mkdir -p runs/md
    mkdir -p runs/nvt
    mkdir -p runs/npt

    # End the run
    echo [Chikaterasu-martini] Run $i finished. Congratulations!

done
