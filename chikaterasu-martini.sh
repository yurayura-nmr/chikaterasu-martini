#!/bin/sh

: '
*************************************************************
Erik Walinda
Kyoto University
Graduate School of Medicine

chikaterasu-martini
tested for martini 2.2 in gromacs 2022 or later
*************************************************************
'

: '
*************************************************************
Manually setup parameters for this run
Most parameters are defined by mdp files

Debug level 

0   full production MD run [NPT]; script will not stop at intermediate steps.
1   topology generation [pdb2gmx]; script will stop after topology made.
2   solvation; script will stop after box solvation.
3   addition of counterions; script will stop after ions added.
4   equilibration; script will stop before production MD. 

Requires the following files to be prepared for the system of interest. For example in the case of K48-linked diubiquitin:

* ./k48.top - topology file for K48-linked diubiquitin [it will #include another itp file]
* ./Protein_A+Protein_B.itp - topology with the two ubiquitin chains [created by martinize script]
* ./chika_mdp/martini.itp - martini forcefield parameters; potentially modified to modify protein-water interactions etc.
* ./gromacs/cord/k48-cg.pdb - the starting conformation PDB file should be placed into gromacs/coord
*************************************************************
'

# == How many independent trajectories shall be simulated ? ==
nruns=10                    # 1 ... 10

# == Debug level ==
debug_level=1               # Manually set debug level. Or give as argument, e.g.: ./chikaterasu 0

# === Box shape and size ===
box_manual=true
box_dim="7.500   7.500   7.500"

# top/itp must be manually prepared and ready in folder
protein_name="k48-CG"
top="k48" # name


: '
*************************************************************
After setting up the above parameters, 
no manual change should be necessary in the script from here
*************************************************************
'
if [ -z "$1" ]
then
    read -p "[Chikaterasu-martini] Command line arguments are empty. Using manually set debug level $debug_level." dummy
else
    read -p "[Chikaterasu-martini] Command line arguments provided. Using first argument as debug level $1." dummy
    debug_level=$1
fi

# Setup directories for the run
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

: '
*************************************************************
Assuming that the topology has already been correctly set up using the martinize.py script.
If not, execute the script before running chikaterasu-martini.

For example for K48-linked diubiquitin (need to merge 2 protein chains into one topology):
python2.7 martinize.py -f 1aar_modified.pdb -o 1aar_modified.top -x 1aar_modified-CG.pdb -dssp dssp -p backbone -merge A,B

Even simpler for monoubiquitin
python2.7 martinize.py -f 1UBQ.pdb -o single-ubq.top -x 1UBQ-CG.pdb -dssp dssp -p backbone
*************************************************************
'

: '
*************************************************************
[editconf]
Generate box around protein and perform 1 vacuum minimization
*************************************************************
'
# Create box with editconf
cd gromacs/coord
gmx editconf -f $protein_name.pdb -d 1.0 -bt triclinic -o $protein_name.gro

if [ "$box_manual" = true ] ; then
    sed -i '$ d' $protein_name.gro
    echo "$box_dim" >> $protein_name.gro
fi

# Use topology to prepare vacuum minimization
cd ../top/
cp ../../chika_mdp/martini.itp .
cp ../../*.top .
cp ../../*.itp .

gmx grompp -p $top.top -f ../../chika_mdp/minimization.mdp -c ../coord/$protein_name.gro -o ../emin/minimization-vac.tpr

# Vacuum minimization
cd ../emin/
gmx mdrun -deffnm minimization-vac -v

if [ "$debug_level" = 1 ] ; then
    echo "[Chikaterasu-martini] Debug level 1 set. Exiting after initial topology generation [pdb2gmx]"
    exit 1
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
echo -n "W\t\t" >> ../top/system_solvated.top
grep -c W system_solvated.gro >> ../top/system_solvated.top

# Solution minimization
cd ../emin
gmx grompp -p ../top/system_solvated.top -c ../solvation/system_solvated.gro -f ../../chika_mdp/minimization.mdp -o minimization-sol.tpr
gmx mdrun -deffnm minimization-sol -v

if [ "$debug_level" = 2 ] ; then
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

if [ "$debug_level" = 3 ] ; then
    echo "[Chikaterasu-martini] Debug level 3 set. Exiting after adding counterions."
    exit 1
fi

cd ../..

: '
*************************************************************
Start the final equilibration and production MD loop
*************************************************************
'

for i in `seq 1 $nruns`;

do
    # Equilibration
    cd runs/
    mkdir -p npt
    cd npt
    cp ../../gromacs/top/*.top .
    cp ../../gromacs/top/*.itp .
    cp ../../gromacs/emin/minimization-sol.gro .

    gmx grompp -p system_solvated.top -c minimization-sol -f ../../chika_mdp/equilibration.mdp -o equilibration.tpr -maxwarn 2
    gmx mdrun -deffnm equilibration -v 

    if [ "$debug_level" = 4 ] ; then
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
