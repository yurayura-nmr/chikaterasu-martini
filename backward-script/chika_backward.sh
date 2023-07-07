"""
Erik Walinda
Kyoto University
Graduate School of Medicine

Required input files
--------------------

* frame-of-interest_to_backcalc.pdb     # the CG structure to be backcalculated [frame of interest out of trajectory]
* protein-name_atomistic.pdb            # Original atomistic coordinates (non-coarse-grained)
* dynamic.tpr                           # additionally needed

Preparation
-----------

1. Run ana_chikaterasu-martini.sh first to obtain frames of the trajectory with PBC removed.
   e.g. edit ana_chikaterasu-martini.sh with the following settings:
   
   traj=true
   dt=5000
   runs=1

   and run it:
   
   ./ana_chikaterasu-martini.sh

   ... which wil yield a PBC-removed coarse-grained trajectory file.
   
2. Extract the frame of interest (as identified by Pymol or distance analysis). 
   Here we are backmapping a single frame, not an entire trajectory.
   So let's assume as an example the frame is at time t=3180 ns.
   
   gmx trjconv -f ./md_fit.xtc -s ./md_target.tpr -b 3180 -e 3180 -tu ns -o 0p55_nm_cyclic.pdb -conect

2. Run the commands outlined below step by step.

3. At initram.sh calculation will start. GROMACS will use the GPU if possible.

In this repository, we adjusted initram.sh for non-ancient gromacs versions since the original initram file writes its own mdp files which will not work in gromacs-2022.

Expected notes from GROMACS:

NOTE 1 [file unknown]:
  You are using constraints on all bonds, whereas the forcefield has been
  parametrized only with constraints involving hydrogen atoms. We suggest
  using constraints = h-bonds instead, this will also improve performance.

NOTE 2 [file 4-mdpr-0.0005.mdp]:
  Removing center of mass motion in the presence of position restraints
  might cause artifacts. When you are using position restraints to
  equilibrate a macro-molecule, the artifacts are usually negligible.

NOTE 3 [file 4-mdpr-0.0005.mdp]:
  You are using a plain Coulomb cut-off, which might produce artifacts.
  You might want to consider using PME electrostatics.
"""

# Run these commands step-by-step (not meant to be executed as a single script!)

# 1. Create gro file of CG model using trjconv.
#    This will use the frame of interest that we extracted from the trajectory in the above explanation.

#               (frame-of-interest to be backmapped)
#              ____________ CG structure ____________
gmx trjconv -f cyclic_extended_frame1_to_backcalc.pdb -s dynamic.tpr -o cyclic_extended_frame1_to_backcalc-CG.gro 

# (1b). Most cases probably can skip this step.
#       In the special case of diubiquitin, LYQ and GLQ were used as isopeptide-linked amino acid names.
#       These need to changed back before backmapping.
#       Convert GLQ back to GLY and LYQ back to LYS.

# vi ... # e.g. k48_closed_atomistic.pdb
# :%s/LYQ/LYS/g
# :%s/GLQ/GLY/g

# 2. Obtain the atomistic topology under the CHARMM27 forcefield for the system of interest.
#    i.e., the top file, not the orginal pdb file.
#    This can be easily done using gmx2pdb or using chikaterasu (the non-martini version) as follows:
#    Run "./chikaterasu 1" to obtain atomistic topology (CHARMM27 forcefield).
./chikaterasu.sh 1

# and copy the resulting atomistic topology to the working folder of the backward script.
cp chikaterasu/gromacs/top/topol.top chikaterasu-martini/backward-script/

# 3. We also need the position restraint file (used by backward).
#    This is also easily obtained by chikaterasu (debug level 5).
./chikaterasu.sh 5
cp chikaterasu/gromacs/top/posre.itp chikaterasu-martini/backward-script/

# 4. Try to run initram:
#    This requires python 2.7, so need to activate or set up an environment:
#    conda create --name py27 python=2.7

conda activate py27
./chika_initram.sh -f frame-of-interest_to_backcalc-CG.gro -o aa_charmm.gro -to charmm36 -p ./topol.top

# If all runs well, we can then visualize the atomistic structure (in gro format9 as usual (vmd or pymol).
