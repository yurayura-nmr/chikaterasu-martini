# For example, if the CG structure to be back-calculated is called "1p30_nm_cyclic-CG_to_backcalc.gro" 

mkdir -p final_results

./chika_initram.sh -f 1p30_nm_cyclic-CG_to_backcalc.gro -o aa_charmm.gro -to charmm36 -p ./topol.top 
mv aa_charmm.gro final_results/1p30_nm_backcalculated.gro

# can add more back calculation runs here
