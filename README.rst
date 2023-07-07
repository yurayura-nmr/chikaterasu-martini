Chikaterasu Martini version
===========================

.. image:: logo.png
   :alt: Chikaterasu logo
   :align: right

Bash script to automate setup of multiple identical MD simulations.
MD parameters are for the martini 2.2 forcefield.

Usage
-----

1. Prepare the topology using the martinize python2.7 script (see folder martinize-script).
2. Prepare the coarse-grained protein structure file in ./gromacs/coord/
3. Run: ./chikaterasu-martini 1 [This will test if topology can be generated correctly in gromacs]
4. Run: ./chikaterasu-martini 2 [This will test if solvation worked fine]
5. Run: ./chikaterasu-martini 3 [This will test if counterions were correctly added (not yet implemented)]
6. Run: ./chikaterasu-martini 4 [This will test if equilibration worked correctly in gromacs]
7. Run: ./chikaterasu-martini 5 [This will test if a short equilibration MD trajectory can be obtained in gromacs]
8. Run: ./chikaterasu-martini 0 [This will assume everything so far went well and go to the actual production MD run]

MD parameters are adjusted by editing the gromacs .mdp files in the chika_mdp directory.

Misc
----

Change log
""""""""""

See github commits

To do
"""""

See issues
