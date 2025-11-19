# Chikaterasu Martini Version

**C**ommand-line **H**elper for **I**nstalling **K**ernel-**A**ccelerated **T**hermodynamic **E**xploration **R**uns **A**nd **S**etting **U**p

<img src="logo.png" alt="Chikaterasu logo" align="right" />

> Bash script to automate setup of multiple identical MD simulations.  
> Compatible with the **Martini 2.2** force field.

**Prerequisite**: Install Anaconda or Miniconda for Python environment management before proceeding.

---

## 📦 Setup Instructions

### 1. Generate the Topology

Use the `martinize` script with Python 2.7. You can set up a dedicated conda environment for this:

```bash
conda create -n py27 python=2.7
conda activate py27
```

Run the `martinize.py` script:

```bash
cd martinize-script
python2.7 martinize.py -f Abeta.pdb -o Abeta.top -x Abeta-CG.pdb -p backbone -elastic -ef 500 -el 0.5 -eu 0.9 -ea 0 -ep 0
```

---

### 2. Prepare GROMACS Input Files

Move the generated files into the expected GROMACS directories:

```bash
cp Abeta-CG.pdb ../gromacs/coord/
cp Abeta.top ..
cp Protein_A.itp ..
```

---

### 3. Configuration

#### Edit Script Variables
Open `chikaterasu-martini.sh` and set the following variables at the top of the script:
- **Protein name**: Set to your PDB filename (without extension)
- **Topology name**: Set to your generated .top filename

#### Determine Box Size
For optimal solvation, determine appropriate box dimensions:
1. Use the all-atom Chikaterasu script temporarily **only for the solvation step**
2. Run: `./chikaterasu.sh 2` (stops after solvation)
3. Check the generated `*_newbox.gro` file in `gromacs/solvation/`
4. Note the box dimensions listed at the bottom of the file
5. Apply these same dimensions to the Martini version simulation parameters


---

## 🚀 Simulation Workflow

Run the automation script in stages:

```bash
./chikaterasu-martini.sh 1
```

- ✅ **Step 1**: Validates topology generation in GROMACS

```bash
./chikaterasu-martini 2
```

- ✅ **Step 2**: Tests solvation

```bash
./chikaterasu-martini 3
```

- 🚧 **Step 3**: Checks counterion addition (not yet implemented)

```bash
./chikaterasu-martini 4
```

- ✅ **Step 4**: Runs a short equilibration trajectory

```bash
./chikaterasu-martini 0
```

- ✅ **Step 0**: Launches the full production MD run, assuming prior steps succeeded

---

## ⚙️ MD Parameters

Simulation parameters can be customized by editing the `.mdp` files in the `chika_mdp/` directory.

---

## 📚 Miscellaneous

### 📝 Change Log

See [GitHub commit history](#) for detailed changes.

### 📌 To-Do

Refer to the [GitHub Issues](#) page for outstanding tasks and feature requests.

---

Let me know if you'd like this formatted into a `README.md`, have links to insert, or want badges added (like license, Python version, etc.)!
