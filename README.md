# Chikaterasu Martini Version

![Chikaterasu Logo](logo.png)

> Bash script to automate setup of multiple identical MD simulations.  
> Compatible with the **Martini 2.2** force field.

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
python2.7 martinize.py -f Abeta.pdb -o Abeta.top -x Abeta-CG.pdb -p backbone
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

## 🚀 Simulation Workflow

Run the automation script in stages:

```bash
./chikaterasu-martini 1
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
