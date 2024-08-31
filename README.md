# tbm_gfs

Project for code related to my PhD work 2012-2016 at Trinity College Dublin.

Google Scholar: https://scholar.google.fi/citations?user=XeiCaFIAAAAJ&hl=en&oi=sra

## Setup
```
python3 -m venv .env
source .env/bin/activate
pip install -r requirements.txt
pip install -e .
```

## Currently implemented
* Bulk lattice Green's Functions
  * Graphene
    * Double integral form
    * Single integral forms (ka, kz)
    * Stationary phase approximation (armchair direction, same sublattice only)
  * Carbon nanotubes
    * Integral form 
    * analytic (exact) form 
* Examples
    * Calculating local density of states for graphene

## TODO list
* Code
  * Stationary Phase Approximation (Zigzag, opposite sublattices)
  * GFs for Hamiltonians with 2nd NN hopping integrals and wavefunction overlap 
  * Bulk GFs for semi-infinite graphene sheets
  * Friedel oscillations 
  * Recursive GFs for nanoribbon and nanotubes
    * Standard recursive methods
    * Rubio-Sancho / Sancho-Rubio method for recursive GFs
    * Kubo formula for conductance
  * Magnetic coupling / RKKY
  * Impurities, dN, dE calculations
  * Sublattice asymmetry (Ising model?)
  * Graphene "sector" method
* Features
  * Docker images
  * Docker build/run commands
