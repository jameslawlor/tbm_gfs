# tbm_gfs

Project for code related to my PhD work 2012-2016 at Trinity College Dublin.

Google Scholar: https://scholar.google.fi/citations?user=XeiCaFIAAAAJ&hl=en&oi=sra

## Setup
```
python3 -m virtualenv .env
source .env/bin/activate
pip install -r requirements.txt
pip install -e .
```

## Currently implemented
* Green's Functions
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
  * GFs for Hamiltonian with 2nd NN overlap 
  * Friedel oscillations 
  * Recursive GFs for nanoribbon and nanotubes using Dyson equation
  * Magnetic coupling / RKKY
  * Impurities, dN, dE calculations
  * Sublattice asymmetry (Ising model?)
* Features
  * Docker images
  * Docker build/run commands
