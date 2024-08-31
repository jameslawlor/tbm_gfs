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
* Other features
  * Dockerised
  * CI pipeline using Github Actions

## TODO list
See [project board](https://github.com/users/jameslawlor/projects/1)
