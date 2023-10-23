# QSFit Euclid NISP AGN 

Repo for easy reproducibility of the QSFit analyis of 'Incident' and 'Simulated' AGN spectra performed in Euclid Collaboration: Lusso et al. (in-prep).


## Code repository

The repo shall contain the following:
  - A Python script (`extract_simulated_fits_from_h5.py`) to extract simulated spectra from the H5 files;
  - A Julia package to read incident and simulated spectra, to run the QSFit analysis, and to generate the final catalogs.


### Run QSFit analysis


#### Clone this repository
```
git clone https://github.com/MattSelwood/QSFit-Euclid-NISP-AGN.git
```
Note: the above command should be executed **only once** to download the necessary files.

#### Install Julia packages

Change into the repository directory and start a Julia session with the commands:
```
cd QSFit-Euclid-NISP-AGN
julia --project=.
```

Install all required packges with:
```julia
using Pkg
Pkg.instantiate()
```
Note: the above commands should be executed **only once** to install the necessary Julia files.

#### Run the analysis

Change into the repository directory and start a Julia session with the commands:
```
cd QSFit-Euclid-NISP-AGN
julia --project=.
```

Run the analysis with:
```julia
using WP9_QSFitAnalysis
run_analysis("Input/CATALOG_WP9_INTERMEDIATE_RUN_v2.fits", # Catalog
             "Input/SPECTRA_WP9_INTERMEDIATE_RUN_v2.fits", # Incident spectra
             "Input/SimulatedSpectra",                     # Simulated spectra
             "output")                                     # Output directory
```

Extract the catalogs with:
```julia
using WP9_QSFitAnalysis
extract_catalogs("output")  # Must be the same directory used for output in the previous step
```

