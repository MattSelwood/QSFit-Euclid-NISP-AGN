import h5py
import numpy as np
from astropy.io import fits

# Input path variables
catalog_path = ""  # Path to CATALOG_WP9_INTERMEDIATE_RUN_v2.fits
output_path = "SimulatedSpectra/"

# Data
files = [
    "EUC_SIR_W-COMBSPEC_25463_20220128T115901.726349Z.h5",
    "EUC_SIR_W-COMBSPEC_25463_20220128T115935.423408Z.h5",
    "EUC_SIR_W-COMBSPEC_25463_20220128T115953.169464Z.h5"
]

# Read in CATALOG_WP9_INTERMEDIATE_RUN_v2.fits
catalog = fits.open(catalog_path)[1].data
source_ids = catalog['SOURCE_ID']


# Generate FITS for each simulated spectrum
for file in files:

    h5file = h5py.File(file, 'r')  # Read current h5 file

    # convert keys_list (SOURCE_IDs) to list, order it ascending and convert back to string to be called as labels
    source_list = list(map(str, sorted(list(map(int, h5file['Combined Spectra'].keys())))))
    
    # Cycle through each source in file
    for ID in source_list: 

        spectrum = h5file['Combined Spectra'][ID]

        # Construct wavelength array
        wstart = spectrum.attrs['Wlen start']
        wend = spectrum.attrs['Wlen end']
        winc = spectrum.attrs['Wlen inc']
        wavelength = np.arange(wstart, wend+10, winc, dtype=np.float32)

        # Get normalization factor
        norm_factor = spectrum['DATA']['COMBINED'].attrs['Normalization factor']

        # Get Flux and Uncertainty
        flux = np.array(spectrum['DATA']['COMBINED']['Science'])[0] * norm_factor
        sig = np.sqrt(np.array(spectrum['DATA']['COMBINED']['Variance'])[0]) * norm_factor

        # Transform to units of 10-17 erg/s/cm2
        flux_17 = flux / 1.E-17
        sig_17 = sig / 1.E-17

        # Create and write fits
        output_fits = fits.HDUList()
        col1 = fits.Column(name='Wavelength', unit='Angstrom', format='E', array=wavelength)
        col2 = fits.Column(name='Flux', unit='10^-17 erg/s/cm^2/A', format='E', array=flux_17)
        col3 = fits.Column(name='Uncertainty', unit='10^-17 erg/s/cm^2/A', format='E', array=sig_17)
        coldefs = fits.ColDefs([col1, col2, col3])
        bintable_hdu = fits.BinTableHDU.from_columns(coldefs)
        output_fits.append(bintable_hdu)
        output_fits.writeto(f"{output_path}{ID}.fits")


