#  README for Optipol Reduction  #
### 09/12/2017

### REQUIREMENTS
 * astropy
 * [ccdproc](http://ccdproc.readthedocs.io/en/latest/ccdproc/install.html "Astropy ccdproc")
 * [pyregion](http://pyregion.readthedocs.io/en/latest/ "pyregion")
 * [joblib](http://pythonhosted.org/joblib/index.html "joblib")

## Example usage

### Construct master bias
```bash
python stack.py NGC4565/raw/*.fit -imgtype bias -o NGC4565/reduced/MasterBias.fits
```

### Construct master darks
```bash
python stack.py NGC4565/raw/*.fit -imgtype dark -exptime 60 -o NGC4565/reduced/MasterDark_60.fits
python stack.py NGC4565/raw/*.fit -imgtype dark -exptime 120 -o NGC4565/reduced/MasterDark_120.fits
```
Optionally, add the `--n` flag to normalize by EXPTIME.

### Construct master flats by HWP position
```bash
python stack.py NGC4565/raw/flat*.fit -o NGC4565/reduced/MasterFlat.fits --normw --wolly -maskfile masks/wolly_mask.reg
```
Optionally, add `-njobs -1` to parallelize.

### Clean images
The array has many bad pixels.  Use `LACosmics` to clean.
```bash
python clean.py NGC4565/raw/NGC4565_R-0* -odir NGC4565/reduced/ -sclip 3
```
