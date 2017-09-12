#  README for Optipol Reduction  #
### 09/12/2017

### REQUIREMENTS
 * astropy
 * [ccdproc](http://ccdproc.readthedocs.io/en/latest/ccdproc/install.html "Astropy ccdproc")

##Example usage

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
