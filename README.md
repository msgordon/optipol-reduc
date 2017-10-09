#  README for Optipol Reduction  #
### 09/12/2017

### REQUIREMENTS
 * astropy
 * [ccdproc](http://ccdproc.readthedocs.io/en/latest/ccdproc/install.html "Astropy ccdproc")
 * [pyregion](http://pyregion.readthedocs.io/en/latest/ "pyregion")
 * [photutils](https://photutils.readthedocs.io/en/latest/index.html "photutils")
 * [image-registration](https://github.com/keflavich/image_registration "image registration)
 * [joblib](http://pythonhosted.org/joblib/index.html "joblib")

## General notes
All scripts defined below use `argparse` for processing command line arguments.  Each program can be run with the `-h` flag to print a description and example usage (e.g. `python stack.py -h`).

Each script also takes the optional `--c` flag to clobber/overwrite files on write.

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
Optionally, add the `--n` flag to normalize by EXPTIME.  This will be performed regardless in `proc.py'.

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

### CCD Process
```bash
python proc.py NGC4565/reduced/NGC4565_R-0*.fit -bias NGC4565/reduced/MasterBias.fits -dark NGC4565/reduced/MasterDark_120.fits -flat NGC4565/reduced/MasterFlat_*.fits -odir NGC4565/proc/ -maskfile masks/bad_pix.fits
```
Specifying `-maskfile` attempts to mask bad pixels during the CCD reduction process.  However, this does not always produce clean images during the dark subtraction step.  Adding the `--fixpix` flag will interpolate over the bad pixels, yielding nicer looking images, but this adds **significant** computation time.  We recommend using the `-njobs` option to process multiple images in parallel.  Avoid using every CPU, since the re-gridding will fill up each CPU.


### Split Wollaston beams
```bash
python split.py NGC4565/proc/*.fit -odir NGC4565/split -maskfile masks/wolly_mask.reg
```
Using the `wolly_mask.reg` ds9 region file, split the images into ordinary (A) and extraordinary (B) beams.
