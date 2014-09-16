#  README for Optipol Reduction  #
### 09/16/2014

## REQUIREMENTS
-numpy
-pyfits
-matplotlib
-pyregion
-PyGuide
-astropy
-LAcosmics (included)


## BIAS SUBTRACT
```bash
python zeroproc.py ../raw/bias* -zout ../reduced/MasterBias.fit
python zeroproc.py ../reduced/MasterBias.fit --data ../raw/NGC3628*.fit -dout ../reduced
python zeroproc.py ../reduced/MasterBias.fit --data ../raw/flat*.fit -dout ../reduced
```

## FLAT FIELD
```bash
python flatproc.py ../reduced/flat* -fout ../reduced
python flatproc.py ../reduced/MasterFlat*.fit --data ../reduced/NGC3628*.fit -dout ../reduced
```

## NORMALIZE
```bash
python expnorm.py ../reduced/NGC3628*.fit
```

## SKY SUB
**_Make reg boxes in ds9, one top and one bottom_**
python skysub.py ../reduced/NGC3628*.fit -reg skyBoxes.reg

## WOLLY SPLIT
python wolly_split.py ../reduced/NGC3628*.fit -o ../reduced/ -prefix NGC3628_R-
python cosmic_clean.py ../reduced/NGC3628*_A.fit
python cosmic_clean.py ../reduced/NGC3628*_B.fit

## ALIGN
python point_align.py ../reduced/NGC3628*_A.fit -coords x y
python point_align.py ../reduced/NGC3628*_B.fit -coords x y
python point_align.py ../reduced/NGC3628*.al.fit --aligned


