# BIAS SUBTRACT
zeroproc.py ../raw/bias* -zout ../reduced/MasterBias.fit
zeroproc.py ../reduced/MasterBias.fit --data ../raw/NGC3628*.fit -dout ../reduced
zeroproc.py ../reduced/MasterBias.fit --data ../raw/flat*.fit -dout ../reduced

# FLAT FIELD
python flatproc.py ../reduced/flat* -fout ../reduced
python flatproc.py ../reduced/MasterFlat*.fit --data ../reduced/NGC3628*.fit -dout ../reduced

# NORMALIZE
python expnorm.py ../reduced/NGC3628*.fit

# SKY SUB
## Make reg boxes in ds9, one top and one bottom
python skysub.py ../reduced/NGC3628*.fit -reg skyBoxes.reg

# WOLLY SPLIT
python wolly_split.py ../reduced/NGC3628*.fit -o ../reduced/ -prefix NGC3628_R-
python cosmic_clean.py ../reduced/NGC3628*_A.fit
python cosmic_clean.py ../reduced/NGC3628*_B.fit


# ALIGN THE SHIT
###python get_regions.py -o centersA3.tsv
python bullshitalign.py centersA4.tsv
#ds9 ../reduced/*_B.fit
##python get_regions.py -o centersB.tsv
python bullshitalign.py centersB2.tsv

# REBIN
rebin.py ../reduced/*_s.fit -size 10 --func sum

# QU PAIRS
## ds9 
python QUpair.py ../reduced/*_s_bin.fit -reg dustBox.reg

#python profalign.py ../reduced/*_A.fit 470 330 200



#python shiftalign.py ../reduced/*_A.fit 492 432
