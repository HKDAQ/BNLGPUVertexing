
paramfile="parameters.txt"
paramfileSK="parameters_SK.txt"
paramfileHK="parameters_HK.txt"

rm ${paramfile}
ln -s ${paramfileSK} ${paramfile}

dir="SK_only_dark_noise"

ln -s ~/files/${dir}/*txt .


### thresholds for SK
#thresholds='50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100'

./daq_nhits
    
mv all_hits_emerald*.txt ~/files/${dir}/
    
rm all_hits_*txt
rm all_pmts.txt
rm detector.txt

rm ${paramfile}
ln -s ${paramfileHK} ${paramfile}
