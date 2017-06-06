
paramfile="parameters.txt"

#coverage='14'
coverage='40'

dir="large_radioactive_${coverage}per"

ln -s ~/files/${dir}/*txt .

### thresholds for 14%
#thresholds='4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30'

### thresholds for 40%
thresholds='10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50'

for T in $thresholds; do

    echo RUNNING THRESHOLD ${T}

    sed -i "s/.*water_like_threshold_number_of_pmts.*/ water_like_threshold_number_of_pmts ${T}/g" ${paramfile}

    ./daq_code

    mv all_hits_emerald*.txt ~/files/${dir}/
    
done

rm all_hits_*txt
rm all_pmts*txt
rm detector.txt

