
paramfile="parameters.txt"

wallikedistance='0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.'

for D in $wallikedistance; do

    echo RUNNING DISTANCE ${D}

    sed -i "s/.*wall_like_distance.*/ wall_like_distance ${D}/g" ${paramfile}

    ./daq_code

    mv all_hits_emerald*.txt outputs/all_hits_emerald_distance_${D}.txt

done
