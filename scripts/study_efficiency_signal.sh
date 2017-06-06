
paramfile="parameters.txt"

cd ~/files/
#list_of_runs=`ls -d em*14per`
list_of_runs=`ls -d em*40per`
cd -


for file in ${list_of_runs}; do

    echo RUNNING FILE ${file}

    E=`echo ${file} | awk -F "_" '{print $2}' | awk -F "MeV" '{print $1}' `
    coverage=`echo ${file} | awk -F "_" '{print $6}' | awk -F "per" '{print $1}' `

    ln -s ~/files/${file}/*txt .


### thresholds for 14%
#    thresholds='4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30'

### thresholds for 40%
    thresholds='10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50'


    for T in $thresholds; do
	
	echo RUNNING FILE ${file} THRESHOLD ${T}
	
	sed -i "s/.*water_like_threshold_number_of_pmts.*/ water_like_threshold_number_of_pmts ${T}/g" ${paramfile}
	
	./daq_code
	
	mv all_hits_emerald*.txt ~/files/${file}/
	
    done

    rm all_hits_*txt
    rm all_pmts.txt
    rm detector.txt

done
