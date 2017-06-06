
paramfile="parameters.txt"
outfile="processing_times.txt"

rm ${outfile}

sed -i "s/.*output_txt.*/ output_txt 0/g" ${paramfile}

make

times='500ns 1mus 10mus 100mus 1ms 10ms'
coverages='14 40'
distances='500. 1000. 2000.'

for c in ${coverages}; do

    cp hits/all_pmts_${c}per.txt all_pmts.txt

    for t in ${times}; do

	cp hits/all_hits_${c}per_${t}.txt all_hits_1.txt

	for d in ${distances}; do

	    sed -i "s/.*distance_between_vertices.*/ distance_between_vertices ${d} /g" ${paramfile}

	    echo DOING COVERAGE ${c} TIME ${t} DISTANCE ${d}

	    ./daq_code
	    T=`./daq_code | grep "total execution" | awk '{print $5}'`

	    echo "C ${c} T ${t} D ${d} --> T ${T}" >> ${outfile}

	done

    done

done
