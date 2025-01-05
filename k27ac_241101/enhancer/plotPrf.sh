for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	echo $bn
	day=`echo $bn | cut -d "." -f 1`
	

	#if [[ $bn =~ "d5"  ]];then

	# samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' |  sed "s/,/ /g" | sed "s/_merged//g" | sed "s/.ms_cpm//g"`

	#else
		#samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' |  sed "s/,/ /g" | sed "s/.MSnorm//g" | sed "s/_H/.H/g"`

	#fi


	plotProfile -m $f -o $bn.prof.pdf --refPointLabel 'Enhancer' \
	--samplesLabel "NSD2i rep1" "NSD2i rep2" "Control rep1" "Control rep2" \
	--colors "#66B8D8" "#C8D746" "#9999A1" "#CDC9C9" \
	--plotWidth 8.5	--perGroup \
	--regionsLabel "H3K27ac Day $day" \
	--yMin 0.04175 --yMax 0.057

done
