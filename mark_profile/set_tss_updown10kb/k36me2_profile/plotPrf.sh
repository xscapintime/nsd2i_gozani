for f in p*.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	# echo $bn
	
	regions=`gunzip -c $f | head -n 1 | cut -d "[" -f 9 |  cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/.bed//g" | sed "s/_prom_1k//g" | cut -d " " -f 1`

	#if [[ $bn =~ "d5"  ]];then

		#samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' |  sed "s/,/ /g" | sed "s/_merged//g" | sed "s/.ms_cpm//g"`

	#else
		samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' |  sed "s/,/ /g" | sed "s/.MSnorm//g" | sed "s/_H/.H/g"`

	#fi


	if [[ $bn =~ "EPI" ]];then
		ymin=3
		ymax=53
	elif [[ $bn =~ "KRAS" ]];then
		ymin=3
		ymax=53
	elif [[ $bn =~ "pan" ]];then
		ymin=3
		ymax=51
	elif [[ $bn =~ "PRC" ]];then
		ymin=1
		ymax=46
	fi

	echo $bn $ymin $ymax

	plotProfile -m $f -o $bn.prof.pdf --refPointLabel 'TSS' \
	--regionsLabel "Gene set" "Controls" \
	--colors "#B3EE3A" "#CDC9C9" \
	--plotWidth 7 --plotTitle $regions \
	--samplesLabel $samples \
	--yMin $ymin --yMax $ymax

done
