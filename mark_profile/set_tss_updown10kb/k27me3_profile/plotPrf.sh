for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	# echo $bn
	
	regions=`gunzip -c $f | head -n 1 | cut -d "[" -f 9 |  cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/.bed//g" | sed "s/_prom_1k//g" | cut -d " " -f 1`

	if [[ $bn =~ "d5"  ]];then

		samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' |  sed "s/,/ /g" | sed "s/_merged//g" | sed "s/.ms_cpm//g"`

	else
		samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' |  sed "s/,/ /g" | sed "s/.MSnorm//g" | sed "s/_H/.H/g"`

	fi


	if [[ $bn =~ "EPI" ]];then
		ymin=0
		ymax=14.3
	elif [[ $bn =~ "KRAS" ]];then
		ymin=0
		ymax=14
	elif [[ $bn =~ "pan" ]];then
		ymin=0
		ymax=12
	elif [[ $bn =~ "PRC" ]];then
		ymin=0
		ymax=22
	fi

	echo $bn $ymin $ymax

	plotProfile -m $f -o $bn.prof.pdf --refPointLabel 'TSS' \
	--regionsLabel "Gene set" "Controls" \
	--colors "#4682B4" "#CDC9C9" \
	--plotWidth 7 --plotTitle $regions \
	--samplesLabel $samples \
	--yMin $ymin --yMax $ymax

done
