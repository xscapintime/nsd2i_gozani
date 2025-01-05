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


	plotProfile -m $f -o $bn.prof.pdf --refPointLabel 'TSS' \
	--samplesLabel "NSD2i" "Control" \
	--colors "#66B8D8" "#9999A1" \
	--plotWidth 8.5	--perGroup \
	--regionsLabel "H3K27ac Day $day" \
	--yMin 0.053 --yMax 0.248

done
