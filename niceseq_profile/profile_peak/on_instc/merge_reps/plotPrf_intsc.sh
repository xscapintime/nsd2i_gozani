for f in *.intsc.2kdist.10bin.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	#target=`echo $bn | cut -d . -f 1`
	samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/MiaPaCa2.//g" | sed "s/.rp1.cpm//g"`
	#regions=`gunzip -c $f | head -n 1 | cut -d "[" -f 9 |  cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/.bed//g"`
		
	#if [[ $bn =~ "day_1" ]]
	#then
	#	max=2.5
	#	min=0
	
	#else
	#	max=6
        #       min=0

	#fi

	#plotHeatmap -m $f -o $bn.png --yMax 0.7 --yMin 0 --zMax 1.1 --colorMap RdBu_r --samplesLabel $samples --regionsLabel NicE-seq
	
	plotProfile -m $f -o $bn.profile.png --numPlotsPerRow 3 --samplesLabel $samples --refPointLabel "NicE-seq peak" --plotWidth 5.5 \
--averageType mean --plotType se --regionsLabel NiCE-seq  --colors  "#5773CC" "#5773CC" "#5773CC" "#8B8989" "#8B8989" "#8B8989" \
--legendLocation upper-left
done
