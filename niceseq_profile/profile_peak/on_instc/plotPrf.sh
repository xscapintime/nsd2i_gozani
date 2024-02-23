for f in *.intsc.unq*.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	#target=`echo $bn | cut -d . -f 1`
	samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/MiaPaCa2.//g" | sed "s/.cpm//g"`
	regions=`gunzip -c $f | head -n 1 | cut -d "[" -f 9 |  cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/.bed//g"`
		
	#if [[ $bn =~ "day_1" ]]
	#then
	#	max=2.5
	#	min=0
	
	#else
	#	max=6
        #       min=0

	#fi

	plotHeatmap -m $f -o $bn.png --yMax 0.75 --yMin 0 --zMax 0.9 --colorMap RdBu_r --samplesLabel $samples --regionsLabel $regions
done
