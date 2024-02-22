for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	#target=`echo $bn | cut -d . -f 1`
	samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/MiaPaCa2.//g" | sed "s/.cpm//g"`
	
	#if [[ $bn =~ "day_1" ]]
	#then
	#	max=2.5
	#	min=0
	
	#else
	#	max=6
        #       min=0

	#fi

	plotHeatmap -m $f -o $bn.png --yMax 0.8 --yMin 0 --zMax 1.45 --colorMap RdBu_r --samplesLabel $samples
done
