for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	target=`echo $bn | cut -d . -f 1,2`
	
	#if [[ $bn =~ "day_1" ]]
	#then
	#	max=2.5
	#	min=0
	
	#else
	#	max=6
        #       min=0

	#fi

	plotHeatmap -m $f -o $bn.png --yMax 1 --yMin 0 --colorMap RdBu_r  --regionsLabel promoter enhancer
done
