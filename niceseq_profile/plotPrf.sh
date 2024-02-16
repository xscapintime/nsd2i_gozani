for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	target=`echo $bn | cut -d . -f 1,2,3`
	
	#if [[ $bn =~ "day_1" ]]
	#then
	#	max=2.5
	#	min=0
	
	#else
	#	max=6
        #       min=0

	#fi

	plotHeatmap -m $f -o $bn.png --regionsLabel "NicE-seq" --yMax 1 --yMin 0 --colorMap RdBu_r --zMax 1.6 --kmeans 3
done
