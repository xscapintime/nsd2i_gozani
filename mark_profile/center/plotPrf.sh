#for f in *.gz
for f in *H3*.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	target=`echo $bn | cut -d . -f 1,2`
	mark=`echo $bn | cut -d . -f 2`
	
	#if [[ $bn =~ "day_1" ]]
	#then
	#	max=2.5
	#	min=0
	
	#else
	#	max=6
        #       min=0

	#fi

	plotHeatmap -m $f -o $target.png --colorMap RdBu_r --samplesLabel "$target.rp1" "$target.rp2"  --xAxisLabel "NicE-seq IDR" --regionsLabel "$mark CUT&RUN" 
done
