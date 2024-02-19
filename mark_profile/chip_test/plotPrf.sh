#for f in *.gz
for f in `ls *H3*.gz | grep -v _D`
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

	plotHeatmap -m $f -o $target.png --colorMap RdBu_r --samplesLabel "$target.D1.rp1" "$target.D1.rp2" "$target.D5.rp1" "$target.D5.rp2" "$target.D9.rp1" "$target.D9.rp2" \
	 --xAxisLabel "NicE-seq IDR" --regionsLabel "$mark CUT&RUN" 
done
