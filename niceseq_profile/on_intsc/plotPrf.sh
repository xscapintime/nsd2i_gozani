for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	target=`echo $bn | cut -d . -f 1`
	
	#if [[ $bn =~ "NSD2i" ]]
	#then
	#	group=NSD2i

	
	#else
	#	group=vehicle

	#fi

	plotHeatmap -m $f -o $target.nice.png --regionsLabel "NicE-seq" --yMax 1 --yMin 0 --colorMap RdBu_r --zMax 1.4 --samplesLabel "$target.D1.rp1" "$target.D1.rp2" "$target.D5.rp1" "$target.D5.rp2" "$target.D9.rp1" "$target.D9.rp2" --xAxisLabel "NicE-seq IDR"
done
