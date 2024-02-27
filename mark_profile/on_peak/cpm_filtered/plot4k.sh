for f in *.4kdist.10bin.gz

do

	target=`echo $f | cut -d . -f 1,2,3`
	samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/MiaPaCa2.//g" | sed "s/.cpm//g" | sed "s/.sorted//g"`

	max=0.15
	min=0.04
	hmax=0.6

	echo plot $target
	plotHeatmap -m $f -o $target.png --colorMap RdBu_r --samplesLabel $samples --regionsLabel "CUT&RUN" --xAxisLabel "NicE-seq peak" --yMin $min --yMax $max --zMax $hmax

done

