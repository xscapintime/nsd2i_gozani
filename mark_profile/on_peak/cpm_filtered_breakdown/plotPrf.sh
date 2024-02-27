for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	target=`echo $bn | cut -d . -f 1,2`
	samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/MiaPaCa2.//g" | sed "s/.cpm//g" | sed "s/.sorted//g"`
	
	if [[ $target =~ "H3K4me1" ]]
	then
		max=0.55
		min=0.08
		hmax=1.4

	
	elif [[ $target =~ "H3K4me3" ]];then
		max=1.7
		min=0.2
		hmax=4.2

	
	elif  [[ $target =~ "H3K27Ac" ]];then
		max=0.3
		min=0.02
                hmax=0.85


	elif  [[ $target =~ "H3K36me2" ]];then
                max=0.25
		min=0.02
                hmax=0.62


	elif  [[ $target =~ "H3K36me3" ]];then
                max=0.1
		min=0.03
                hmax=0.43


	fi


	plotHeatmap -m $f -o $target.png --colorMap RdBu_r --samplesLabel $samples --regionsLabel promoter enhancer --xAxisLabel "NicE-seq peak" --yMin $min --yMax $max --zMax $hmax

done
