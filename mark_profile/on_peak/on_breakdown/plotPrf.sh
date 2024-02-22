for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	target=`echo $bn | cut -d . -f 1,2`
	samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/MiaPaCa2.//g" | sed "s/.cpm//g" | sed "s/.sorted//g"`
	
	if [[ $target =~ "H3K4me1" ]]
	then
		max=1
		hmax=2.5

	
	elif [[ $target =~ "H3K4me3" ]];then
		max=3.5
		hmax=8

	
	elif  [[ $target =~ "H3K27Ac" ]];then
		max=0.6
                hmax=1.8


	elif  [[ $target =~ "H3K36me2" ]];then
                max=0.4
                hmax=1


	elif  [[ $target =~ "H3K36me3" ]];then
                max=0.3
                hmax=0.75


	fi


	plotHeatmap -m $f -o $target.png --colorMap RdBu_r  --regionsLabel promoter enhancer --samplesLabel $samples --yMin 0 --yMax $max --zMax $hmax --xAxisLabel "NicE-seq peak"
done
