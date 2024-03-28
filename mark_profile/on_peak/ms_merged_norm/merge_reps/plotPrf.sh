for f in NSD2i*.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	target=`echo $bn | cut -d . -f 1 | cut -d _ -f 2`
	samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/MiaPaCa2.//g" | sed "s/.cpm//g" | sed "s/.sorted//g" |sed "s/.ms//g" |  tr -s ' ' '\n' | cut -d _ -f 1,2 | tr -s '\n' ' '`
	
	if [[ $target =~ "H3K4me1" ]]
	then
		max=0.3
		min=0.12
		hmax=1.3

	
	elif [[ $target =~ "H3K4me3" ]];then
		max=1.2
		min=0.05
		hmax=3.7

	
	elif  [[ $target =~ "H3K27Ac" ]];then
		max=0.17
		min=0.05
                hmax=0.8


	elif  [[ $target =~ "H3K36me2" ]];then
                max=0.15
		min=0.04
                hmax=0.6


	elif  [[ $target =~ "H3K36me3" ]];then
                max=0.092
		min=0.035
                hmax=0.43


	fi


	#plotHeatmap -m $f -o $target.png --colorMap RdBu_r --samplesLabel $samples --regionsLabel "$target CUT&RUN" --xAxisLabel "NicE-seq peak" --missingDataColor "#FFF6EB"  #--yMin $min --yMax $max --zMax $hmax

	plotProfile -m $f -o $target.profile.png --numPlotsPerRow 3 --samplesLabel $samples --refPointLabel "NicE-seq peak" --plotWidth 5.5 \
--averageType mean --plotType se --regionsLabel $target --colors  "#5773CC" "#5773CC" "#5773CC" "#8B8989" "#8B8989" "#8B8989" \
--legendLocation upper-left
	

done
