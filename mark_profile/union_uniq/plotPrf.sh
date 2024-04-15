for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	target=`echo $bn | cut -d . -f 1 | cut -d _ -f 2`
	samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' |  sed "s/,/ /g" |  tr -s ' ' '\n' | cut -d _ -f 1,2 | tr -s '\n' ' '`
	bedlist=`gunzip -c $f | head -n 1 | cut -d "[" -f 9 | cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" |  tr -s ' ' '\n'  | cut -d . -f 2 | tr -s '\n' ' '`

	plotHeatmap -m $f -o $target.png --colorMap RdBu_r --samplesLabel $samples --regionsLabel $bedlist  --xAxisLabel "peak" --missingDataColor "#FFF6EB"  #--yMin $min --yMax $max --zMax $hmax

	#plotProfile -m $f -o $target.profile.png --numPlotsPerRow 3 --samplesLabel $samples --refPointLabel "NicE-seq peak" --plotWidth 5.5 \
#--averageType mean --plotType se --regionsLabel $target --colors  "#5773CC" "#5773CC" "#5773CC" "#8B8989" "#8B8989" "#8B8989" \
#--legendLocation upper-left
	

done
