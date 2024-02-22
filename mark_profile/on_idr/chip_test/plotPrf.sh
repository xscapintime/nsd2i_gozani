for f in *.gz
#for f in `ls *H3*.gz | grep -v _D`
do

	bn=`basename $f | sed "s/.gz//g"`
	target=`echo $bn | cut -d . -f 1`
	#mark=`echo $bn | cut -d . -f 2`
	samples=`gunzip -c $f | head -n 1 | cut -d "[" -f 11 | cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/miapaca.//g" | sed "s/.cpm//g"`
	
	#if [[ $bn =~ "day_1" ]]
	#then
	#	max=2.5
	#	min=0
	
	#else
	#	max=6
        #       min=0

	#fi

	plotHeatmap -m $f -o $target.png --colorMap RdBu_r --samplesLabel $samples --xAxisLabel "NicE-seq IDR" 
done
