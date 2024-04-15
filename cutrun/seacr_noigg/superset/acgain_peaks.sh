for bar in 100000 110000 111000 111111
do

	cat H3K27Ac.binarycode.bed | awk -v var="$bar" '$5 == var' > H3K27Ac.$bar.bed
	echo  $bar
done
