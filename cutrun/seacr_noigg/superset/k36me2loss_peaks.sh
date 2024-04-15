for bar in 000111 100111 110111 111111
do

	cat H3K36me2.binarycode.bed | awk -v var="$bar" '$5 == var' > H3K36me2.$bar.bed
	echo  $bar
done

