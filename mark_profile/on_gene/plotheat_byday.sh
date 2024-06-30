for f in *.5kdist.1mbreg.10bin.gz
do

	
	echo subset by day for $bn



computeMatrixOperations subset -m H3K27me3.genebody.5kdist.1mbreg.10bin.gz -o H3K27me3.d1.genebody.5kdist.1mbreg.10bin.gz --samples "NSD2i_D1_H3K27me3" "Vehicle_D1_H3K27me3"
