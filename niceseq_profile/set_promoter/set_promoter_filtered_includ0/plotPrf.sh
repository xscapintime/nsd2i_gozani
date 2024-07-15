for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	echo $bn
	
	regions=`gunzip -c $f | head -n 1 | cut -d "[" -f 9 |  cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/.bed//g" | sed "s/_prom_1k//g" | cut -d " " -f 1`

	ymin=0

	if [[ $bn =~ "E2F" ]];then
                ymax=1.05
                zmax=1.43

	elif [[ $bn =~ "EPITHELIAL" ]];then
                ymax=0.77
                zmax=1.25

	elif [[ $bn =~ "INTERFERON" ]];then
                ymax=0.77
                zmax=1.15

        elif [[ $bn =~ "KRAS_SIGNALING_DN" ]];then
                ymax=0.66
                zmax=0.85

        elif [[ $bn =~ "KRAS_SIGNALING_UP" ]];then
                ymax=0.72
                zmax=1.2

	elif [[ $bn =~ "MYC_V1" ]];then
                ymax=1.05
                zmax=1.65

	elif [[ $bn =~ "MYC_V2" ]];then
                ymax=1.1
                zmax=1.45


     	elif [[ $bn =~ "pan" ]];then
                ymax=0.7
                zmax=1.1


     	elif [[ $bn =~ "PRC" ]];then
                ymax=0.57
                zmax=0.75

	fi

	echo $bn $ymax $zmax

	#plotHeatmap -m $f -o $bn.png --colorMap Blues --refPointLabel 'Promoter' --yMin $ymin --yMax $ymax --zMax $zmax --regionsLabel $regions #--colorList "#9F79EE,#8B8989"
	plotProfile -m $f -o $bn.prof.pdf --refPointLabel 'TSS' --yMin $ymin --yMax $ymax --regionsLabel "Gene set" "Controls" --colors "#7CCD7C" "#CDC9C9" --plotWidth 5.5 --plotTitle $regions



done
