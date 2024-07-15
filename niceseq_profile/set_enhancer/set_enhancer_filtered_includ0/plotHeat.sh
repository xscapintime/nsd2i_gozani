for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	echo $bn
	
	regions=`gunzip -c $f | head -n 1 | cut -d "[" -f 9 |  cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/.bed//g" | sed "s/_prom_1k//g"`


	if [[ $bn =~ "E2F" ]];then
                ymax=0.093
                ymin=0.043

	elif [[ $bn =~ "EPITHELIAL" ]];then
                ymax=0.08
                ymin=0.04

        elif [[ $bn =~ "HEDGEHOG" ]];then
                ymax=0.075
                ymin=0.037

	elif [[ $bn =~ "INTERFERON" ]];then
                ymax=0.076
                ymin=0.04

        elif [[ $bn =~ "KRAS_SIGNALING_DN" ]];then
                ymax=0.076
                ymin=0.04

        elif [[ $bn =~ "KRAS_SIGNALING_UP" ]];then
                ymax=0.073
                ymin=0.037

	elif [[ $bn =~ "MYC_V1" ]];then
                ymax=0.11
                ymin=0.044

	elif [[ $bn =~ "MYC_V2" ]];then
                ymax=0.12
                ymin=0.044

        elif [[ $bn =~ "pan" ]];then
                ymax=0.075
                ymin=0.035
        
     	elif [[ $bn =~ "PRC" ]];then
                ymax=0.07
                ymin=0.038


	fi

	echo $bn $ymax $zmax

	plotHeatmap -m $f -o $bn.pdf --colorMap Blues --refPointLabel 'Enhancer' --regionsLabel $regions --yMin $ymin --yMax $ymax --heatmapHeight 20 #--colors "#7CCD7C" "#CDC9C9"  # --zMax $zmax --regionsLabel $regions #--colorList "#9F79EE,#8B8989"
	#plotProfile -m $f -o $bn.prof.pdf --refPointLabel 'Promoter' --yMin $ymin --yMax $ymax --regionsLabel "Gene set" "Controls" --colors "#7CCD7C" "#CDC9C9" --plotWidth 5.5 --plotTitle $regions



done
