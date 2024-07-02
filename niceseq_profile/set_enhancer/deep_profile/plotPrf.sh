for f in *.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	echo $bn
	
	regions=`gunzip -c $f | head -n 1 | cut -d "[" -f 9 |  cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/.bed//g" | sed "s/_enhancer//g" | cut -d " " -f 1`

	# ymin=0


	if [[ $bn =~ "E2F" ]];then
                ymax=0.091
                ymin=0.04

	elif [[ $bn =~ "EPITHELIAL" ]];then
                ymax=0.075
                ymin=0.035

        elif [[ $bn =~ "HEDGEHOG" ]];then
                ymax=0.067
                ymin=0.031

	elif [[ $bn =~ "INTERFERON" ]];then
                ymax=0.071
                ymin=0.04

        elif [[ $bn =~ "KRAS_SIGNALING_DN" ]];then
                ymax=0.071
                ymin=0.04

        elif [[ $bn =~ "KRAS_SIGNALING_UP" ]];then
                ymax=0.07
                ymin=0.032

	elif [[ $bn =~ "MYC_V1" ]];then
                ymax=0.11
                ymin=0.046

	elif [[ $bn =~ "MYC_V2" ]];then
                ymax=0.12
                ymin=0.041

     	elif [[ $bn =~ "PRC" ]];then
                ymax=0.063
                ymin=0.036


	fi

	echo $bn $ymax $ymin

	#plotHeatmap -m $f -o $bn.png --colorMap Blues --refPointLabel 'Promoter' --yMin $ymin --yMax $ymax --ymin $ymin --regionsLabel $regions #--colorList "#9F79EE,#8B8989"
	plotProfile -m $f -o $bn.prof.pdf --refPointLabel 'Enhancer' --regionsLabel "Gene set" "Control" --colors "#9F79EE" "#8B8989" --plotWidth 5.5 --plotTitle $regions --yMin $ymin --yMax $ymax



done
