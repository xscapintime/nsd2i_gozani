for f in PRC*.gz
do

	bn=`basename $f | sed "s/.gz//g"`
	echo $bn
	
	regions=`gunzip -c $f | head -n 1 | cut -d "[" -f 9 |  cut -d "]" -f 1 | sed 's/"//g' | sed "s/,/ /g" | sed "s/.bed//g" | sed "s/_enhancer//g" | cut -d " " -f 1`

	# ymin=0

	if [[ $bn =~ "E2F" ]];then
                ymax=0.093
                ymin=0.041

	elif [[ $bn =~ "EPITHELIAL" ]];then
                ymax=0.078
                ymin=0.041

        elif [[ $bn =~ "HEDGEHOG" ]];then
                ymax=0.075
                ymin=0.037

	elif [[ $bn =~ "INTERFERON" ]];then
                ymax=0.077
                ymin=0.041

        elif [[ $bn =~ "KRAS_SIGNALING_DN" ]];then
                ymax=0.08
                ymin=0.039

        elif [[ $bn =~ "KRAS_SIGNALING_UP" ]];then
                ymax=0.073
                ymin=0.039

	elif [[ $bn =~ "MYC_V1" ]];then
                ymax=0.105
                ymin=0.044

	elif [[ $bn =~ "MYC_V2" ]];then
                ymax=0.115
                ymin=0.044
        
        # elif [[ $bn =~ "pancreatic_neoplasia_suppressed" ]];then
        #         ymax=0.073
        #         ymin=0.043

        elif [[ $bn =~ "pan" ]];then
                ymax=0.077
                ymin=0.042


     	elif [[ $bn =~ "PRC" ]];then
                ymax=0.074
                ymin=0.038



	fi

	echo $bn $ymax $ymin


	plotProfile -m $f -o $bn.prof.pdf --refPointLabel 'Enhancer' --regionsLabel "Gene set" "Controls" --colors "#87CEFA" "#CDC9C9" --plotWidth 5.5 --plotTitle $regions --yMin $ymin --yMax $ymax



done
