#!/bin/bash

for f in *.gz; do
    bn=$(basename "$f" .gz)  # Remove the .gz extension

    regions=$(gunzip -c "$f" | head -n 1 | cut -d "[" -f 9 | cut -d "]" -f 1 | sed 's/"//g' | sed 's/,/ /g' | sed 's/.bed//g' | sed 's/_prom_1k//g')

    ymin=0
    ymax=0
    zmax=0

    if [[ $bn =~ "E2F" ]]; then
        ymax=1.05
        zmax=1.43
    elif [[ $bn =~ "INTERFERON" ]]; then
        ymax=0.65
        zmax=1.25
    elif [[ $bn =~ "KRAS_SIGNALING_DN" ]]; then
        ymax=0.5
        zmax=0.85
    elif [[ $bn =~ "KRAS_SIGNALING_UP" ]]; then
        ymax=0.6
        zmax=1.2
    elif [[ $bn =~ "MYC_V1" ]]; then
        ymax=1.05
        zmax=1.65
    elif [[ $bn =~ "MYC_V2" ]]; then
        ymax=1.05
        zmax=1.45
    elif [[ $bn =~ "PRC" ]]; then
        ymax=0.38
        zmax=0.75
    fi

    echo "$bn $ymax $zmax"

    plotHeatmap -m "$f" -o "$bn.png" --colorMap Blues --refPointLabel 'Promoter' --yMin "$ymin" --yMax "$ymax" --zMax "$zmax" --samplesLabel "$regions"

done
