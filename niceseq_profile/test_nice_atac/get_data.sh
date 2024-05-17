rsync -av --include="*/" --include="*Results.html" --exclude="*" --exclude='*/' 'lyshi@narval.alliancecan.ca:/home/lyshi/scratch/nsd2i_gozani_2402/niceseq_profile/test_nice_atac/' .
rsync -av --include="*/" --include="knownResults.txt" --exclude="*" --exclude='*/' 'lyshi@narval.alliancecan.ca:/home/lyshi/scratch/nsd2i_gozani_2402/niceseq_profile/test_nice_atac/' .
rsync -av 'lyshi@narval.alliancecan.ca:/home/lyshi/scratch/nsd2i_gozani_2402/niceseq_profile/test_nice_atac/data/wget.sh' data/
