#! /bin/bash

#start='date +%s'
cd ~/data/jkerriga/128DayOutput/fhd_61
vis_data=$(du -sh * |wc -l)
echo $vis_data
while [ $vis_data -ne 9 ];
do
    vis_data=$(du -sh *|wc -l)
done

#end='date+%s'
#runtime=$((end-start))
#echo $runtime