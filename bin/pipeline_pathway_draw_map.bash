#!/usr/bin/bash




R --vanilla --no-save --args $1 $2 < extra_log2_value.r > ../source/temp.log

python kegg_link.py $2 > ../result/web_link.txt 



