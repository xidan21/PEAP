#!/bin/bash




R --vanilla --no-save --args $1 < extract_genes_title.r > temp.log

R --vanilla --no-save --args ../source/ranked_ensembl_ids.txt < template.r $2 > ../result/mes.log 
