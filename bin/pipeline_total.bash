#!/bin/bash



rm ../result/*.*
rm ../source/*.*

bash pipeline_pathway_gsea.bash $1 $2 &
bash pipeline_pathway_draw_map.bash $1 $2 & 

wait
bash pipeline_combine.bash $2

















