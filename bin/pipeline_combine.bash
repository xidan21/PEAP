#!/bin/bash




less ../result/mes.log | grep "\[1\] \"there" | sed 's/^.*path://' | sed 's/\//	/g' | sed 's/"//g' > ../result/gsea.txt

R --vanilla --no-save --args ../result/gsea.txt $1 < combine.r > temp.log		

