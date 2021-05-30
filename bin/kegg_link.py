#!/usr/bin/env python

import re
import os
import sys
import operator
from os import listdir

# 22273 ENSMUSG00000107235 NA
# 67003 ENSMUSG00000030884 -0.00470605258785215





def extract_file_titles(species, head_file):

	files = listdir(("../lib/" + species + "/"))
	pathway_files = []

	for i in xrange(0, len(files)):

		 if re.search(head_file, files[i]):	
	
			pathway_files.append(files[i][:-4])
	return pathway_files
	
def url_link(extra_log2_value_file):
	path_id = str(extra_log2_value_file)[10:18] # mmu00190
	print path_id + "\t",
	# http://www.kegg.jp/kegg-bin/show_pathway?hsa04530/default%3white/hsa:6709%09%23ff1493/hsa:8531%09%230000ff,red

	x = open(extra_log2_value_file,'r') # extra_log2_value.r

	kegg_line = "http://www.kegg.jp/kegg-bin/show_pathway?"+path_id+"/default%3white"

	for line in x:

		if re.search("^\d+\s+\w+\d+\s+.*$", line):

			#print line

			sub = re.search("^(\d+)\s+\w+\d+\s+(.*)$", line)











			entrez_id = sub.group(1)

			log_value = sub.group(2)

			if operator.ne(log_value, "NA"):



				if operator.gt(float(log_value),0):

					# 255,255,255 -> white
					# 255,0,0 -> red
					# 255,20,147 -> pink
					
					color_code_1_red = hex(255).split('x')[1]
					color_code_2_red = hex(255 - int(255*float(log_value))).split('x')[1]

					color_code_total_red = str(color_code_1_red)+str(color_code_2_red)+str(color_code_2_red)
	#				print "/mmu:"+str(entrez_id)+"%09%23ff1493"
				
					kegg_line = kegg_line+"/hsa:"+str(entrez_id)+"%09%23"+str(color_code_total_red)
				#	kegg_line = kegg_line+"/mmu:"+str(entrez_id)+"%09%23ff1493"



				elif operator.lt(float(log_value),0):

					
					# 255,255,255 -> white
					# 0,0,255 -> blue

					color_code_2_blue = hex(255).split('x')[1]
					color_code_1_blue = hex(255 - abs(int(255*float(log_value)))).split('x')[1]

					color_code_total_blue = str(color_code_1_blue)+str(color_code_1_blue)+str(color_code_2_blue)

					kegg_line = kegg_line+"/hsa:"+str(entrez_id)+"%09%23"+str(color_code_total_blue)
				#	kegg_line = kegg_line+"/mmu:"+str(entrez_id)+"%09%230000ff,red"

			else:
				kegg_line = kegg_line+"/hsa:"+str(entrez_id)
			    
		else:

			kegg_line = " "


	print str(kegg_line)

	#webbrowser.open_new(str(kegg_line))

	#os.system("python -m webbrowser -t '%s'" %(str(kegg_line.rstrip())))











if __name__=='__main__':


	if operator.eq(sys.argv[1],"Human"):

		pathway_files = extract_file_titles("human", "hsa")


	elif operator.eq(sys.argv[1],"Mouse"):

		pathway_files = extract_file_titles("mouse", "mmu")

	elif operator.eq(sys.argv[1],"Rat"):

		pathway_files = extract_file_titles("Rat", "mmu")


	for i in xrange(0, len(pathway_files)):

		url_link("../source/" + pathway_files[i] + "_entrez_id_with_log2_value.txt")
	


