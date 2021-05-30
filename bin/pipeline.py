#!/usr/bin/env python

import os
import sys
import re
import os
import subprocess

ranked_gene_list = sys.argv[1]
gene_type_gene_list = sys.argv[2]
pat = sys.argv[3]
start = sys.argv[4]
end = sys.argv[5]
#r_path = sys.argv[6]

######################################################################

def convert_chr (key):
	dict = {
	
		'1':'I',
		'2':'II',
		'3':'III',
		'4':'IV',
		'5':'V',
		'6':'VI',
		'7':'VII',
		'8':'IIX',
		'9':'IX',
		'10':'X',
		'11':'XI',
		'12':'XII',
		'13':'XIII',
		'14':'XIV',
		'15':'XV',
		'16':'XVI',
		'17':'XVII',
		'18':'XVIII',
		'19':'XIX',
		'20':'XX',
		'21':'XXI',
		'22':'XXII',
		'23':'XXIII',
		'24':'XXIV',
		'25':'XXV',
		'26':'XXVI',
		'27':'XXVII',
		'28':'XXVIII',
		'29':'XXIX',
		'30':'XXX',
		'31':'XXXI',
		'32':'XXXII',
		'33':'XXXIII',
		'34':'XXXIV',
		'35':'XXXV',
		'36':'XXXVI',
		'37':'XXXVII',
		'38':'XXXVIII',
		'39':'XXXIX'
			
		}

	return dict[key]

######################################################################

def translate (key):	
	dict = {
		
"thorhynchus_anatinus":"oanatinus_gene_ensembl",
"Taeniopygia_guttata":"tguttata_gene_ensembl",
"Cavia_porcellus":"cporcellus_gene_ensembl",
"Gasterosteus_aculeatus":"gaculeatus_gene_ensembl",
"Loxodonta_africana":"lafricana_gene_ensembl",
"Ictidomys_tridecemlineatus":"itridecemlineatus_gene_ensembl",
"Myotis_lucifugus":"mlucifugus_gene_ensembl",
"Homo_sapiens":"hsapiens_gene_ensembl",
"Choloepus_hoffmanni":"choffmanni_gene_ensembl",
"Ciona_savignyi":"csavignyi_gene_ensembl",
"Felis_catus":"fcatus_gene_ensembl",
"Rattus_norvegicus":"rnorvegicus_gene_ensembl",
"Gallus_gallus":"ggallus_gene_ensembl",
"Tupaia_belangeri":"tbelangeri_gene_ensembl",
"Pelodiscus_sinensis":"psinensis_gene_ensembl",
"Xenopus_tropicalis":"mfuro_gene_ensembl",
"Equus_caballus":"xtropicalis_gene_ensembl",
"Callithrix_jacchus":"ecaballus_gene_ensembl",
"Pongo_abelii":"cjacchus_gene_ensembl",
"Danio_rerio":"pabelii_gene_ensembl",
"Xiphophorus_maculatus":"drerio_gene_ensembl",
"Tetraodon_nigroviridis":"xmaculatus_gene_ensembl",
"Tursiops_truncatus":"tnigroviridis_gene_ensembl",
"Latimeria_chalumnae":"ttruncatus_gene_ensembl",
"Saccharomyces_cerevisiae":"lchalumnae_gene_ensembl",
"Ailuropoda_melanoleuca":"scerevisiae_gene_ensembl",
"Caenorhabditis_elegans":"amelanoleuca_gene_ensembl",
"Macaca_mulatta":"celegans_gene_ensembl",
"Pteropus_vampyrus":"mmulatta_gene_ensembl",
"Monodelphis_domestica":"pvampyrus_gene_ensembl",
"Vicugna_pacos":"mdomestica_gene_ensembl",
"Anolis_carolinensis":"vpacos_gene_ensembl",
"Oreochromis_niloticus":"acarolinensis_gene_ensembl",
"Tarsius_syrichta":"oniloticus_gene_ensembl",
"Otolemur_garnettii":"tsyrichta_gene_ensembl",
"Takifugu_rubripes":"trubripes_gene_ensembl",
"Drosophila_melanogaster":"dmelanogaster_gene_ensembl",
"Petromyzon_marinus":"pmarinus_gene_ensembl",
"Erinaceus_europaeus":"eeuropaeus_gene_ensembl",
"Microcebus_murinus":"mmurinus_gene_ensembl",
"Oryzias_latipes":"olatipes_gene_ensembl",
"Pan_troglodytes":"ptroglodytes_gene_ensembl",
"Echinops_telfairi":"etelfairi_gene_ensembl",
"Ciona_intestinalis":"cintestinalis_gene_ensembl",
"Ochotona_princeps":"oprinceps_gene_ensembl",
"Gorilla_gorilla":"ggorilla_gene_ensembl",
"Dipodomys_ordii":"dordii_gene_ensembl",
"Nomascus_leucogenys":"nleucogenys_gene_ensembl",
"Sus_scrofa":"sscrofa_gene_ensembl",
"Mus_musculus":"mmusculus_gene_ensembl",
"Oryctolagus_cuniculus":"ocuniculus_gene_ensembl",
"Meleagris_gallopavo":"mgallopavo_gene_ensembl",
"Gadus_morhua":"gmorhua_gene_ensembl",
"Sorex_araneus":"saraneus_gene_ensembl",
"Dasypus_novemcinctus":"dnovemcinctus_gene_ensembl",
"Procavia_capensis":"pcapensis_gene_ensembl",
"Bos_taurus":"btaurus_gene_ensembl",
"Macropus_eugenii":"meugenii_gene_ensembl",
"Sarcophilus_harrisii":"sharrisii_gene_ensembl",
"Canis_familia":"cfamiliaris_gene_ensembl"
	}

	return dict[key]

########################################################################################

def R_path_name(path):
	

	os.system("whoami > ../source/getwd.log")

	x = open("../source/getwd.log",'r')

	for line in x:
        	if not re.search("^$",line):
	                print line
	                sub = re.search("^.*\\\(.*)\n",line)
	                user = sub.group(1)

	path_version = path[2:]
	path_version = path_version[:-2]

	#print path_version

	path_lib = "C:/Users/"+user+"/Documents/R/win-library/"+path_version+"/"

	return path_lib

########################################################################################

if __name__ == "__main__":

#	cmd = "\"c:/Program Files/R/"+r_path+"/bin/Rscript.exe\""

	#print cmd

#	r_lib = R_path_name(r_path)

	os.system("rm ../result/*")

	roman_names = ["Saccharomyces_cerevisiae","blabla"]

	for spe_name in roman_names:
		if re.search(species,spe_name):
			os.system("R --vanilla --no-save --args %s %s %s %s < extra_gene.r >../source/temp.log" %(translate(species),convert_chr(chr_num),start,end))	
			
			break
		else:
			continue	
	else:
		os.system("R --vanilla --no-save --args %s %s %s %s < extra_gene.r >../source/temp.log" %(translate(species),chr_num,start,end))

	os.system("python detect_snps_in_utr_region.py %s %s %s %s" %(chr_num,start,end,snps_file))
	os.system("python convert_to_cds.py")

	os.system("python detect_snps_in_flanking_region.py %s %s %s %s" %(chr_num,start,end,snps_file))
	os.system("python generating_two_seqs.py")
	os.system("python detect_non_synomynous.py %s %s %s %s" %(chr_num,start,end,snps_file))

	os.system("cat ../source/snps_in_flanking_regions.txt ../source/snps_in_utr_region.txt ../source/snps_in_coding_region.txt | sort -n -k2 > ../result/total.txt")

	os.system("python sort.py > ../result/scored_aa.txt")
	os.system("sort -n -k17 -r ../result/scored_aa.txt > ../result/ranked_aa.txt")
