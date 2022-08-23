######Covid mRNA vaccine, repertoire analysis
#####Figure3: Search SARS-Cov-2 Spike reactive clones using public database and tcrdist3

###Preprocess data using vdjtools
##Run the following docker container, and perform analysis in the "tmp" directory
#haokimriid/vdjtools
##Usage
#docker run -v /mnt/"Your working directory"/:/tmp/ -it haokimriid/vdjtools:0.0.1 /bin/bash
#cd ..
#cd tmp
#(Ex) docker run -v /mnt/c/Users/Owner/Desktop/COVID_paper/:/tmp/ -it haokimriid/vdjtools:0.0.2 /bin/bash
#Prepare dataset of bulk TCRseq for searchin metaclonotypes
sample=(001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 041)
Tcell=(CD8)
for k in "${sample[@]}"
do
for l in "${Tcell[@]}"
do
java -Xmx32G -jar vdjtools-1.2.1/vdjtools-1.2.1.jar PoolSamples -i aaVJ\
 Fig1_Fig2/raw.data/NaraCOVID_TP1_"${l}"_"${k}"_TCR.txt\
 Fig1_Fig2/raw.data/NaraCOVID_TP2_"${l}"_"${k}"_TCR.txt\
 Fig1_Fig2/raw.data/NaraCOVID_TP3_"${l}"_"${k}"_TCR.txt\
 Fig1_Fig2/raw.data/NaraCOVID_TP4_"${l}"_"${k}"_TCR.txt\
 Fig3/Pool/NaraCOVID_"${l}"_"${k}"_TCR
done
done
rm Fig3/Pool/*pool*.summary.txt

###Search SARS-Cov-2 reactive clones using tcrdist3
##Run the following docker container
#haokimriid/covid_tcrdist3.0
##Usage
#docker run -v /mnt/"Your working directory"/:/mydata/ -it haokimriid/covid_tcrdist3.0
#(Ex) docker run -v /mnt/c/Users/Owner/Desktop/COVID_paper/Fig3/:/mydata/ -it haokimriid/covid_tcrdist3.0
##Rscripts are run outside docker container (RStudio etc)

#Prepare data for tcrdist3 analysis from public dataset
Rscript mydata/Scripts_for_tcrdist3/Prepare_tcrdist3_input_Minervina.R #For Minervina dataset
Rscript mydata/Scripts_for_tcrdist3/Prepare_tcrdist3_input_Francis.R #For Francis dataset
#Generate metaclonotypes
#FIles are output in hla_restricted_meta_clonotypes 
python3 -c "from tcrsampler.setup_db import install_all_next_gen; install_all_next_gen(dry_run = False)" #For Snyder's dataset
ipython mydata/Scripts_for_tcrdist3/tcrdist3_metaclonotype_covid_Minervina.py #For Minervina dataset
ipython mydata/Scripts_for_tcrdist3/tcrdist3_metaclonotype_covid_Francis.py #For Francis dataset

#Prepare tcrdist3 output for metaclonotype search
Rscript mydata/Scripts_for_tcrdist3/Prepare_Metaclonotype_Search_input_Minervina.R
Rscript mydata/Scripts_for_tcrdist3/Prepare_Metaclonotype_Search_input_Francis.R

##Metaclonotype search on bulk TCR repertoire data
mkdir mydata/Scripts_for_tcrdist3/tcrdist3_output
sh mydata/Scripts_for_tcrdist3/metaclonotype_search.sh

###Exit docker container, then run "Fig3/Fig3_master.R"



