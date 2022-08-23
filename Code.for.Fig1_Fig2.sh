######Covid mRNA vaccine, repertoire analysis
#####Figure1 and Figure2: Bulk repertoire analysis on CD4 and CD8 T cells

###Run the following docker container, and perform analysis in the "tmp" directory
#haokimriid/vdjtools
##Usage
#docker run -v /mnt/"Your working directory"/:/tmp/ -it haokimriid/vdjtools:0.0.1 /bin/bash
#cd ..
#cd tmp
#(Ex) docker run -v /mnt/c/Users/Owner/Desktop/COVID_paper/:/tmp/ -it haokimriid/vdjtools:0.0.2 /bin/bash

###download raw.data
wget --load-cookies /tmp/cookies.txt "https://drive.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate\
 'https://drive.google.com/uc?export=download&id=1330po_xfL7Tj8amLnOUThTAqQ2Q557ZG' -O-\
 | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1330po_xfL7Tj8amLnOUThTAqQ2Q557ZG" -O Fig1_Fig2/raw.data.tar.gz && rm -rf /tmp/cookies.txt
###Unzip rawdata
tar -zxvf Fig1_Fig2/raw.data.tar.gz

###Fig1C-E: WGCNA analysis for variation patterns in clone frequency
#Preparation: Pool 4 timepoints in each individual
sample=(001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 041)
Tcell=(CD4 CD8)
for k in "${sample[@]}"
do
for l in "${Tcell[@]}"
do
java -Xmx32G -jar vdjtools-1.2.1/vdjtools-1.2.1.jar JoinSamples -x 1 -i strict\
 Fig1_Fig2/raw.data/NaraCOVID_TP1_"${l}"_"${k}"_TCR.txt\
 Fig1_Fig2/raw.data/NaraCOVID_TP2_"${l}"_"${k}"_TCR.txt\
 Fig1_Fig2/raw.data/NaraCOVID_TP3_"${l}"_"${k}"_TCR.txt\
 Fig1_Fig2/raw.data/NaraCOVID_TP4_"${l}"_"${k}"_TCR.txt\
 Fig1_Fig2/Join/NaraCOVID_"${l}"_"${k}"_TCR
done
done
rm Fig1_Fig2/Join/*join*.summary.txt

###Exit docker container, then run "Fig1_Fig2/Fig1_Fig2_master.R"
