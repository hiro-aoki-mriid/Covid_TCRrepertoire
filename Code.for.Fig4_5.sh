######Covid mRNA vaccine, repertoire analysis
#####Figure 4, 5: Repertoire analysis on Spike reactive CD4 and CD8 T cells determined by AIM assay

###Run the following docker container, and perform analysis in the "tmp" directory
#haokimriid/vdjtools
##Usage
#docker run -v /mnt/"Your working directory"/:/tmp/ -it haokimriid/vdjtools:0.0.1 /bin/bash
#cd ..
#cd tmp
#(Ex) docker run -v /mnt/c/Users/Owner/Desktop/COVID_paper/:/tmp/ -it haokimriid/vdjtools:0.0.2 /bin/bash

###download raw.data
wget --load-cookies /tmp/cookies.txt "https://drive.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate\
 'https://drive.google.com/uc?export=download&id=1nJ-L32WtsgZ2dfIyPfzt8ytIWDFfKK82' -O-\
 | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1nJ-L32WtsgZ2dfIyPfzt8ytIWDFfKK82" -O Fig4_Fig5/raw.data.AIM.tar.gz && rm -rf /tmp/cookies.txt
###Unzip rawdata
tar -zxvf Fig4_Fig5/raw.data.AIM.tar.gz

###Fig4_Fig5B: Determine AIM+ clones 
##Create Scatter plot of AIM-NonNV overlapping clones and define AIM+ clones by their enrichment into AIM relative to NonNV
#Generate data table by VDJtools
sample=(001 002 003 004 005 008 009 012 021 023 027 031 034 035 038 040)
Tcell=(CD4 CD8)
mkdir Fig4_Fig5/Join
for k in "${sample[@]}"
do
for l in "${Tcell[@]}"
do
java -Xmx32G -jar vdjtools-1.2.1/vdjtools-1.2.1.jar JoinSamples -x 1 -i strict\
 Fig4_Fig5/raw.data.AIM/COVID_AIM_AIM_"${l}"_"${k}"_TCR.txt\
 Fig4_Fig5/raw.data.AIM/COVID_AIM_NonNV_"${l}"_"${k}"_TCR.txt\
 Fig4_Fig5/Join/COVIDAIM_"${l}"_"${k}"_TCR
java -Xmx32G -jar vdjtools-1.2.1/vdjtools-1.2.1.jar OverlapPair -i strict\
 Fig4_Fig5/raw.data.AIM/COVID_AIM_AIM_"${l}"_"${k}"_TCR.txt\
 Fig4_Fig5/raw.data.AIM/COVID_AIM_NonNV_"${l}"_"${k}"_TCR.txt\
 Fig4_Fig5/Join/COVIDAIM_"${l}"_"${k}"_TCR
done
done
rm Fig4_Fig5/Join/*join*.summary.txt
cat Fig4_Fig5/Join/*paired.strict.summary.txt > Fig4_Fig5/AIM_NonNV_OL_summary.txt #Used for S3A
rm Fig4_Fig5/Join/*paired.strict*.txt

###Exit docker container, then run "Fig4_Fig5/Fig4_Fig5_master.R"

