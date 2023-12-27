###Script for "CD8+ T-cell memory induced by successive SARS-CoV-2 mRNA vaccinations is characterized by clonal replacement" by Aoki et al.

##This program run on docker container. Please pull the following docker images
#haokimriid/r-base
#haokimriid/beta_binomial
#haokimriid/vdjtools:0.0.1
#haokimriid/singlecell
#haokimriid/covid_tcrdist3.0

###Download data
##For bulk TCRseq analysis
wget --load-cookies /tmp/cookies.txt "https://drive.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.google.com/uc?export=download&id=1OWhxqi2GrcMxny0tvxfu5DLUxLLjFDkv' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1OWhxqi2GrcMxny0tvxfu5DLUxLLjFDkv" -O raw.data.tar.bz2 && rm -rf /tmp/cookies.txt
tar -jxvf raw.data.tar.bz2
mv raw.data data/original
##For beta-binomial tests
wget --load-cookies /tmp/cookies.txt "https://drive.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.google.com/uc?export=download&id=1z1ce87jLHPUFOtOFw45Qtpevu_QV-KH7' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1z1ce87jLHPUFOtOFw45Qtpevu_QV-KH7" -O differential_abundance.tar.bz2 && rm -rf /tmp/cookies.txt
tar -jxvf differential_abundance.tar.bz2
mv differential_abundance data/public
##For SCT analysis
wget --load-cookies /tmp/cookies.txt "https://drive.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.google.com/uc?export=download&id=1B8dI20Z4XY7RTkJjFO_GlEZPlXq5VeJ6' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1B8dI20Z4XY7RTkJjFO_GlEZPlXq5VeJ6" -O SCT.tar.bz2 && rm -rf /tmp/cookies.txt
tar -jxvf SCT.tar.bz2
mv SCT data/original

########################## Supplementary Tables: Summarize clonotype counts ####################################
##Sammarize clonotype counts for each dataset
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/SupTables/Repertoire_info_summary.R

########################## Figure 1: Temporal analysis on vaccine responding clonotypes ####################################
##Convert to Immunoseq format
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig1_beta-binomial/Convert_Immunoseq.R

###Run beta-binomial using merged training data
docker run -v ./:/home/ haokimriid/beta_binomial python2.7 home/data/public/differential_abundance/rundiffabBatch_2017_09_24.py --batchfile home/metadata/batchfile_covid.tsv --config home/data/public/differential_abundance/configuration001.ini --train home/data/public/differential_abundance/TrainingTSVs/replicates_Standard.csv --tsvDir home/result/intermediate/1_beta-binomial/Immunoseq --outDir home/result/intermediate/1_beta-binomial/binomial_test --parallel

#Tracking by JoinSamples
sample=(001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 041)
Tcell=(CD4 CD8)
for k in "${sample[@]}"
do
for l in "${Tcell[@]}"
do
docker run -v ./:/data/mydata/ haokimriid/vdjtools:0.0.1 java -Xmx32G -jar vdjtools-1.2.1/vdjtools-1.2.1.jar JoinSamples -x 2 -i strict\
 mydata/data/original/raw.data/NaraCOVID_TP1_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP2_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP3_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP4_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP5_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP6_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP8_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP9_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP10_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP11_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/result/intermediate/1_beta-binomial/JoinTP/NaraCOVID_"${l}"_"${k}"_TCR
done
done

#ID18 did not have TP11 (because they have 4th vaccination before TP11)
# -> JoinSamples on TP1-TP10, add dummy column for TP11
sample=(018)
Tcell=(CD4 CD8)
for k in "${sample[@]}"
do
for l in "${Tcell[@]}"
do
docker run -v ./:/data/mydata/ haokimriid/vdjtools:0.0.1 java -Xmx32G -jar vdjtools-1.2.1/vdjtools-1.2.1.jar JoinSamples -x 2 -i strict\
 mydata/data/original/raw.data/NaraCOVID_TP1_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP2_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP3_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP4_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP5_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP6_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP8_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP9_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP10_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/result/intermediate/1_beta-binomial/JoinTP/NaraCOVID_"${l}"_"${k}"_TCR
done
done
rm result/intermediate/1_beta-binomial/JoinTP/*join*.summary.txt
#Append dummy data
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig1_beta-binomial/ID18_dummy_add.R

###Append the result of Binomial test
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig1_beta-binomial/Append_betabinomial.R

###Summarize results
#Fig1C,D: Extent and diversity of responders
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig1_beta-binomial/Summarize_Responder.R
#Fig1E: Example of response patterns
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig1_beta-binomial/Combine_Responders.R
#Fig1F,G: Summarize response patterns
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig1_beta-binomial/Summarize_TP12Exp_TotalFreq.R

########################## Figure 2: scTCRseq analysis ####################################

###Run the following docker container, and perform analysis in the "tmp" directory
#haokimriid/singlecell
##Usage
#docker run --rm -p 8888:8888 -v /mnt/"Your working directory"/:/work haokimriid/singlecell jupyter notebook --allow-root
#docker run --rm -p 8888:8888 -v /mnt/d/COVID/:/work haokimriid/singlecell jupyter notebook --allow-root

##Following analysis are performed on jupyter notebook on docker images
#Preprocessing single-cell data: scTCR_preprocess.ipynb
#Create graph: scTCR_plotting.ipynb

########################## Figure 3: Analysis of response to the 3rd shot ####################################

#Summarize
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig3_P8/Summarize_TP12Exp_3rdResp.R

#Generate historgam for 2nd responder clones
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig3_P8/MainResp_Histogram.R

########################## Figure 4: Correlation analysis ####################################

#Calculate clonality at TP1
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig4_metaanalysis/Clonality_All.R

##Correlation analysis between repertoire parameters
#Fig4A and 4B: Create heatmap
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig4_metaanalysis/Metaanalysis_Correlation.R
#Merge repertoire parameters based on correlation heatmap
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig4_metaanalysis/Merge_repertoire_params.R
#Fig4C: Metaanalysis between CD4 and CD8 repertoire parameters
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig4_metaanalysis/Metaanalysis_Correlation_CD4vsCD8.R

#FigS4B-D: Metaanalysis between repertoire and clinical parameters
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig4_metaanalysis/Metaanalysis_vsClinical.R

########################## Figure5,6: TCR repertoire analysis on Tetramer+ T cells ####################################
#Combine the results of Tetramer TCRseq
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig5_Tetramer/Append_Tet.R

#Fig5: calculate total frequency of tetramer+ clonotypes
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig5_Tetramer/Tet_positivity_timepoint.R

#Fig6: Summarize response patterns of tetramer+ clonotypes
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig5_Tetramer/tetramer_positivity_TP12ExpFreq.R

########################## Figure7: TCR motif analysis for Spike specific clonotypes ####################################
### Preprocess
mkdir result/intermediate/7_tcrdist
#Prepare dataset of bulk TCRseq for searching metaclonotypes
sample=(001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 041)
Tcell=(CD8)
for k in "${sample[@]}"
do
for l in "${Tcell[@]}"
do
docker run -v ./:/data/mydata/ haokimriid/vdjtools:0.0.1 java -Xmx32G -jar vdjtools-1.2.1/vdjtools-1.2.1.jar PoolSamples -i aaVJ\
 mydata/data/original/raw.data/NaraCOVID_TP1_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP2_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP3_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP4_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP5_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP6_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP8_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP9_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP10_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP11_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/result/intermediate/7_tcrdist/Pool/NaraCOVID_"${l}"_"${k}"_TCR
done
done

#ID18
sample=(018)
Tcell=(CD8)
for k in "${sample[@]}"
do
for l in "${Tcell[@]}"
do
docker run -v ./:/data/mydata/ haokimriid/vdjtools:0.0.1 java -Xmx32G -jar vdjtools-1.2.1/vdjtools-1.2.1.jar PoolSamples -i aaVJ\
 mydata/data/original/raw.data/NaraCOVID_TP1_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP2_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP3_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP4_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP5_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP6_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP8_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP9_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/data/original/raw.data/NaraCOVID_TP10_"${l}"_"${k}"_TCR.pool.strict.table.txt\
 mydata/result/intermediate/7_tcrdist/Pool/NaraCOVID_"${l}"_"${k}"_TCR
done
done
rm result/intermediate/7_tcrdist/Pool/*pool*.summary.txt

###Search SARS-Cov-2 reactive clones using tcrdist3

#Prepare data for tcrdist3 analysis
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig7_tcrdist/Prepare_tcrdist3_input.R

#Generate metaclonotypes
#Files are output in hla_restricted_meta_clonotypes 
docker run -v ./:/mydata/ haokimriid/covid_tcrdist3.0 ipython mydata/code/Fig7_tcrdist/tcrdist3_metaclonotype_covid.py

#Prepare tcrdist3 output for metaclonotype search
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig7_tcrdist/Combine_tcrdist3_output.R
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig7_tcrdist/Prepare_Search_input_KMayerB.R

#Extract Metaclonotype query for participants, considering their HLA
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig7_tcrdist/Prepare_Search_input_participants.R

#Metaclonotype search on bulk TCR repertoire data
mkdir result/intermediate/7_tcrdist/tcrdist3_output
##Metaclonotype search on bulk TCR repertoire data
sample=(001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 035 036 037 038 039 040 041)
for k in "${sample[@]}"
do
docker run -v ./:/mydata/ haokimriid/covid_tcrdist3.0 ipython mydata/code/Fig7_tcrdist/tcrdist3_tabulating.metaclonotype.py\
 mydata/result/intermediate/7_tcrdist/Pool/NaraCOVID_CD8_"${k}"_TCR.pool.aaVJ.table.txt\
 mydata/result/intermediate/7_tcrdist/Metaclonotype_query/participants/"${k}".HLAI.KmayerB_metaclonotypes.tsv\
 tcrdist3_output/KmayerB
docker run -v ./:/mydata/ haokimriid/covid_tcrdist3.0 ipython mydata/code/Fig7_tcrdist/tcrdist3_tabulating.metaclonotype.py\
 mydata/result/intermediate/7_tcrdist/Pool/NaraCOVID_CD8_"${k}"_TCR.pool.aaVJ.table.txt\
 mydata/result/intermediate/7_tcrdist/Metaclonotype_query/participants/"${k}".HLAI.vdjdb_metaclonotypes.tsv\
 tcrdist3_output/vdjdb
docker run -v ./:/mydata/ haokimriid/covid_tcrdist3.0 ipython mydata/code/Fig7_tcrdist/tcrdist3_tabulating.metaclonotype.py\
 mydata/result/intermediate/7_tcrdist/Pool/NaraCOVID_CD8_"${k}"_TCR.pool.aaVJ.table.txt\
 mydata/result/intermediate/7_tcrdist/Metaclonotype_query/participants/"${k}".HLAI.Tet_metaclonotypes.tsv\
 tcrdist3_output/Tet
docker run -v ./:/mydata/ haokimriid/covid_tcrdist3.0 ipython mydata/code/Fig7_tcrdist/tcrdist3_tabulating.metaclonotype.py\
 mydata/result/intermediate/7_tcrdist/Pool/NaraCOVID_CD8_"${k}"_TCR.pool.aaVJ.table.txt\
 mydata/result/intermediate/7_tcrdist/Metaclonotype_query/participants/"${k}".HLAI.iedb_metaclonotypes.tsv\
 tcrdist3_output/iedb
done

#Merge metaclonotype information from tcrdist3
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig7_tcrdist/Append_tcrdist.R

#Export Epitope-specific clonotypes
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig7_tcrdist/tcrdist_positivity_summary.R

###Plotting figures
#FigS7A: Create Venn diagram for the correspondence between tcrdist and Tet-TCRseq
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig7_tcrdist/tcrdist_venn_diagram.R

#Fig7B: calculate total frequency of Epitope-specific clonotypes
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig7_tcrdist/tcrdist_positivity_timepoint.R

#Fig7C: Summarize response patterns of Epitope-specific clonotypes
docker run -v ./:/tmp/ haokimriid/r-base Rscript tmp/code/Fig7_tcrdist/tcrdist_positivity_TP12ExpFreq.R
