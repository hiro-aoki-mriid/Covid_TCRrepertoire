##Metaclonotype search on bulk TCR repertoire data
files="mydata/Pool/*.txt"
for filepath in ${files}
do
i=`basename ${filepath} .txt`
ipython mydata/tcrdist3_tabulating.metaclonotype_Minervina.py mydata/Pool/"${i}".txt
ipython mydata/tcrdist3_tabulating.metaclonotype_Snyder.py mydata/Pool/"${i}".txt
ipython mydata/tcrdist3_tabulating.metaclonotype_Francis.py mydata/Pool/"${i}".txt
done
