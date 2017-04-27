import pysam
import argparse
import csv
import sys


from scipy.stats import spearmanr, kendalltau, pearsonr


ap = argparse.ArgumentParser()
ap.add_argument('exacFreq', help='vcf from exac')
ap.add_argument('convergeFreq', help='Mapped reads in bam format')
ap.add_argument('out', help='for example, /u/home/s/serghei/project/CONVERGE/clinvar/legend/chr1.legend')


args = ap.parse_args()





#1,13380,.,C,G,0,1692,0.0
#------

dictExac={}
exacSet=set()
print "Read ExAC freq ...", args.exacFreq


#1,69395,.,G,A,0,39426,1.0

with open(args.exacFreq) as csvfile:
    readCSV = csv.reader(csvfile)
    for row in readCSV:
        if row!=[]:
            pos=row[1]
            rsid=row[2]
            REF=row[3]
            ALT=row[4]
            REF_COUNT=row[5]
            ALT_COUNT=row[6]
            
            MAF=row[7]
            dictExac[pos]=(rsid,REF,ALT,MAF,REF_COUNT,ALT_COUNT)
            exacSet.add(pos)


print "Number of variants from ExAC", len(exacSet)









#1       10388   2       21280   A:21126 T:154
print "Read CONVERGE freq ...", args.convergeFreq


freqList_exac=[]
freqList_converge=[]

freqList=[]

inconsistent_exac_converge=0

with open(args.convergeFreq) as csvfile:
    readCSV = csv.reader(csvfile, delimiter="\t")
    for row in readCSV:
        if row!=[] and "POS" !=row[1]:
            pos=row[1]
            allele_1=row[4].split(":")[0]
            allele_2=row[5].split(":")[0]
            
            allele_1_count=int(row[4].split(":")[1])
            allele_2_count=int(row[5].split(":")[1])
            
            freqCONVERGE=allele_2_count/21280.0
            

            if pos in exacSet:
                
                freqEXAC=1-float(dictExac[pos][3])
                
                REF_EXAC_COUNT=dictExac[pos][5]
                ALT_EXAC_COUNT=dictExac[pos][4]
                
                #('rs75062661', 'A', 'G', '0.0268551907383')
                if allele_1==dictExac[pos][1] and allele_2==dictExac[pos][2]:
                    #print pos,allele_1,allele_2,freqEXAC,freqCONVERGE, dictExac[pos], row
                    freqList_exac.append(freqEXAC)
                    freqList_converge.append(freqCONVERGE)
                    freqList.append((freqEXAC,freqCONVERGE,pos,allele_1,allele_2,REF_EXAC_COUNT,ALT_EXAC_COUNT,allele_1_count,allele_2_count))
                
                else:
                    inconsistent_exac_converge+=1
                
                #print row
                #print allele_1,allele_1_count
                #print allele_2,allele_2_count
                #print dictExac[pos]
                #sys.exit(1)


print "inconsistent_exac_converge",inconsistent_exac_converge


print pearsonr(freqList_converge,freqList_exac)



out=open(args.out,"w")


#header
out.write("freqEXAC,freqCONVERGE,pos,allele_1,allele_2,REF_EXAC_COUNT,ALT_EXAC_COUNT,allele_1_count,allele_2_count")
out.write("\n")

for i in freqList:
    out.write(str(i[0])+","+str(i[1])+","+str(i[2])+","+str(i[3])+","+str(i[4])+","+str(i[5])+","+str(i[6])+","+str(i[7])+","+str(i[8]))
    out.write("\n")



out.close()

print "Done!"







