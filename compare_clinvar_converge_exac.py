import pysam
import argparse
import csv
import sys


from scipy.stats import spearmanr, kendalltau, pearsonr


ap = argparse.ArgumentParser()
ap.add_argument('exacFreq', help='vcf from exac')
ap.add_argument('convergeFreq', help='Mapped reads in bam format')
ap.add_argument('chr', help='Mapped reads in bam format')


ap.add_argument('out', help='for example, /u/home/s/serghei/project/CONVERGE/clinvar/legend/chr1.legend')


args = ap.parse_args()



clinvarVCF="/u/home/s/serghei/project/CONVERGE/clinvar/new/clinvar/output/b37/single/clinvar_alleles.single.b37.vcf"

#clinvar
#1	976059	.	C	T	.	.	MEASURESET_TYPE=Variant;MEASURESET_ID=210111;RCV=RCV000195231;ALLELE_ID=206691;SYMBOL=AGRN;HGVS_C=NM_198576.3:c.526C>T;HGVS_P=NP_940978.2:p.Leu176_eq_;MOLECULAR_CONSEQUENCE=NM_198576.3:c.526C>T:synonymous_variant;CLINICAL_SIGNIFICANCE=Likely_benign|Uncertain_significance;PATHOGENIC=0;BENIGN=1;CONFLICTED=0;REVIEW_STATUS=criteria_provided..conflicting_interpretations;GOLD_STARS=1;ALL_SUBMITTERS=Genetic_Services_Laboratory..University_of_Chicago|PreventionGenetics..PreventionGenetics;ALL_TRAITS=not_specified|not_specified|NOT_SPECIFIED;ALL_PMIDS=25741868;ORIGIN=germline;XREFS=MedGen:CN169374


dictLinvar={}
posClinvar=set()

with open(clinvarVCF) as csvfile:
    readCSV = csv.reader(csvfile,delimiter="\t")
    next(readCSV)
    for row in readCSV:
        
        if row[0]==args.chr:
            pos=row[1]
            ref=row[3]
            alt=row[4]
            BENIGN=row[7].split("BENIGN=")[1].split(";")[0]
            PATHOGENIC=row[7].split("PATHOGENIC=")[1].split(";")[0]
            CONFLICTED=row[7].split("CONFLICTED=")[1].split(";")[0]
            
            
            #if CONFLICTED!="1":
            dictLinvar[pos+ref+alt]=(ref,alt,row[7],PATHOGENIC,BENIGN,CONFLICTED)
            posClinvar.add(pos+ref+alt)






print "Number of CLINVAR variants=",len(posClinvar)








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
            
            MAF=float(row[7])
            if (1-MAF)<0.05:
                dictExac[pos+ref+alt]=(rsid,REF,ALT,1-MAF,REF_COUNT,ALT_COUNT)
                exacSet.add(pos+ref+alt)


print args.exacFreq, "Number of minor variants from ExAC ", len(exacSet)







#CONVERGE ------------------------

#1       10388   2       21280   A:21126 T:154
print "Read CONVERGE freq ...", args.convergeFreq




freqList=[]

freqList_BENIGN=[]

freqList_notEXAC_CONVERGE=[]

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
            
            #freqList_notEXAC_CONVERGE
            
            
            
          
            
            if (pos+allele_1+allele_2) not in exacSet and (pos+allele_1+allele_2) in posClinvar:
                if allele_1==dictLinvar[pos+allele_1+allele_2][0] and allele_2==dictLinvar[pos+allele_1+allele_2][1]:
                    print pos+allele_1+allele_2,freqCONVERGE
                    
                    if freqCONVERGE>0.05:
                            if dictLinvar[pos+allele_1+allele_2][3]=="1":
                                freqList_notEXAC_CONVERGE.append((freqCONVERGE,pos,allele_1,allele_2,allele_1_count,allele_2_count,dictLinvar[(pos,allele_1,allele_2)][0],dictLinvar[(pos,allele_1,allele_2)][1],dictLinvar[(pos,allele_1,allele_2)][2]))
                                print (freqCONVERGE,pos,allele_1,allele_2,allele_1_count,allele_2_count,dictLinvar[(pos,allele_1,allele_2)][0],dictLinvar[(pos,allele_1,allele_2)][1],dictLinvar[(pos,allele_1,allele_2)][2])
        
        
        
                
                
                
            elif (pos+allele_1+allele_2) in exacSet and (pos+allele_1+allele_2) in posClinvar:
                freqEXAC=1-float(dictExac[pos+allele_1+allele_2][3])
                
                REF_EXAC_COUNT=dictExac[pos+allele_1+allele_2][5]
                ALT_EXAC_COUNT=dictExac[pos+allele_1+allele_2][4]
                
                if 1-freqEXAC <0.05 and freqCONVERGE>0.04:
                                    
                
                
                        
                        
                        
                    if dictLinvar[pos+allele_1+allele_2][3]=="1" and dictLinvar[pos+allele_1+allele_2][5]=="0":
                            freqList.append((1-freqEXAC,freqCONVERGE,pos,allele_1,allele_2,REF_EXAC_COUNT,ALT_EXAC_COUNT,allele_1_count,allele_2_count,dictLinvar[pos+allele_1+allele_2][0],dictLinvar[pos+allele_1+allele_2][1],dictLinvar[pos+allele_1+allele_2][2]))
                    elif dictLinvar[pos+allele_1+allele_2][4]=="1" and dictLinvar[pos+allele_1+allele_2][5]=="0":
                            freqList_BENIGN.append((1-freqEXAC,freqCONVERGE,pos,allele_1,allele_2,REF_EXAC_COUNT,ALT_EXAC_COUNT,allele_1_count,allele_2_count,dictLinvar[pos+allele_1+allele_2][0],dictLinvar[pos+allele_1+allele_2][1],dictLinvar[pos+allele_1+allele_2][2]))
                




print args.chr,len(freqList_notEXAC_CONVERGE)

print freqList_notEXAC_CONVERGE

print "-----"


#PATHOGENIC
out=open(args.out+".PATHOGENIC","w")



out.write("chr,freqEXAC,freqCONVERGE,pos,allele_1,allele_2,REF_EXAC_COUNT,ALT_EXAC_COUNT,allele_1_count,allele_2_count")
out.write("\n")

for i in freqList:
    out.write(args.chr+","+str(i[0])+","+str(i[1])+","+str(i[2])+","+str(i[3])+","+str(i[4])+","+str(i[5])+","+str(i[6])+","+str(i[7])+","+str(i[8])+","+str(i[9])+","+str(i[10])+","+str(i[11])   )
    out.write("\n")



out.close()



#BENIGN

out=open(args.out+".BENIGN","w")


#header
out.write("chr,freqEXAC,freqCONVERGE,pos,allele_1,allele_2,REF_EXAC_COUNT,ALT_EXAC_COUNT,allele_1_count,allele_2_count")
out.write("\n")

for i in freqList_BENIGN:
    out.write(args.chr+","+str(i[0])+","+str(i[1])+","+str(i[2])+","+str(i[3])+","+str(i[4])+","+str(i[5])+","+str(i[6])+","+str(i[7])+","+str(i[8])+","+str(i[9])+","+str(i[10])+","+str(i[11])   )
    out.write("\n")



out.close()




#freqList_notEXAC_CONVERGE

out=open(args.out+".notEXAC_CONVERGE","w")


#header
out.write("chr,freqCONVERGE,pos,allele_1,allele_2,allele_1_count,allele_2_count")
out.write("\n")

for i in freqList_notEXAC_CONVERGE:
    out.write(args.chr+","+str(i[0])+","+str(i[1])+","+str(i[2])+","+str(i[3])+","+str(i[4])+","+str(i[5])+","+str(i[6])+","+str(i[7])+","+str(i[8])  )
    out.write("\n")



out.close()

print "Done!"







