
import csv
import sys
import argparse
from collections import Counter
import numpy


ap = argparse.ArgumentParser()
ap.add_argument('inCLINVAR', help='CLINVAR file downloaded from https://github.com/macarthur-lab/clinvar')
ap.add_argument('inCONVERGE', help='CLINVAR file downloaded from https://github.com/macarthur-lab/clinvar')
ap.add_argument('inEXAC', help='CLINVAR file downloaded from https://github.com/macarthur-lab/clinvar')
ap.add_argument('chr', help='chr')
ap.add_argument('out', help='out')


#ap.add_argument('--m', action='store_true',help='Save multi-mapped reads')

#cmd https://gist.github.com/daler/ec481811a44b3aa469f3

args = ap.parse_args()



# Open CONVERGE

#1       10388   2       21280   A:21126 T:154


with open(args.inCONVERGE,'r') as f:
    reader=csv.reader(f,delimiter="\t")
    reader.next()
    for line in reader:





pathogenic_list=[]
benign_list=[]
conflicted_list=[]


clinvarSet=set()






dict={}
    
with open(args.inCLINVAR,'r') as f:
    reader=csv.reader(f,delimiter="\t")
    reader.next()
    for line in reader:
        
        chr=line[0]
        pos=line[1]
        allele=""
        if chr==args.chr:
            clinvarSet.add(pos)
            dict[pos]=line


        if line[4]=="ALT":
            allele=line[3]
        else:
            allele=line[2]
        
        
        
        
        
        
        pathogenic=line[8]
        benign=line[9]
        conflicted=line[10]
        
        
        
        
        
        pathogenic_list.append(pathogenic)
        benign_list.append(benign)
        conflicted_list.append(conflicted)




print len(clinvarSet)

print '97564154' in clinvarSet


#print Counter(pathogenic_list)
#print Counter(benign_list)
#print Counter(conflicted_list)


#--------

pathogenicCONVERGE=[]
pathogenicLikelyCONVERGE=[]
benignCONVERGE=[]
conflicted=[]

dict_freq_fullInfo={}


#1	55464869	G	T	ALT	4389	BSND	Pathogenic	no assertion criteria provided	0	NM_057176.2:c.10G>T	NP_476517.1:p.Glu4Ter	OMIM	Sensorineural deafness with mild renal dysfunction;SENSORINEURAL DEAFNESS WITH MILD RENAL DYSFUNCTION	19646679				germline	MedGen:C2748440


#

dictRSID_CLINVAR={}

with open("/u/home/s/serghei/project/CONVERGE/clinvar/clinvar_20170130.vcf") as f:
    reader=csv.reader(f,delimiter="\t")
    for line in reader:
        if "#" not in line[0]:
            if line[0]==args.chr:
                allele_1=line[3]
                allele_2=line[4]
                pos=line[1]
                rsid=line[2]
                dictRSID_CLINVAR[pos]=rsid



dictRSID_CONVERGE={}

with open(args.legend) as f:
    reader=csv.reader(f,delimiter="\t")
    reader.next()
    for line in reader:
        
        rsid=line[2]
        pos=line[0]
        dictRSID_CONVERGE[pos]=rsid


print "Open freq file"
#1       10388   2       21280   A:21126 T:154


with open(args.inCONVERGE,'r') as f:
    reader=csv.reader(f,delimiter="\t")
    reader.next()
    for line in reader:
        pos=line[1]
        
        
        if pos in clinvarSet:
            
            #print line,pos
            #print pos, dict[pos]
            
            pathogenic=dict[pos][8]
            benign=dict[pos][9]
            conflicted=dict[pos][10]
            
            status=dict[pos][7]
            
            allele=dict[pos][4]
            
            
            if pathogenic=='1' and conflicted=='0' and status!="Likely pathogenic":
                
                
                if allele=="ALT":
                    
                    rsidCLINVAR=dictRSID_CLINVAR[pos]
                    rsidCONVERGE=dictRSID_CONVERGE[pos]
                    
                    
                    aleleCount=int(line[5].split(":")[1])
                    if aleleCount!=0:
                    
                        if rsidCLINVAR ==rsidCONVERGE:
                            freq=int(line[5].split(":")[1])/21280.0
                            pathogenicCONVERGE.append(freq)
                            dict_freq_fullInfo[freq]=(line,dict[pos])
                        #else:
                            #print rsidCLINVAR,rsidCONVERGE


                else:
                    print line
                    print allele
                    print "EXIT"
                    sys.exit(1)
        

            elif pathogenic=='1' and conflicted=='0' and status=="Likely pathogenic":
                
                
                if allele=="ALT":
                    
                    rsidCLINVAR=dictRSID_CLINVAR[pos]
                    rsidCONVERGE=dictRSID_CONVERGE[pos]
                    
                    
                    aleleCount=int(line[5].split(":")[1])
                    if aleleCount!=0:
                        
                        if rsidCLINVAR ==rsidCONVERGE:
                            freq=int(line[5].split(":")[1])/21280.0
                            pathogenicCONVERGE.append(freq)
                            dict_freq_fullInfo[freq]=(line,dict[pos])
                        #else:
                            #print rsidCLINVAR,rsidCONVERGE
            
            
                else:
                    print line
                    print allele
                    print "EXIT"
                    sys.exit(1)
                
                


fileOut_pathogenicCONVERGE=open(args.out+".pathogenicCONVERGE","w")
fileOut_pathogenicLikelyCONVERGE=open(args.out+".pathogenicLikelyCONVERGE","w")

fileOut_pathogenicCONVERGE_summary=open(args.out+".pathogenicCONVERGE_summary","w")
fileOut_pathogenicLikelyCONVERGE_summary=open(args.out+".pathogenicLikelyCONVERGE_summary","w")
fileOut_benignCONVERGE_summary=open(args.out+".pathogenicLikelyCONVERGE_summary","w")



npathogenicCONVERGE_commom=0
npathogenicLikelyCONVERGE_commom=0
nbenignCONVERGE_commom=0



for i in pathogenicCONVERGE:
    if float(i) >0.01:
        pos=dict_freq_fullInfo[i][0][1]
        rsid=dictRSID_CLINVAR[pos]
        fileOut_pathogenicCONVERGE.write(str(i)+","+rsid+","+str(dict_freq_fullInfo[i][0][0])+","+str(dict_freq_fullInfo[i][0][1])+","+str(dict_freq_fullInfo[i][0][2])+","+str(dict_freq_fullInfo[i][0][3])+","+str(dict_freq_fullInfo[i][0][4])+","+str(dict_freq_fullInfo[i][0][5])+","+dict_freq_fullInfo[i][1][11]+","+dict_freq_fullInfo[i][1][16])
        fileOut_pathogenicCONVERGE.write("\n")
        
        npathogenicCONVERGE_commom+=1
        print i, dict_freq_fullInfo[i]


print "========="
print "Number of pathogenic in CONVERGE", len(pathogenicCONVERGE)
print "Mean freq", numpy.mean(pathogenicCONVERGE)
print "STD freq", numpy.std(pathogenicCONVERGE)
print "Number of common from CONVERGE which are pathogenic",npathogenicCONVERGE_commom

fileOut_pathogenicCONVERGE_summary.write(args.chr+","+str(len(pathogenicCONVERGE))+","+str(npathogenicCONVERGE_commom)+","+str(numpy.mean(pathogenicCONVERGE))+","+str(numpy.std(pathogenicCONVERGE)))
fileOut_pathogenicCONVERGE_summary.write("\n")
#===============



for i in pathogenicLikelyCONVERGE:
    if float(i) >0.01:
        os=dict_freq_fullInfo[i][0][1]
        rsid=dictRSID_CLINVAR[pos]
        fileOut_pathogenicCONVERGE.write(str(i)+","+rsid+","+str(dict_freq_fullInfo[i][0][0])+","+str(dict_freq_fullInfo[i][0][1])+","+str(dict_freq_fullInfo[i][0][2])+","+str(dict_freq_fullInfo[i][0][3])+","+str(dict_freq_fullInfo[i][0][4])+","+str(dict_freq_fullInfo[i][0][5])+","+dict_freq_fullInfo[i][1][11]+","+dict_freq_fullInfo[i][1][16])

        fileOut_pathogenicLikelyCONVERGE.write("\n")

        print i, dict_freq_fullInfo[i]
        npathogenicLikelyCONVERGE_commom+=1




print "========="
print "Number of LIKELY pathogenic in CONVERGE", len(pathogenicLikelyCONVERGE)
print "Mean freq", numpy.mean(pathogenicLikelyCONVERGE)
print "STD freq", numpy.std(pathogenicLikelyCONVERGE)
print "Number of common from CONVERGE which are LIKELY pathogenic",npathogenicLikelyCONVERGE_commom

fileOut_pathogenicLikelyCONVERGE_summary.write(args.chr+","+str(len(pathogenicLikelyCONVERGE))+","+str(npathogenicLikelyCONVERGE_commom)+","+str(numpy.mean(pathogenicLikelyCONVERGE))+","+str(numpy.std(pathogenicLikelyCONVERGE)))
fileOut_pathogenicLikelyCONVERGE_summary.write("\n")
#===================


sys.exit(1)

for i in pathogenicLikelyCONVERGE:
    if float(i) >0.01:
        print i
        nbenignCONVERGE_commom+=1

print "Number of benign >1%", nbenignCONVERGE_commom


fileOut_pathogenicCONVERGE.close()

print "========="
print "Number of benign in CONVERGE", len(benignCONVERGE)
print "Mean freq", numpy.mean(benignCONVERGE)
print "STD freq", numpy.std(benignCONVERGE)
print "Number of common from CONVERGE which are benign",nbenignCONVERGE_commom

fileOut_benignCONVERGE_summary.write(args.chr+","+str(len(benignCONVERGE))+","+str(nbenignCONVERGE_commom)+","+str(numpy.mean(benignCONVERGE))+","+str(numpy.std(benignCONVERGE)))
fileOut_benignCONVERGE_summary.write("\n")


