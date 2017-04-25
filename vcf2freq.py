import pysam
import argparse
import csv
import sys

#updated March 8, 2017


ap = argparse.ArgumentParser()
ap.add_argument('inVCF', help='vcf from exac')
ap.add_argument('out', help='')
ap.add_argument('chr', help='')
ap.add_argument('i', help='EUR NFE 8,18 EAS 4,16')
ap.add_argument('j', help='')
ap.add_argument('n', help='Number of individuals the varinat to be present in')



args = ap.parse_args()


out=open(args.out,"w")


#1       13372   .       G       C
#------
print "Read variants from ExAC vcf ...", args.inVCF
number3rdAlele=0

i=int(args.i)
j=int(args.j)


numberFiltered10K=0

with open(args.inVCF) as csvfile:
    readCSV = csv.reader(csvfile, delimiter='\t')
    for row in readCSV:
        if row!=[]:
            if "#" not in row[0] and row[0]==args.chr:
                
                if not((len(row[3])>1 and row[3].count(",")==0) or  (len(row[4])>1 and row[4].count(",")==0)):
                    
                    
                    if len(row[3])==1 and len(row[4])==1:
                    
                    

                        AC=int(row[7].split(";")[i].split("=")[1])
                    
                        AN=int(row[7].split(";")[j].split("=")[1])
                        
                        
                        #print row
                        #print row[7].split(";")[i].split("=")[0], row[7].split(";")[j].split("=")[0]
                        
                        
                        
                        
                        if (AN/2)>int(args.n):
                            
                            
                            
                            
                            
                            
                            #print row, AN_FIN+AN_NFE
                            #MAF=float(AC_FIN+AC_NFE)/(AN_FIN+AN_NFE)
                            
                            MAF=float(AC)/(AN)
                            #print AC_NFE,AN_NFE,AC_NFE+AN_NFE, MAF
                        
                        
                            out.write(row[0]+","+row[1]+","+row[2]+","+row[3]+","+row[4]+","+str(AC)+","+str(AN)+","+str(1-MAF))
                            out.write("\n")
                
                        else:
                            numberFiltered10K+=1

#if AC_FIN.count(",")>0:
#                            print AC_FIN
#                            sys.exit(1)
                    else:
                        #print row,row[3],row[4],len(row[3]),len(row[4])
                        #sys.exit(1)
                        number3rdAlele+=1

                #AN_FIN=int(row[7].split(";")[17].split("=")[1])
#AN_NFE=int(row[7].split(";")[18].split("=")[1])
                
                
                    #print row[7].split(";")[17],row[7].split(";")[18]

                    #print AC_FIN,AC_NFE,AN_FIN,AN_NFE
                    #sys.exit(1)




print "Number of variants with 3rd allele,", number3rdAlele
print "Number of variants with less then specified threshold,", numberFiltered10K





