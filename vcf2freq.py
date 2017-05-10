import pysam
import argparse
import csv
import sys

#updated March 8, 2017


ap = argparse.ArgumentParser()
ap.add_argument('inVCF', help='vcf from exac')
ap.add_argument('out', help='')
ap.add_argument('chr', help='')
ap.add_argument('type', help='NFE or EAS')

ap.add_argument('n', help='Number of individuals the varinat to be present in')



args = ap.parse_args()



if args.type!="NFE" and args.type!="EAS":
    print "ERROR : type needs to be NFE or EAS"
    sys.exit(1)



out=open(args.out,"w")


#1       13372   .       G       C
#------
print "Read variants from ExAC vcf ...", args.inVCF
number3rdAlele=0




numberFiltered10K=0

with open(args.inVCF) as csvfile:
    readCSV = csv.reader(csvfile, delimiter='\t')
    for row in readCSV:
        if row!=[]:
            if "#" not in row[0] and row[0]==args.chr:
                
                
                
                
                
                
                    
                    
                    
                if len(row[3])==1 and len(row[4])==1:
                    
                        AC=int(row[7].split("AC_"+args.type)[1].split(";")[0].split("=")[1])
                        AN=int(row[7].split("AN_"+args.type)[1].split(";")[0].split("=")[1])
                        

                        
                      
                        
                        if (AN/2)>int(args.n):
                            
                            
                            
                            if row[1]=="45409472":
                                print "1:",row[1],row[3],row[4]
                        

                            
                            
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
                    
                        number3rdAlele+=1
                        
                        
                        
                        acList=[]
                        acList[:]=[]
                        
                        anList=[]
                        anList[:]=[]
                        
                        
                        
                        
                        for i in row[3].split(","):
                            for j in row[4].split(","):
                                acList.append(i)
                                anList.append(j)
                    
                        k=0
                        for i in row[7].split("AN_"+args.type+"=")[1].split(";")[0].split(","):
                            for j in row[7].split("AC_"+args.type+"=")[1].split(";")[0].split(","):
                                i=int(i)
                                j=int(j)
                                if i==0:
                                    out.write(row[0]+","+row[1]+","+row[2]+","+acList[k]+","+anList[k]+","+str(i)+","+str(j)+",0.0")

                                else:
                                    out.write(row[0]+","+row[1]+","+row[2]+","+acList[k]+","+anList[k]+","+str(i)+","+str(j)+","+str(1-float(int(j))/int(i)) )
                                    out.write("\n")
                                k+=1

                        

                #AN_FIN=int(row[7].split(";")[17].split("=")[1])
#AN_NFE=int(row[7].split(";")[18].split("=")[1])
                
                
                    #print row[7].split(";")[17],row[7].split(";")[18]

                    #print AC_FIN,AC_NFE,AN_FIN,AN_NFE
                    #sys.exit(1)




print "Number of variants with 3rd allele,", number3rdAlele
print "Number of variants with less then specified threshold,", numberFiltered10K





