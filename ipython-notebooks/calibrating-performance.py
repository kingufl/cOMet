import subprocess as sp
import matplotlib
import Bio.SeqIO
import Bio.Restriction
import scipy
#import matplotlib.pyplot as plt
import math
import random
import numpy as np
from timeit import timeit
import time

def FindSandTScores(infile_e, infile_w_e, loc = "/s/oak/b/nobackup/darshanw/valuev_optmap_alignment-master/ovlp/"):     
    prg = loc + "ovlp"
    output = loc+ "TempOutput"
    lines_e = [line for line in open(infile_e)]
    lines_w_e = [line for line in open(infile_w_e)]
    temp_file = loc + "temp"

    Naligned = 0
    Sscores = []
    Tscores = []
    #for i in range(0, 2):
    for i in range(0, len(lines_e)/3):
        t = open(temp_file, 'w')
        for j in range(0,3):
            t.write(lines_e[i*3+j])
        for j in range(0,3):
            t.write(lines_w_e[i*3+j])
        t.close()
        rvalue = sp.check_output([prg, temp_file, output, "0", "1"], stderr=sp.STDOUT)
        rvalue = rvalue.strip().split("\n")
        scores = rvalue[-1].strip().split(":")
        if(len(scores) > 2):
            Naligned += 1
            continue
        else:
            Tscores.append([float(scores[1]), i])
            scores = rvalue[-2].strip().split(":")
            Sscores.append([float(scores[1]), i])
    
    return([Sscores, Tscores, Naligned])

def getCorrespondingValAndPlot(A, B, k, m, d, path, output):
    out = []
    if(len(A) == 0 or len(B) == 0):
        return(out);
    Avalues = list(zip(*A)[0])
    Aindex = list(zip(*A)[1])
    Bvalues = list(zip(*B)[0])
    Bindex = list(zip(*B)[1])
    
    for i in range(0, len(Bindex)):
        try:
            indx = Aindex.index(Bindex[i])
        except:
            continue
        out.append([Avalues[indx], Bvalues[i], Bindex[i]])

    #plt.plot(zip(*out)[0],zip(*out)[1], "*")
    #plt.plot([0, 400], [0, 400], 'r-')
    #plt.xlabel("S-score before correction")
    #plt.ylabel("S-score after correction")
    #plt.savefig((path+'k'+str(k)+'m'+str(m)+'d'+str(d)+'.png'), bbox_inches='tight')
    #plt.close()

    improv = []
    for v in out:
        improv.append(v[1]/v[0])
    output.write("S-score is improved for: " + str(len([v for v in out if v[1]>v[0]])) + " out of "+str(len(out))+" reads\n")
    #plt.hist(improv, 10)
    #plt.xlabel("Ratio of S-score after and before correction")
    #plt.ylabel("Frequency")
    #plt.tight_layout()    
    #plt.savefig((path+'k'+str(k)+'m'+str(m)+'d'+str(d)+'-hist'+'.png'), bbox_inches='tight')
    #plt.close()

def findTrueErrorLocation(E):
    deletion = E.strip().split(":")[0].strip().split(",")
    insertion = E.strip().split(":")[1].strip().split(",")
    insertion = [ (int(insertion[i])+i) for i in range(0, len(insertion))]
    deletion = [ (int(deletion[i])-i) for i in range(0, len(deletion))]
    for z in range(0,len(deletion)):
        less = [i for i in insertion if i < deletion[z]]
        deletion[z] = deletion[z] + len(less)        
    return(deletion, insertion)

def howManyCorrected(dactual, iactual, dcorrected, icorrected):
    correctlyCorrectedDeletions = 0
    correctlyCorrectedInsertions = 0
    incorrectlyCorrectedDeletions = 0
    incorrectlyCorrecteInsertions = 0
    for d in dcorrected:
        for delloc in d[1]:
            if delloc in dactual[d[0]][1]:
                correctlyCorrectedDeletions+=1
            else:
                incorrectlyCorrectedDeletions+=1                
    for i in icorrected:
        for insloc in i[1]:
            if insloc in iactual[i[0]][1]:
                correctlyCorrectedInsertions+=1
            else:
                incorrectlyCorrecteInsertions+=1                
    return(correctlyCorrectedDeletions,correctlyCorrectedInsertions,incorrectlyCorrectedDeletions,incorrectlyCorrecteInsertions)

def findInsertionDeletion(factual, fcorrected, output):
    dactual = []
    iactual = []
    dcorrected = []
    icorrected = []
    deletions = 0
    insertions = 0
    lines = [line for line in open(factual)]
    for i in range(0, len(lines)):
        R = findTrueErrorLocation(lines[i])
        dactual.append([i, R[0]])
        iactual.append([i, R[1]])
    lines = [line for line in open(fcorrected)]
    dloc = []
    iloc = []    
    for i in range(0, len(lines)):        
        if(i%2 == 0):
            line = lines[i].strip().split(" ")
            dloc = []
            iloc = []
            for e in line:
                if(e[0] == '-'):                    
                    e = e[1:]
                    dloc.append(int(e))
                elif(e[0] == '+'):
                    e = e[1:]
                    iloc.append(int(e))
        else:            
            deletions += int(lines[i].strip().split(" ")[2])
            insertions += int(lines[i].strip().split(" ")[6])
            readNo = int(lines[i].strip().split(" ")[0])
            dcorrected.append([readNo, dloc])
            icorrected.append([readNo, iloc])
    
    ccd, cci, icd, ici = howManyCorrected(dactual, iactual, dcorrected, icorrected)
    output.write("Deletions: "+ str(deletions)+"\n")
    output.write("Insertions: "+ str(insertions)+"\n")
    output.write("Correctly-Corrected-Deletions: "+ str(ccd)+"\n")
    output.write("Incorrectly-Corrected-Deletions: "+ str(icd)+"\n")
    output.write("Correctly-Corrected-Insertion: "+ str(cci)+"\n")
    output.write("Incorrectly-Corrected-Insertion: "+ str(ici)+"\n")

def getAvg(S):
    add = 0.0
    for x in S:
        add=add+x[0]
    avg = add/len(S)
    return(avg)

def seperateFile(inf, outf, elocf):
    inlines = [line for line in open(inf)]
    outlines = [line for line in open(outf)]
    of = open(outf, "w")
    ef = open(elocf, "w")
    for i in range(0, (len(outlines)-len(inlines))):
        ef.write(outlines[i])
    for i in range((len(outlines)-len(inlines)), (len(outlines))):
        of.write(outlines[i])
    of.close()
    ef.close()

copies = 200
prg_path = "/s/oak/b/nobackup/darshanw/COmap/"
prg = prg_path + "bin/COmap"
in_file = prg_path + "test/sim_single_molecule_100_newDel"
eloc_file = prg_path + "test/sim_single_molecule_100_elocations"
e_free = prg_path + "test/sim_single_molecule_100_efree"
output = open(prg_path + "test/output/output.txt", 'w', 0)
#comb = [[2,2,2],[3,2,2],[4,2,2],[5,2,2],[6,2,2],[7,2,2],[8,2,2], [3,3,2], [3,4,2], [3,5,2], [3,6,2], [3,7,2], [3,8,2], [3,2,3], [3,2,4], [3,2,5], [3,2,6], [3,2,7], [3,2,8]]

comb = [[3,2,2]]
for k,m,d in comb:
    out_file = prg_path + "test/output/sim_single_molecule_100_corrected"
    out_file = out_file+"_k"+str(k)+"m"+str(m)+"d"+str(d)+"_"+str(copies)
    eloc_corrected = prg_path + "test/output/corrected"            
    eloc_corrected = eloc_corrected+"_k"+str(k)+"m"+str(m)+"d"+str(d)+"_"+str(copies)                        
    
    output.write("k: "+str(k)+" | m: "+str(m) + " | d : " + str(d)+ " > \n")            
    time = (timeit(stmt = "subprocess.call(['"+prg+"', '-x','-z', '-k', '"+str(k)+"', '-m', '"+str(m)+"', '-d', '"+str(d)+"', '-f','"+in_file+"'], stdout = open('"+out_file+"', 'wb'), stderr=subprocess.STDOUT)", setup = "import subprocess", number = 1))
    output.write("Time taken to execute: "+ str(time)+" sec\n")
    seperateFile(in_file, out_file, eloc_corrected)
    S1, T1, Na1 = FindSandTScores(in_file, e_free)          
    output.write("aligned without correction :" + str(len(S1)) + "("+ str(len(S1)+Na1) + ")\n")
    output.write("S-score average before correction: "+ str(getAvg(S1))+"\n")
    S2, T2, Na2 = FindSandTScores(out_file, e_free)
    output.write("aligned after correction :" + str(len(S2)) + "("+ str(len(S2)+Na2) + ")\n")
    output.write("S-score average after correction: "+ str(getAvg(S2))+"\n")
    getCorrespondingValAndPlot(S1,S2, k, m, d, (prg_path+'test/images/'), output)
    #sp.call([prg, "-x", "-k", str(k), "-m", str(m), "-d", str(d), "-f",in_file], stdout = open(eloc_corrected, "wb"), stderr=sp.STDOUT)            
    findInsertionDeletion(eloc_file, eloc_corrected, output)
output.close()
