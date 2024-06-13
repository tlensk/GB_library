#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 20:28:38 2023

@author: tatianalenskaia
"""

import core_methods as cm


def CheckGB_complete(path, fInName, label):
    
    
    fIn = open(path+fInName,"r")
    lines = fIn.readlines()
    fIn.close()
    
    n = len(lines)

    
    
    
    # Count LOCUS entries
    fOut = open(path+label+"_LOCUS.txt", "w")
    ct = 0
    for i in range(n):
        line = lines[i]
        
        if line[0:5] == "LOCUS":
        #if line == "//\n":
            ct += 1   
            print(ct, i, line.strip(), sep = "\t", file = fOut)
    fOut.close()





    # Count // entries
    fOut = open(path+label+"_end.txt", "w")
    ct = 0
    for i in range(n):
        line = lines[i]
        

        if line == "//\n":
            ct += 1
            print(ct, i, line.strip(), sep = "\t", file = fOut)
            
            '''
            st = ""
            if lines[min(i+1,n-1)] == "\n":
                st = lines[min(i+2,n-1)]
            else:
                st = lines[min(i+1,n-1)]
                
            print(ct, i, st.strip(), sep = "\t", file = fOut)
            '''
    fOut.close()
    
    
    
    
    # Count ORIGIN entries
    
    fOut = open(path+label+"_ORIGIN", "w")
    ct = 0
    for i in range(n):
        line = lines[i]
        
        if "ORIGIN" in line:
            ct += 1
            
            print(ct, i, line.strip(), sep = "\t", file = fOut)
    fOut.close()
    
    return


#========================================================




#path = "/Users/tatianalenskaia/_Data/PAT_755phs/efetch_755phs_seq1/"
path = "/Users/tatianalenskaia/_Data/PAT_755phs/BatchEntrez_755phs/"

#fInName = "755phs_seq1.gb"
fInName = "sequence.gb"

label = "log"


#CheckGB_complete(path, fInName, label)



tt = []

t_record = []

fIn = open(path+fInName, "r")
lines = fIn.readlines()
fIn.close()



for line in lines:
    #print(line)
    if line != "\n":
        tt.append(line)
    if line == "//\n":
        if tt != []:
            #print(len(tt))
            t_record.append(tt)
            tt = []

    

print(len(t_record))


# --------------- Cut sequence.gb into individual gb files


sep = "\t"

subdir = "755phs_gb/"

fLog = open(path+subdir+"flist.txt","w")


for item in t_record:
    t_lc = item[0].strip().split()
    #print(t_lc)
    g_acc = t_lc[1]
    

    j_v = -1
    for ln in item:
        if ln[0:7] == "VERSION":
            j_v = item.index(ln)
            break

    
    g_ver = item[j_v].strip().split()[-1]
 
    
    fLog.write(g_ver+"\n")
    
    fOut = open(path+subdir+g_ver+".gb","w")
    for it in item:
        fOut.write(it)
    fOut.close()
    #break
fLog.close()




#------------------- end of cut sequence.db-------------




"""

# ----------------------Record output for sequence size check

sep = "\t"

fOutName = fInName.split(".")[0]+"_sizecheck.txt"
fOut = open(path+fOutName, "w") 
print("id","source","rec_size","gbase","gtype","download_size","size_mismatch", sep = sep, file = fOut)

c_mm = 0
for item in t_record:
    t_lc = item[0].strip().split()
    #print(t_lc)
    g_acc = t_lc[1]
    g_sz = int(t_lc[2])
    g_base = t_lc[4]
    g_type = t_lc[5]
    
    
    
    j = -1
    j_v = -1
    j_s = -1
    for ln in item:
        if ln[0:7] == "VERSION":
            j_v = item.index(ln)
        
        if ln[0:6] == "SOURCE":
            j_s = item.index(ln)
        
        if "ORIGIN" in ln:
            j = item.index(ln)
            #print(item[j-1:j+2])
            break
    t_seq = item[j+1:-1]
    #print(t_seq[0]+t_seq[-1])
    
    g_ver = item[j_v].strip().split()[-1]
    g_source = " ".join(item[j_s].strip().split()[1:])
    
    #print(g_source)
    
    #print(g_acc, g_ver)

    seq = ""
    for l in t_seq:
        t_l = l.strip().split()
        st = "".join(t_l[1:])
        #print(l)
        #print(st)
        seq = seq+st
    n_seq = len(seq)
    if n_seq != g_sz:
        fl = 1
        c_mm += 1
    else:
        fl = 0
    #print(g_sz, len(seq))    
    print(g_ver, g_source, g_sz,g_base, g_type,n_seq,fl, sep = sep, file = fOut)
    #print(item[0].strip().split()[1],item[j].strip(),item[-1].strip())
    #break    

fOut.close()

print("Size mismatch:", c_mm)

#---------------------------------------------- end of sequence size check
"""
