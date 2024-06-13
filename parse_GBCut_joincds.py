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



def GetGBRecords(path, fInName):

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
    
        
    
    return t_record

def GetField(t,i,s='"\n'):
    if i <0:
        return []
    j = -1
    #print("----------")
    #print(t,i, j)
    #print()
    for it in t[i:]:
        #print(it)
        if s in it:
            j = t.index(it)
            break
    #print(j)    
    return t[i:j+1]


def ExtractField(t):
    st = ""
    
    if t != []:
        n = len(t)
        
        #print(t)
        #print(n)
    
        
        if n == 1:
            ln0 = t[0]
            t_ln0 = ln0.strip().split('=')
            st = t_ln0[-1][1:-1]
        
        else:
            l1 = t[0]
            l1 = l1.strip()
            #print(l1)
            st = l1.split("=")[-1][1:]
            
            #print(st)
            
            for ln in t[1:]:
                ln = ln.strip()
                #print(ln)
                t_ln = ln.strip().split()
                if len(t_ln) > 1:
                    st = st+" "+ln
                else:
                    st = st+ln
            
            st = st[:-1]
             
    return st




#========================================================


#path = "/Users/tatianalenskaia/_Data/PAT_755phs/efetch_755phs_seq1/"
#fInName = "755phs_seq1.gb"


path = "/Users/tatianalenskaia/_Data/PAT_755phs/BatchEntrez_755phs/"
fInName = "sequence.gb"
#subdir = "755phs_gb_cds/"
subdir = ""


'''

# Example: 97 genes and 97 cds, no extra features
path = "/Users/tatianalenskaia/_Data/PAT_755phs/"
fInName = "NC_004166.2.gb"
subdir = ""


path = "/Users/tatianalenskaia/_Data/PAT_755phs/BatchEntrez_755phs/755phs_gb/"
fInName = "NC_001416.1.gb"



path = "/Users/tatianalenskaia/_Data/PAT_755phs/"
fInName = "NC_009799.3.gb"
subdir = ""

'''



label = "log"
#CheckGB_complete(path, fInName, label)
t_record = GetGBRecords(path, fInName)

print("Number of gb records:", len(t_record))

sep = "\t"

fOut = open(path+"_features.txt","w")
print("acc","genes","cds","(ct-1)/2", sep = sep, file = fOut)








fOut3 = open(path+fInName.split(".")[0]+"_cds.txt","w")
print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut3)



fOut4 = open(path+fInName.split(".")[0]+"_cds.faa","w")

for item in t_record:
    t_lc = item[0].strip().split()
    #print(t_lc)
    g_acc = t_lc[1]
    g_sz = int(t_lc[2])
    g_base = t_lc[4]
    g_type = t_lc[5]
    
    
    j = -1
    j_f = -1
    j_v = -1
    j_s = -1
    for ln in item:
        

        if ln[0:7] == "VERSION":
            j_v = item.index(ln)
        
        if ln[0:6] == "SOURCE":
            j_s = item.index(ln)

        if ln[0:8] == "FEATURES":
            j_f = item.index(ln)
        
        
        if "ORIGIN" in ln:
            j = item.index(ln)
            #print(item[j-1:j+2])
            break


    g_ver = item[j_v].strip().split()[-1]
    g_source = " ".join(item[j_s].strip().split()[1:])


    t = item[j_f+1:j]
    
    t_features = []
    tt = []
    for line in t:
        if line[5] != " ":
            #print(line.strip())
            if tt != []:
                t_features.append(tt)
                # The bug is caught (prev: tt = [line])
                tt = []
        tt.append(line)
    
    
    if tt != []:
        t_features.append(tt)
        
    
    #print(len(t_features))  
    #print(t_features[-1])
    
    
    ct_cds = 0
    ct_gene = 0
    
    ct = 0
    
    t_cds_index = []
    t_cds = []
    
    
    for it in t_features:
        ln = it[0]
        #print(ln.strip())
        ct += 1

        if ln.strip()[0:4] == "gene":
            ct_gene += 1
            #print(ln.strip())

        if ln.strip()[0:3] == "CDS":
            ct_cds += 1
            
            t_cds_index.append(t_features.index(it))
            t_cds.append(it)
            
            # Check if the last field in cds feature is always translation
            # Not always, exception: NC_009799.3.gb, CDS             9713..14527 /db_xref="GeneID:5580369"
            i_l = -1
            for l in it:
                if l.strip()[0] == "/":
                    i_l = it.index(l)
                    
            if it[i_l].strip()[0:13] != "/translation=":
                print(g_acc, ln.strip(), it[i_l].strip())
            #print(ln.strip())

    '''
    print("Number of genes:", ct_gene)
    print("Number of CDS:", ct_cds)
    
    print((ct-1) / 2)
    '''
    print(g_ver, ct_gene, ct_cds, (ct-1)/2, sep = sep, file = fOut)
    
 
    
 
    # Processing of CDS list   
 
    #print(len(t_cds_index))    
    #print(len(t_cds))




    #fOut1 = open(path+subdir+g_ver+"_cds.txt","w")
    #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut1)


    #fOut2 = open(path+subdir+g_ver+"_cds.faa","w")
    #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut1)



    j = 0



    for it in t_cds:
        #print(t_cds.index(it)+1, it[0].strip())
        
        
        j_pr = -1
        j_pid = -1
        j_tr = -1
        
        cds_loc = it[0].strip().split()[-1]
        
       
        #print(cds_loc)  


        for ln in it:


            if ln.strip()[0:9] == "/product=":
                j_pr = it.index(ln)
                #print(ln.strip())
            if ln.strip()[0:12] == "/protein_id=":
                j_pid = it.index(ln)
                #print(ln.strip())
            if ln.strip()[0:13] == "/translation=":
                j_tr = it.index(ln)
                #print(ln.strip())
                

            
            
            
        t_pr = GetField(it,j_pr)
        st_pr = ExtractField(t_pr)
        
        
        t_pid = GetField(it,j_pid)
        st_pid = ExtractField(t_pid)
        
        
        t_tr = GetField(it,j_tr)
        st_tr = ExtractField(t_tr)
        
        
        '''   
        if cds_loc == "9713..14527":
            print(j_tr)
            print(t_tr)
            print(it[0], it[1])
            break
        '''
                      
        
        #print(st_pr, st_pid, st_tr)
        
        
        print(g_ver, st_pid, st_pr, cds_loc, st_tr, len(st_tr), sep = sep, file=fOut3)
        
        
        sp = "|"
        header = ">"+g_ver+sp+st_pid +sp + st_pr +sp+ cds_loc +sp+str(len(st_tr))+ " aa"
        
        if st_tr != "":
            fOut4.write(header+"\n")
            fOut4.write(st_tr+"\n")
        else:
            print("missing translation!", g_ver, st_pid )
        
        
        
        '''
        if len(t_tr) > 1:
            #print(t_record.index(item), item[0], t_pr)
            print(t_tr)
         ''' 


        
        
        #print(t_pr)
        #break

    '''
    for i in t_cds_index:
        print(t_features[i][0].strip())

    '''
    #fOut1.close()
    #fOut2.close()
    #break
    
    
    

fOut3.close()
fOut4.close()

fOut.close()

    


#===============================================================




"""

#print(t_record)

ct_cds = 0
ct_gene = 0
ct_product = 0
ct_pid = 0
ct_tr = 0
for item in t_record:
    for ln in item:
        #if "CDS" in ln:
            
        if ln.strip()[0:4] == "gene":
            ct_gene += 1

        if ln.strip()[0:3] == "CDS":
            ct_cds += 1
            print(ln.strip())
        if ln.strip()[0:9] == "/product=":
            ct_product += 1
            print(ln.strip())
        if ln.strip()[0:12] == "/protein_id=":
            ct_pid += 1
            print(ln.strip())
        if ln.strip()[0:13] == "/translation=":
            ct_tr += 1
            print(ln.strip())
            
            
            

            

print("Number of genes:", ct_gene)
print("Number of CDS:", ct_cds)
print("Number of products:", ct_product)

"""











"""
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

"""


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
