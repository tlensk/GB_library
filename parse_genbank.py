#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 20:28:38 2023

@author: tatianalenskaia
"""

import core_methods as cm



def Parse_record(tt):
    if len(tt) > 1:
        ln = tt[0]
        t_ln = ln[1:].split(" ",1)
        nm = t_ln[0]
        
        
        t_ln = t_ln[1].split(" [")
        
        # ------Extract desc and loc
        
        i_pr = -1
        i_loc = -1
        
        for item in t_ln:
            if "location=" in item:
                i_loc = t_ln.index(item)
            if "protein=" in item:
                i_pr = t_ln.index(item)
        
        if i_loc >= 0:
            loc = t_ln[i_loc]
        else:
            loc = ""
            
        if i_pr >= 0:
            desc = t_ln[i_pr]
        else:
            desc = ""
        
        # ======end Extract desc and loc
        
        seq = ""
        for item in tt[1:]:
            seq = seq+item
                
    else:
        return []
    
    return [nm, seq, desc, loc]



path = "/Users/tatianalenskaia/HMM/755phs/"
fInName = "NC_004166.2.gb"

# Incomplete record for at least 1 phage
#fInName = "755phs_seq.gbk"
fInName = "755phs_seq.gb"




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

#print(t_record[-1][0]+t_record[-1][-1])


'''
i = 0
for it in t_record:
    i = i+1
    print(i, it[0])

'''

print(len(lines))
print(lines[-1])


n = len(lines)

ct = 0


fOut = open(path+"log_new_start.txt", "w")
for i in range(n):
    line = lines[i]
    
    if line[0:5] == "LOCUS":
    #if line == "//\n":
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














'''

p = 723801
p1 = 721991

#print(lines[p-100:p+1])


print(lines[p1-3:p1])

#print(t_record[255][0])
'''




"""







i = 0
tt = []


d = {}
dd = {}
# the keys in d dictionary should match the names of the sequences in the first column of the .out file

fIn = open(path+fInName, "r")
for line in fIn:
    line = line.strip()
    #print(line)
    i = i+1
    if line[0] == ">":
        if tt != []:
            # Change spot 1 of 2
            res = Parse_record(tt)
            
            if res != []:    
                nm = res[0]
                seq = res[1]
                desc = res[2]
                loc = res[3]
                
                
                ttt = nm[4:].split("_prot_")
                #print(ttt)
                
                gn = ttt[0]
                pt = ttt[1]
                
                if gn not in dd:
                    dd[gn] = []
                dd[gn].append([pt, seq, desc, loc])
                
                if nm not in d:
                    d[nm] = [seq,desc,loc]
                else:
                    print(nm, "duplicate")
            
        tt = []
        tt.append(line)
            
    else:
        tt.append(line)

    #break

if tt != []:
    # Change spot 2 of 2
    res = Parse_record(tt)
    
    if res != []:
        nm = res[0]
        seq = res[1]
        desc = res[2]
        loc = res[3]
        
        
        ttt = nm[4:].split("_prot_")
        #print(ttt)
        
        
        if gn not in dd:
            dd[gn] = []
        dd[gn].append([pt, seq, desc, loc])
        
        if nm not in d:
            d[nm] = [seq, desc, loc]
        else:
            print(nm, "duplicate")
                
print(len(d))


#print(nm, d[nm])

print(len(dd))


print(dd[gn][0])

sep = "\t"

fOutName = fInName.split(".")[0]+"_info.txt"
fOut = open(path+fOutName,"w")
print("Genome","Protein","Sequence", "Description","Location", sep = sep, file = fOut)



fOutName1 = fInName.split(".")[0]+"_stat.txt"
fOut1 = open(path+fOutName1,"w")
print("Genome","Number_of_proteins", sep = sep, file = fOut1)


for it in dd:
    tt = dd[it]
    print(it, len(tt), sep = sep, file = fOut1)
    for it in tt:
        print(it, it[0], it[1], it[2], it[3], sep = sep, file = fOut)
fOut.close()
fOut1.close()






fListName = "Rset_Custom_275.txt"

#fList = open(path+fListName, "r")

t_list = cm.GetListFromFile(path+fListName)
print(len(t_list))




#t_list = ["PF11123.11.out"]

#t_list = ["PF07278.14"]


t_list = ["aPAT_Pfam_ST.out"]

#path = "/Users/tatianalenskaia/HMM/755phs/"
#fInName = "U_Afp14-C.out"
#fInName = "U_A118_TT.out"





for it in t_list:
    
    #fInName = "PF11123.11.out"
    #fInName = it+".out"
    fInName = it
    
    
    fOutName = fInName.split(".")[0]+"_res.txt"
    fOut = open(path+"/parse_hmmsearch/res/"+fOutName,"w")
    
    
    
    
    
    #fOutName1 = fInName.split(".")[0]+"_seqs_e-50.faa"
    #fOutName1 = fInName.split(".")[0]+"_seqs.faa"
    fOutName1 = fInName.split(".")[0]+"_seqs.faa"
    fOut1 = open(path+"/parse_hmmsearch/faa/"+fOutName1,"w")
    
    
    i = 0
    c = 0
    
    d_hmm = {}
    
    fIn = open(path+fInName, "r")
    lines = fIn.readlines()
    for line in lines:
        
        if line[0] != "#":
            line = line.strip()
            t_line = line.split()
            ln = t_line[0]
            
            '''
            ln = "lcl|NC_004740.1_prot_NP_835581.1_65"
            print(ln,"!!!!!", d[ln])
            ln = "lcl|NC_001604.1_prot_NP_042007.1_54"
            print(ln,"!!!!!", d[ln])
            break
            '''
            
            #hmm = t_line[3]
            hmm_acc = t_line[4]
            hmm_name = t_line[3]
            
            e = t_line[6]
            '''
            t_ln = ln.split("##")
            f_nm = t_ln[0]
            s_nm = t_ln[1]
            '''
            p_st = int(t_line[17])-1
            p_fn = int(t_line[18])-1        
            i = i+1
            #print(ln,p_st, p_fn)
            
            
            st = ""
            for i in range(len(d[ln][0])):
                if (i >= p_st) and (i <= p_fn):
                    st = st+d[ln][0][i].upper()
                else:
                    st = st+d[ln][0][i].lower()
            #print(ln)
            
            if ln in d:
                #print(d[ln][p_st:p_fn+1])
                if hmm_acc not in d_hmm:
                    d_hmm[hmm_acc] = []
                #d_hmm[hmm_acc].append([hmm_name, ln, e, p_st,p_fn, d[ln][0][p_st:p_fn+1]])
                d_hmm[hmm_acc].append([hmm_name, ln, d[ln][1], d[ln][2], e, p_st,p_fn, st])
                
                
                e_f = float(e)
                
                if e_f < 100:
                    fOut1.write(">"+ln+" "+e+" "+str(p_st+1)+" "+str(p_fn+1)+" "+str(p_fn-p_st+1)+"\n")
                    #fOut1.write(d[ln]+"\n")
                    fOut1.write(d[ln][0][p_st:p_fn+1]+"\n")
                    c = c+1
                    

        
    fIn.close()

    print(fInName, c)
    
    
    #print(d_hmm)
    
    print("hmm_acc", "hmm_name","Genome_protein", "Protein_desc", "Protein_loc","E-value", "align_start", "align_stop","protein_sequence",file = fOut, sep = "\t" )
    
    for it in d_hmm:
        tt = d_hmm[it]
        for item in tt:
            print(it, item[0],item[1],item[2], item[3], item[4],item[5],item[6],item[7], file = fOut, sep = "\t" )
        
        
    fOut.close()
    
    
    
    fOut1.close()
    #break
"""