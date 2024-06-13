#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 13:12:43 2024

@author: tatianalenskaia
"""

import GB_lib as gbl






def CreateDictLocD_upper(text, m, gtp = "l"):
    if (m <= 0) or (m > len(text)):
        print("(n = "+str(m)+") is not a valid window size for this genome!!!");
        return {}
	
    d_g = dict()
    nn = len(text);
    gtype = gtp.lower();
    gtype = gtype[0];
	
    if gtype == "c":
        text = text + text[0:(m-1)];
        lastpos = nn;
    elif gtype == "l":
        lastpos = nn-m+1;
    else:
        print("Is this genome linear or circular?");
        return d_g;
		
    for ii in range (lastpos):
        bl = text[ii:(ii+m)]; 
        bl = bl.upper();
        if bl in d_g :
            dd = d_g[bl];
            dd[ii] = ""
            d_g[bl] = dd;
        else:
            dd = {}
            dd[ii] = ""
            d_g[bl] = dd;
    return d_g;	



def FindIntersection(d_g11, d_g22):
	
	t_ints = list()
	nn1 = len(d_g11);
	nn2 = len(d_g22);
	
	if nn1 <= nn2:
		d_first = d_g11;
		d_last = d_g22;
	else:
		d_first = d_g22;
		d_last = d_g11;
	
	for s in d_first:
		if s in d_last:
			t_ints.append(s);
	return t_ints;




def PrintMatrix(path, k, pref, row_index, col_index, matrix, sep = ","):
    fOut = open(path+str(k)+pref+".csv","w");
    fOut.write(str(k)+sep+sep.join(str(x) for x in col_index)+"\n")
    for i in range(len(matrix)):
        s = str(row_index[i])
        for el in matrix[i]:
            s = s+sep+str(el)
        fOut.write(s+"\n")
        
    fOut.close()
    return

def PrintSubMatrix(path, k, pref, row_names, row_index, col_names, col_index, matrix, sep = ","):
    if row_index != [] and col_index != []:
        fOut = open(path+str(k)+"_"+pref+".csv","w");
        #fOut.write(str(k)+sep+sep.join(str(x) for x in col_index)+"\n")
        s = str(k)
        for jj in range(len(col_index)):
            s = s+sep+col_names[col_index[jj]]
        fOut.write(s+"\n")

        
        for ii in range(len(row_index)):
            i = row_index[ii]
            s = str(row_names[row_index[ii]])
            
            for jj in range(len(col_index)):
                j = col_index[jj]
                
                if i == j:
                    sign = "-"
                else:
                    sign = ""
                    
                s = s+sep+sign+str(matrix[i][j])
            
            fOut.write(s+"\n")
            
        fOut.close()
    return







path = "/Users/tatianalenskaia/HMM/755phs/"


fInName = "755phs_genomes.fasta"



path = "/Users/tatianalenskaia/UT/Alexa/Alexa_227_borders/"
fInName = "227_genomes.fasta"
label = "_277"

fInName = "Pphs_899.fasta"
label = "_899"




res = gbl.GetFasta_header(path+fInName)

t = res[-1]

d_fa = res[0]

print(len(d_fa))

print(t[0],len(d_fa[t[0]][-1]))

m = 1000


#t = t[0:3]


#def ComputeSimilarityMatrix (pref, m, d_fa, t = []):
    

if t == []:
    t = list(d_fa.keys())


t = t[0:100]

print(t)

n = len(t)

print(t)
print(len(t))



matrix = [[0 for i in range(n)] for j in range(n)] 

for i in range(n):
    textA = d_fa[t[i]][-1]
    #print(fInNameA, len(textA))

    dictA =  CreateDictLocD_upper(textA, m)
    n_a = len(dictA)
    
    j = i-1
    print(i)
    while(j < n-1):
        j = j+1

        textB = d_fa[t[j]][-1]
        dictB =  CreateDictLocD_upper(textB, m)
        
        t_ints = FindIntersection(dictA,dictB)
        n_ints = len(t_ints)
        
        n_b = len(dictB)
        
        
        val = round(n_ints/min(n_a,n_b),2)
        
        print(n_a,n_b,n_ints,val)
        
        if n_ints != 0:
            #print(i,j,len(t_ints))
            matrix[i][j] = val
            matrix[j][i] = val
        #print(fInNameB, len(textB))
        #print(i,j)
        dictB.clear()
    #print("\n")
    dictA.clear()

#print(matrix)

PrintMatrix(path, m,label, t, t, matrix, sep = ",")


