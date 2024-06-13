#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 18:19:48 2023

@author: tatianalenskaia
"""

import GB_lib as gbl
import core_methods as cm


path = "/Users/tatianalenskaia/HMM/Results/From Sarah/Re_ Streptomyces phage genome/"

fInName = "flist_updated.txt"

t_list = cm.GetListFromFile(path+fInName)


print(t_list)



for f in t_list:
    #t_record = gbl.GetGBRecords(path+f)
    #print(t_record)
    
    gbl.SeqDB_joinCDS(f, path)




t_new = []

for f in t_list:
    f = f.split(".")[0]+"_cds.faa"
    t_new.append(f)
    
    
print(t_new)

fOutName = path+"cds_list_updated.txt"
fOut = open(fOutName, "w")
for it in t_new:
    fOut.write(it+"\n")
fOut.close()


cm.ConcatFiles(fOutName, "phiScoe_all12_updated_cds.faa",path)
    
