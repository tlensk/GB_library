#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 21:24:59 2024

@author: tatianalenskaia
"""

import GB_lib as gbl


path = "/Users/tatianalenskaia/UT/Alexa/227/"
fInName = "Alexa_227.gb"



path = "/Users/tatianalenskaia/UT/Vero/Vero_pr/"
fInName = "Vero_30phs.gb"


t_record = gbl.GetGBRecords(path+fInName)

t_fts = gbl.GetGBFeatures(path+fInName)

print(len(t_fts))

n = len(t_fts)




for j in range(n):
    print(j+1)
    print(t_record[j][0])
    print(t_fts[j][0])

    t_features = t_fts[j]
    t = []
    i_ft = -1
    n_ft = len(t_features)
    for i in range(n_ft):
        ft = t_features[i]
        for it in ft:
            if "integrase" in it:
                print(i)
                i_ft = i
                
    print(i_ft)
    if i_ft != -1:
        print(t_features[i_ft])
    print()



'''
for j in range(n):
    j = 2
    #print(t_features[0])
    t_features = t_fts[j]
    t = []
    i_ft = -1
    n_ft = len(t_features)
    for i in range(n_ft):
        ft = t_features[i]
        for it in ft:
            if "YP_001293398.1" in it:
                print(i)
                i_ft = i
                
    print(i_ft)
    print(t_features[i_ft-10:i_ft+12])
    break
'''
