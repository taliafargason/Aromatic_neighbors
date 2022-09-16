# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import re
domains = open('RRM_any_manual_assertion_reviewed_only_organism9606_existence1_fragmentfalse_reviewedtrue_2022.09.03.txt','r')
domain_list = []
for line in domains:
    a = re.split('\s',line)
    domain_list.append(a)
out = open('domains.txt','w')

for line in domain_list:
    #out.write(str(line))
    if len(line)>2:
        b = line[2].split('DOMAIN')
        for i in b:
            if len(i)>1:
                c = i.split('"')
                d = c[0].split('..')
                e = d[1].split(';')
                out.write("%s %s %s %s %s\n"%(line[0],d[0],e[0],c[1],line[3]))

    
out.close()