import re
import os
y = open('RRM_list_uniprot_3.txt','r')
v = open("queries.txt","w")

RRM_list=[]
for line in y:
	q = line.split(' ')
	RRM_list.append(q[0])
	v.write(q[0])
	v.write('\n')
y.close()

#so files can be read by Xplor
extrastuff = ["<",">","DirEntry","'"]
r = re.compile(r"\b|\b".join(extrastuff))#bug somewhere in here telling it to stop after the first thing that is not a hit


directory = '/home/tfarg/pdb_zips/alphafold/alphafold_in'

 
#create a list of pdb files scanned 
a = open("all_files.txt","w")
z = open("RRM_files.txt","w")

for filename in os.scandir(directory):
	#list of all files
	b = (r.sub("", str(filename))) 
	c = b.split('>')
	d = c[0].split("'")
	a.write('%s'%d[1])
	a.write('\n')	
	e = d[1].split("AF-")
	f = e[1].split("-F1-")
	v.write(f[0])
	v.write('\n')
	if f[0] in RRM_list:
		z.write(f[0])
		z.write('\n')
		#create calcSA files
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CB"'), str('-o %s_CB_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CG"'), str('-o %s_CG_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CD1"'), str('-o %s_CD1_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CD2"'), str('-o %s_CD2_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name ND1"'), str('-o %s_ND1_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name OH"'), str('-o %s_OH_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CE1"'), str('-o %s_CE1_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CE2"'), str('-o %s_CE2_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CE3"'), str('-o %s_CE3_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CZ"'), str('-o %s_CZ_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CH2"'), str('-o %s_CH2_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CZ2"'), str('-o %s_CZ2_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name CZ3"'), str('-o %s_CZ3_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name NE1"'), str('-o %s_NE1_sa.txt'%d[1])))
		os.system ('calcSA {} {} {}' .format(str('%s'%d[1]), str('-reportSelection="name NE2"'), str('-o %s_NE2_sa.txt'%d[1])))

a.close()
z.close()
