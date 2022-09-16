#this script has two parts:
#1) read in sa file outputs from surface area calculation script and create a list of surface exposed aromatic RRM residues for each pdb structure
#2) iterate through individual pdb structures and determine whether aromatic residues possess neighboring positive/negative charges, out of these which are phase separated


import re
import os
import math


#directory = '/home/tfarg/pdb_zips/alphafold/distance_of_charges_from_sa_aromatics_test/'
directory = '/home/tfarg/pdb_zips/alphafold/F1_models_only/'

###Set thresholds
dist = 6 #angstroms from aromatic
thresh = 0.6 #threshold for available surface exposed (fraction)
neighbor_threshold = 1 #number of charged residues that count as a "cluster" of neighbors
His_SA_poss = 94.8697
Trp_SA_poss = 130.2387
Tyr_SA_poss = 119.7647
Phe_SA_poss = 108.6047

aromatics = ['PHE','TYR','TRP']#currently removing H to resolve charge concerns, may want to change to make Tyr a special case, easiest way to do this is to separate folders
negatives = ['ASP','GLU']#if want to count TYR and CYS as negatives, add them here
positives = ['ARG','LYS']#if want to count HIS as positive, add it here
#currently we want to know if something is close to the sidechain atoms of the aromatic residues, not just any atoms.  Therefore:
query_atom_type = ['NE1','NE2','ND1','CB','CG','CD1','CD2','CE1','CE2','CE3','CH2','CZ','CZ2','CZ3']#currently only including residues close to or in the aromatic group and not the tyrosine OH.Changing these atoms too much could affect the accuracy of the count if the code is not adjusted (see below)
subject_atom_type = ['N','NE','NE1','NE2','ND','ND1','ND2','NH','NH1','NH2','C','CA','CB','CG','CG1','CG2','CD','CD1','CD2','CE','CE1','CE2','CE3','CZ','CZ1','CZ2','CZ3','O','OD1','OD2','OE1','OE2','OXT','OG','OG1','OG2','OH']

#to correct for columns bumping into each other when number of residues in the protein exceeds 1000
chainID = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']#add letters if chainIDs are different.  May also have to adjust if in different locations.

##Read in parameters
first = open("RRM_list_uniprot_3.txt","r")

#phase separation databases
phasep = open("phasepdb2_MLOhtpb.txt","r")
DrLLPS = open("DrLLPS.txt","r")
llpsprotein = open('llps_proteinc.txt',encoding = 'windows-1252')#I don't know why but I had to change the encoding to get these lists to read in properly on python3.7.  Did not have this problem on python3.10
combo = open('combo.txt',encoding = 'windows-1252')

##Create outfiles
ab = open("by_residue_surface_area_6_0.6.txt","w")
ab.write("ID pdb resid resn SA fraction")
ab.write('\n')
z = open("test2b_6_0.6.txt","w")
qh = open('by_residue_neighbors_6_0.6.txt','w')
test4 = open('test4b_6_0.6.txt','w')
zqb = open('test5b_6_0.6.txt','w')
test7 = open('test7_b_6_0.6.txt','w')
totalxyz = open('overall_counts_6_0.6.txt','w')
simple_counts = open('simple_counts_any_residue_in_an_RRM_6_0.6.txt','w')
simple_counts.write("Counts of all residues within RRM domains regardless of whether they are surface exposed")

ID_ab = []
pdb_ab = []
resid_ab = []
resn_ab = []
SA_ab = []
fraction_ab = []
domain_ID_ab = []

rrmlist = []


for line in first:
    rrmlist.append(line.split(' '))

first.close()


phasep_plist = []
phasep_ident = []
for line in phasep:
    phasep_plist.append(line.split(' '))
phasep.close()


DrLLPs = []
vitro_prot = []
combined = []
a = []
b = []
c = []

for line in llpsprotein:
	result = re.split('\s+',line)
	a.append(result)
for line in combo:
	result = re.split('\s',line)
	b.append(result)
for line in DrLLPS:
	result = re.split('\s+',line)
	c.append(result)

for line in phasep_plist:
	if len(line)>1:
		phasep_ident.append(line[1])
for i in range(len(a)):
	vitro_prot.append(a[i][0])
for j in range(len(b)):
	combined.append(b[j][0])
for k in range(len(c)):
	DrLLPs.append(c[k][1])
llpsprotein.close()
combo.close()
DrLLPS.close()

###PART ONE: Read in surface area calculations

#start simple counts
total_proteins = 0
s_counts = []
SA_aromatic = 0
#so files can be read
extrastuff = ["<",">","DirEntry","'"]
r = re.compile(r"\b|\b".join(extrastuff))#bug somewhere in here telling it to stop after the first thing that is not a hit
 
for filename in os.scandir(directory):
	#list of all files
	b = (r.sub("", str(filename))) 
	c = b.split('>')
	d = c[0].split("'")#d[1] = sa file name
	split_string = re.split(r'_CB|_CD|_CG|_CE|_CH|_CZ|_ND|_NE|_OH',d[1])
	pdbfile = split_string[0]
	e = d[1].split("AF-")
	f = e[1].split("-F")#f[0] = uniprot ID 
	#f = e[1].split("-F1")#f[0] = uniprot ID some proteins have multiple models.  We are sticking to model 1
	#ab.write("\n %s"%f[0])
	gh = f[1].split("pdb")#gi[1] = atom type
	gi = gh[1].split("_")

	total_aromatic = 0
	total_neg = 0
	total_pos = 0
	templistA = [] 
	if gi[1] == 'CG':
		total_proteins +=1
		fi = open('%s'%d[1],'r')
		templistB = []
		for line in fi:
			result = re.split('\s+',line)
			templistB.append(result)#templistB contains lines in the individual sa file
		fi.close()
		for fj in templistB:
#fj[2] = res#, fj[3] = resname, fj[4] = atom fj[5] = surface area
#select only aromatic residues for further study
			#elif fj[3] in aromatics:
				#total_aromatic +=1
				for i in range(0, len(rrmlist)):
					if len(rrmlist[i])>=4:
						if rrmlist[i][0]==f[0]:
#identify the RRM domain range so that only residues in RRM domains will be analyzed in later steps (**make sure multiple RRM domains are included)
							ID = rrmlist[i][0] 	
							low = int(rrmlist[i][4])
							high = int(rrmlist[i][2])#keeping in mind that some proteins have several RRM domains, need a new output for each domain, not each protein
#use for loop to iterate through residue numbers in sa file and an if statement to isolate residues that are <=high and >= low
							#for p in templistB:
							#	if int(p[2])<=high and int(p[2])>=low:
							if int(fj[2])<=high and int(fj[2])>=low:
								if fj[3] in negatives:
									total_neg +=1
								if fj[3] in positives:
									total_pos +=1
#fj[2] = res#, fj[3] = resname, fj[4] = atom fj[5] = surface area
#select only aromatic residues for further study
								if fj[3] in aromatics:
									total_aromatic +=1
									templistA.append(fj)
									sum_sa = float(fj[5])#create a list of floats so that surface areas of atoms can be added together for each residue
									for filename in os.scandir(directory):
										#list of all files
										bb = (r.sub("", str(filename))) 
										cb = bb.split('>')
										db = cb[0].split("'")#d[1] = file name
										eb = db[1].split("AF-")
										fb = eb[1].split("-F")#f[0] = uniprot ID some proteins have multiple models.  We are sticking to model 1
										#fb = eb[1].split("-F1")#f[0] = uniprot ID
										ghb = fb[1].split("pdb")
										gib = ghb[1].split("_")#gi[1] = atom type
										if fb[0] == f[0]:
											templistC = []
											zh = open('%s'%db[1],'r')
											for line in zh:
												result = re.split('\s+',line)
												templistC.append(result)
											for i in templistC:
												if i[2]==fj[2]:
													sum_sa = sum_sa+float(i[5])
									ab.write(" %s"%str(f[0]))
									ID_ab.append(str(f[0]))
									ab.write(" %s"%str(pdbfile))
									pdb_ab.append(str(pdbfile))
									#ab.write(" %s"%this_domain)
									#domain_ID_ab.append(this_domain)
									ab.write(" %i"%int(fj[2]))#where ab is a list of sum_sas #not sure why it is writing multiple copies of each line but seems to work otherwise. Need to create an intermediate list somewhere
									resid_ab.append(int(fj[2]))
									ab.write(" %s"%fj[3])
									resn_ab.append(str(fj[3]))
									ab.write(" %s"%str(sum_sa))
									SA_ab.append(str(sum_sa))
									if fj[3]=='HIS':
										ab.write(" %f"%(float(sum_sa)/His_SA_poss))
										fraction_ab.append((float(sum_sa)/His_SA_poss))
										sa1 = float(sum_sa)/His_SA_poss
									elif fj[3]=='TRP':
										ab.write(" %f"%(float(sum_sa)/Trp_SA_poss))
										fraction_ab.append((float(sum_sa)/Trp_SA_poss))
										sa1 = float(sum_sa)/Trp_SA_poss
									elif fj[3]=='TYR':
										ab.write(" %f"%(float(sum_sa)/Tyr_SA_poss))
										fraction_ab.append((float(sum_sa)/Tyr_SA_poss))
										sa1 = float(sum_sa)/Tyr_SA_poss
									elif fj[3]=='PHE':
										ab.write(" %f"%(float(sum_sa)/Phe_SA_poss))
										fraction_ab.append((float(sum_sa)/Phe_SA_poss))
										sa1 = float(sum_sa)/Phe_SA_poss
									ab.write("\n")
									if float(sa1) >= thresh:
										SA_aromatic += 1


		ab.write('\n')
		simple_counts.write("\n protein: %s total negatives: %s total positives: %s total aromatics: %s SA_aromatic %s"%(f[0], total_neg, total_pos, total_aromatic, SA_aromatic))
		#s_counts +=("\n protein: %s total negatives: %s total positives: %s total aromatics: %s surface-exposed aromatics %s"%(f[0], total_neg, total_pos, total_aromatic, SA_aromatic))
		SA_aromatic = 0
ab.close()

totalxyz.write("total_#_of_proteins: %s \n"%total_proteins)	
####PART TWO: find charge clusters neighboring surface exposed aromatics
aromatics = []
neighbors = []
nonredundant_neighbors = ["dummy dummy","dummy dummy"]
per_pdb_totals = [] #one line = each total for each pdb file
for i in range(0, len(pdb_ab)): #iterates through list of pdb files
	neg = 0	
	pos = 0
	aromatic = 0
	SA_aromatic = 0
	pdb_doc = open("%s"%pdb_ab[i],"r")#opens a single pdb file
	pdb = []
	for j in pdb_doc:
		result = re.split('\s+',j)
		
		pdb.append(result)
	used_residues = []
	for k in range(0, len(pdb)): #iterates through a single pdb file
		if pdb[k][0]=="ATOM":
#if residue number is greater than 1000 for the alphafold structures, chain ID bumps into residue number, so splitting by space makes line off
#correct for this by rewriting line to split manually
			if pdb[k][4][0] in chainID:#should be true for all alphafold structures
				split_string = re.split(r'A|B|C|D|E|F|G|H|I|J|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z',pdb[k][4])
				if len(split_string[1])>0:#this should only be true if columns are bumping into each other
					try:
						test4.write("%s %s %s \n"%(split_string,len(split_string),ID_ab[i]))
						pdb[k] = (pdb[k][0], pdb[k][1], pdb[k][2], pdb[k][3], pdb[k][4][0], pdb[k][4][2], pdb[k][5], pdb[k][6], pdb[k][7], pdb[k][8], pdb[k][9], pdb[k][10], pdb[k][11]) 
					except IndexError:#if this problem coincides with the next, there will not be enough columns to split this way
						try:
							test4.write("%s %s %s \n"%(split_string,len(split_string),ID_ab[i]))
							pdb[k] = (pdb[k][0], pdb[k][1], pdb[k][2], pdb[k][3], pdb[k][4][0], pdb[k][4][2], pdb[k][5], pdb[k][6], pdb[k][7], pdb[k][8], pdb[k][9], pdb[k][10])
						except IndexError:
							test4.write("%s %s %s \n"%(split_string,len(split_string),ID_ab[i]))
							pdb[k] = (pdb[k][0], pdb[k][1], pdb[k][2], pdb[k][3], pdb[k][4][0], pdb[k][4][2], pdb[k][5], pdb[k][6], pdb[k][7], pdb[k][8], pdb[k][9]) 
					 
			else:#if there is no chainID here
				pdb[k] = (pdb[k][0], pdb[k][1], pdb[k][2], pdb[k][3], "dummy_chainID", pdb[k][4], pdb[k][5], pdb[k][6], pdb[k][7], pdb[k][8], pdb[k][9], pdb[k][10])
			#for some structures, the x/y/z coordinates overlap.  A better way would have been to use biopython to read in the pdb files, but this is how I corrected for this:
			try:
				float(pdb[k][6])
			except ValueError:
				split_string = re.split(r'-',pdb[k][6])
				test7.write("x: %s len%s"%(split_string,len(split_string)))
				if len(split_string)==3:
					a = 0-float(split_string[1])
					b = 0-float(split_string[2])
					pdb[k] = [pdb[k][0], pdb[k][1], pdb[k][2], pdb[k][3], pdb[k][4], pdb[k][5], a, b, pdb[k][7], pdb[k][8], pdb[k][9], pdb[k][10]]
					test7.write("%s %s \n"%(ID_ab[i],pdb[k]))

				elif len(split_string)==2:
					a = float(split_string[0])
					b = 0-float(split_string[1])
					pdb[k] = [pdb[k][0], pdb[k][1], pdb[k][2], pdb[k][3], pdb[k][4], pdb[k][5], a, b, pdb[k][7], pdb[k][8], pdb[k][9], pdb[k][10]]   
					test7.write("%s %s \n"%(ID_ab[i],pdb[k]))
				else:
					test7.write("fell through the cracks")
			try:
				float(pdb[k][7])
			except ValueError:
				split_string = re.split(r'-',pdb[k][7])
				test7.write("y: %s len%s \n"%(split_string,len(split_string)))
				if len(split_string)==3:
					a = 0-float(split_string[1])
					b = 0-float(split_string[2])
					pdb[k] = [pdb[k][0], pdb[k][1], pdb[k][2], pdb[k][3], pdb[k][4], pdb[k][5], pdb[k][6], a, b, pdb[k][8], pdb[k][9], pdb[k][10]]
					test7.write("%s %s \n"%(ID_ab[i],pdb[k]))
				elif len(split_string)==2:
					a = float(split_string[0])
					b = 0-float(split_string[1])
					pdb[k] = [pdb[k][0], pdb[k][1], pdb[k][2], pdb[k][3], pdb[k][4], pdb[k][5], pdb[k][6], a, b, pdb[k][8], pdb[k][9], pdb[k][10]]                    
					test7.write("%s %s \n"%(ID_ab[i],pdb[k]))
				else:
					test7.write("fell through the cracks")

			
			zqb.write("%s %s %s \n"%(ID_ab[i],resid_ab[i],pdb[k][5]))
			if int(resid_ab[i]) == int(pdb[k][5])and float(fraction_ab[i]) >= float(thresh): #find lines in the pdb file that correspond to the aromatic residue identified in the resid list (will only take one residue at a time until i becomes i+1)  
				if pdb[k][2] in query_atom_type:
				#SA_aromatic = SA_aromatic+1
				#for x in query_atom_type:
				#	if bool(re.search("%s"%x,pdb[k][2])) == True:	#not working, supposed to only show atoms whose interactions we want to observe.  May need to remove wildcard and type all values to make work
						x_q = float(pdb [k][6])
						y_q = float(pdb [k][7])
						z_q = float(pdb [k][8])
						
						#now that we have the 'query' coordinates, we need to iterate through the pdb file 'subject' coordinates again to identify hits
						for h in range (0, len(pdb)):
							if pdb[h][0]=="ATOM":
								if str("%s %s"%(pdb[h][3],pdb[h][5])) in neighbors:
									pass
								else:
									if pdb[h][4][0] in chainID:
										split_string = re.split(r'A|B|C|D|E|F|G|H|I|J|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z',pdb[h][4])
										if len(split_string[1])>0:
											test4.write("%s %s %s \n"%(split_string,len(split_string),ID_ab[i]))
											pdb[h] = (pdb[h][0], pdb[h][1], pdb[h][2], pdb[h][3], pdb[h][4][0], pdb[h][4][2], pdb[h][5], pdb[h][6], pdb[h][7], pdb[h][8], pdb[h][9], pdb[h][10], pdb[h][11]) 
									else:
										pdb[h] = (pdb[h][0], pdb[h][1], pdb[h][2], pdb[h][3], "dummy_chainID", pdb[h][4], pdb[h][5], pdb[h][6], pdb[h][7], pdb[h][8], pdb[h][9], pdb[h][10])
									#for some structures, the x/y/z coordinates overlap and cannot be split by space.  This seems to only happen when either x or y is a very negative number.  Therefore to correct for this:
									try:
										float(pdb[h][6])
									except ValueError:
										split_string = re.split(r'-',pdb[h][6])
										test7.write("x: %s len%s"%(split_string,len(split_string)))
										if len(split_string)==3:
											a = 0-float(split_string[1])
											b = 0-float(split_string[2])
											pdb[h] = [pdb[h][0], pdb[h][1], pdb[h][2], pdb[h][3], pdb[h][4], pdb[h][5], a, b, pdb[h][7], pdb[h][8], pdb[h][9], pdb[h][10]]
											test7.write("%s %s \n"%(ID_ab[i],pdb[h]))
										elif len(split_string)==2:
											a = float(split_string[0])
											b = 0-float(split_string[1])
											pdb[h] = [pdb[h][0], pdb[h][1], pdb[h][2], pdb[h][3], pdb[h][4], pdb[h][5], a, b, pdb[h][7], pdb[h][8], pdb[h][9], pdb[h][10]]   
											test7.write("%s %s \n"%(ID_ab[i],pdb[h]))
										else:
											test7.write("fell through the cracks")
									try:
										float(pdb[h][7])
									except ValueError:
										split_string = re.split(r'-',pdb[h][7])
										test7.write("y: %s len%s \n"%(split_string,len(split_string)))
										if len(split_string)==3:
											a = 0-float(split_string[1])
											b = 0-float(split_string[2])
											pdb[h] = [pdb[h][0], pdb[h][1], pdb[h][2], pdb[h][3], pdb[h][4], pdb[h][5], pdb[h][6], a, b, pdb[h][8], pdb[h][9], pdb[h][10]]
											test7.write("%s %s \n"%(ID_ab[i],pdb[h]))
										elif len(split_string)==2:
											a = float(split_string[0])
											b = 0-float(split_string[1])
											pdb[h] = [pdb[h][0], pdb[h][1], pdb[h][2], pdb[h][3], pdb[h][4], pdb[h][5], pdb[h][6], a, b, pdb[h][8], pdb[h][9], pdb[h][10]]                    
											test7.write("%s %s \n"%(ID_ab[i],pdb[h]))
										else:
											test7.write("fell through the cracks")

									if pdb[h][2] in subject_atom_type:	
											x_s = float(pdb[h][6])
											y_s = float(pdb[h][7])
											z_s = float(pdb[h][8])			
											distance = math.sqrt( (x_q-x_s)**2 + (y_q-y_s)**2 + (z_q-z_s)**2 )
											if distance <= dist:
												#qh.write("%s"%pdb[h][3])
												#qh.write("_%s"%pdb[h][5])
												#qh.write("_%s "%pdb[h][2])
												neighbors.append(str("%s %s"%(pdb[h][3],pdb[h][5])))
						nonredundant_neighbors = []
						for qm in neighbors:
							if str(qm) in nonredundant_neighbors:
								pass
							else:
								nonredundant_neighbors.append(str(qm))
								z.write("%s"%str(qm))
								z.write("\n")


						
						try:
							if pdb[k][5] == pdb[k+1][5]:#we don't want to start counting if the next atom is part of the same residue because the lists are not yet complete
								pass
							else:#we don't want to start counting if the next atom is part of the same residue because the lists are not yet complete
								negative_neighbors = 0
								positive_neighbors = 0
								other_neighbors = 0
								for line in nonredundant_neighbors:
									e = line.split(" ")
									if str(e[0]) in negatives: 
										negative_neighbors += 1
									elif str(e[0]) in positives: 
										positive_neighbors += 1
									else:
										other_neighbors += 1
								qh.write(" %s %s %s %s negative neighbors: %i positive neighbors: %i other neighbors: %i "%(ID_ab[i], resid_ab[i], resn_ab[i], pdb[k][2], int(negative_neighbors),int(positive_neighbors),int(other_neighbors)))
								if ID_ab[i] in phasep_ident:
									qh.write("phasep: y ")
									a = 'phasep:y'
								else:
									qh.write("phasep: n ")
									a = 'phasep:n'
								if ID_ab[i] in vitro_prot:
									qh.write("vitro_prot: y ")
									b = 'vitro_prot:y'
								else:
									qh.write("vitro_prot: n ")
									b = 'vitro_prot:n'
								if ID_ab[i] in combined:
									qh.write("combined: y ")
									c = 'combined:y'
								else:
									qh.write("combined: n ")
									c = 'combined:n'
								if ID_ab[i] in DrLLPs:
									qh.write("DrLLPs: y ")
									d = 'DrLLPs:y'
								else:
									qh.write("DrLLPs: n ")
									d = 'DrLLPs:n'
								aromatics.append("%s %s %s %s %s %s %s %s %s %s %s"%(ID_ab[i], resid_ab[i], resn_ab[i], pdb[k][2], negative_neighbors,positive_neighbors,other_neighbors,a,b,c,d))
								for line in neighbors:
									qh.write(" %s,"%line)
								qh.write('\n')
								neighbors = []
								#neighbors.append("new neighbors %s_%s_%s: "%(pdb[k][5],pdb[k][3],pdb[k][2]))
						except IndexError:
								if k == (len(pdb)-1): #if we are looking at the last residue in the pdb file, the above will throw an error
									negative_neighbors = 0
									positive_neighbors = 0
									other_neighbors = 0
									for line in nonredundant_neighbors:
										e = line.split(" ")
										if str(e[0]) in negatives: 
											negative_neighbors += 1
										elif str(e[0]) in positives: 
											positive_neighbors += 1
										else:
											other_neighbors += 1
									qh.write(" %s %s %s %s negative neighbors: %i positive neighbors: %i other neighbors: %i "%(ID_ab[i], resid_ab[i], resn_ab[i], pdb[k][2], int(negative_neighbors),int(positive_neighbors),int(other_neighbors)))
									if ID_ab[i] in phasep_ident:
										qh.write("phasep: y ")
										a = 'phasep:y'
									else:
										qh.write("phasep: n ")
										a = 'phasep:n'
									if ID_ab[i] in vitro_prot:
										qh.write("vitro_prot: y ")
										b = 'vitro_prot:y'
									else:
										qh.write("vitro_prot: n ")
										b = 'vitro_prot:n'
									if ID_ab[i] in combined:
										qh.write("combined: y ")
										c = 'combined:y'
									else:
										qh.write("combined: n ")
										c = 'combined:n'
									if ID_ab[i] in DrLLPs:
										qh.write("DrLLPs: y ")
										d = 'DrLLPs:y'
									else:
										qh.write("DrLLPs: n ")
										d = 'DrLLPs:n'
									aromatics.append("%s %s %s %s %s %s %s %s %s %s %s"%(ID_ab[i], resid_ab[i], resn_ab[i], pdb[k][2], negative_neighbors,positive_neighbors,other_neighbors,a,b,c,d))
									for line in neighbors:
										qh.write(" %s,"%line)
									qh.write('\n')
									neighbors = []
					
								else:
									qh.write("unknown index error: index: %i len(pdb): %i pdb[k]: %s ID_ab[i]: %s pdb[k][5]: %s pdb[k+1][5]: %s resid_ab[i]: %s \n"%(i,len(pdb),pdb[k],ID_ab[i], pdb[k][5], pdb[k+1][5],  resid_ab[i]))

totalxyz.close()		
