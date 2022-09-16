import re as re
numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'  #to allow residue to be read in as a float
rx = re.compile(numeric_const_pattern, re.VERBOSE)

first = open("domains.txt","r")
zhang = open("Dr_Zhang_phase_separate.txt","r")
zhanglist = [a for a in zhang]
phasep = open("phasepdb2_MLOhtp.txt","r")
combo_ = open("combo.txt","r")
llpsprotein = open("llps_proteinc.txt","r")
llps_RNA = open("llps_protein+RNAc.txt","r")
DrLLPS = open("drllpsc.txt","r")

list = []
for line in first:
    list.append(line.split(' '))
first.close()
combo = [re.split('\s',a)[0] for a in combo_]
plist = []
ident = []
location = []
for line in phasep:
    if len(line)>2:
        try:
            a = re.split('\s+',line)
            location.append(a)
            ident.append(a[1])
        except IndexError:
            print("failed to load phasep line: %s"%str(line))

vitro_prot = []
vitro_RNP = []
a = []
b = []
c = []
for line in llpsprotein:
    a.append(line.split('\n'))
for line in llps_RNA:
    b.append(line.split('\n'))

for i in range(len(a)):
    vitro_prot.append(a[i][0])
for j in range(len(b)):
    vitro_RNP.append(b[j][0])
llpsprotein.close()
llps_RNA.close()


DrLLPs = []
comma_split_DrLLPs = []
for line in DrLLPS:
	result = re.split('\s+',line)
	c.append(result)
for k in range(len(c)):
	DrLLPs.append(c[k][1])
DrLLPS.close()

ID_list = []
RRM_list = []
RRM_list_out = open("RRM_list_uniprot_3.txt","w")
for i in range(0, len(list)):
    if len(list[i])>3:
        if 'RRM' in list[i][3]:
            ID = list[i][0] 
            low = int(list[i][1])
            high = int(list[i][2])
            seq = list[i][4]
            RRMseq = seq[low:high]
            nRK = RRMseq.count('R')+RRMseq.count('K')
            nDE = RRMseq.count('D')+RRMseq.count('E')
            #nFYHW = RRMseq.count('F')+RRMseq.count('Y')+RRMseq.count('H')+RRMseq.count('W')
            nFYW = RRMseq.count('F')+RRMseq.count('Y')+RRMseq.count('W')
            nH = RRMseq.count('H')
            RRM_list.append(RRMseq)
            ID_list.append(ID)
            RRM_list_out.write(ID)
            RRM_list_out.write(" ")
            RRM_list_out.write("High %s "%high)
            RRM_list_out.write("Low %s "%low)
            RRM_list_out.write(RRMseq)
            RRM_list_out.write(" RK=%f" %nRK)
            RRM_list_out.write(" DE=%f" %nDE)
            RRM_list_out.write(" FYW=%f" %nFYW)
            RRM_list_out.write(" H=%f" %nH)
            
            if ID in ident:
                RRM_list_out.write(" phasepdb= y")
            else:
                RRM_list_out.write(" phasepdb= n")

            

            if ID in vitro_prot:
                RRM_list_out.write(" llbsdb_prot= y")
            else:
                RRM_list_out.write(" llbsdb_prot= n")
                
                
            if ID in combo:
                RRM_list_out.write(" combo= y")
            else:
                RRM_list_out.write(" combo= n")

            if ID in DrLLPs:
                RRM_list_out.write(" DRLLPS= y")
            else:
                RRM_list_out.write(" DRLLPS= n")
            


                
            RRM_list_out.write(" location_per_phasepdb: ")
            for m in location:
                if m[1]==ID:
                    RRM_list_out.write("%s,"%str(m[2]))
            if ID in ident:
                pass
            else:
                RRM_list_out.write("-")

            RRM_list_out.write(" location_per_DrLLPS: ")
            for x in c:
                if ID == x[1]:
                    RRM_list_out.write(x[5])
            if ID in DrLLPs:
                pass
            else:
                RRM_list_out.write("-")
            if ID in zhanglist:
                RRM_list_out.write(" zhang= y")
            else:
                RRM_list_out.write(" zhang= n")
            if ID in vitro_RNP:
                RRM_list_out.write(" llbsdb_RNA+prot= y")
            else:
                RRM_list_out.write(" llbsdb_RNA+prot= n")
                
                
            RRM_list_out.write('\n')
RRM_list_out.close()
        
