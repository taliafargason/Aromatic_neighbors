import re
import statistics
test = open("raw_data_for_violin.txt","w")
test2 = open("test_2_2.txt","w")
test4 = open("test4_2.txt","w")
by_residue_neighbors = open("by_residue_neighbors_6_0.6.txt","r")

pdbs = []
simple = open("simple_counts_any_residue_in_an_RRM_6_0.6.txt","r")
for line in simple:
    result = re.split('\s+',line)
    pdbs.append(result[2])

by_residue = []
IDs_by_residue = []
for line in by_residue_neighbors:
        a = re.split('\s+',line)
        by_residue.append(a)
        IDs_by_residue.append(a[1])

by_residue_neighbors.close()
mer2RS = open("report_SR.txt","r")
RS2 = [re.split('\s+',line)for line in mer2RS]
RS2ident = [line[0] for line in RS2]
mer4RS = open("report_SRSR.txt","r")
RS4 = [re.split('\s+',line)for line in mer4RS]
RS4ident = [line[0] for line in RS4] #if len(line)==4]
mer6RS = open("report_SRSRSR.txt","r")
RS6 = [re.split('\s',line)for line in mer6RS]
RS6ident = [line[0] for line in RS6]# if len(line)==4]
mer8RS = open("report_SRSRSRSR.txt","r")
RS8 = [re.split('\s',line)for line in mer8RS]
RS8ident = [line[0] for line in RS8]# if len(line)==4]
for line in RS4:
        test.write(str(line))
        a = line[0] in RS4ident
        test.write(line[0])
        test.write(str(a))
RRM_list = []
RRM = open("RRM_list_uniprot_3.txt","r")
for line in RRM:
        a = re.split('\s+',line)
        if len(a)>2:
            RRM_list.append(a)
RRM.close()

sa_aromatics_counts = []
aromatic_neighbor_counts = []
determine_max_aromatic = []
determine_max_neg_neighbor = []
determine_max_pos_neighbor = []
out1 = open("Aromatics_meeting_criteria_by_domain_out_1_THRESHOLD=1neighbors.txt","w")

def tally_up( line ):
        if line[0] in pdbs:
            sa_aromatics = 0
            aromatics_w_neg_neighbors = 0
            aromatics_w_pos_neighbors = 0
            if line[0] in IDs_by_residue:
                for i in by_residue:
                        if i[1] == line[0] and float(i[2]) <= float(line[2]) and float(i[2]) >= float(line[4]):
                                sa_aromatics +=1
                                if int(i[7])>0:
                                        aromatics_w_neg_neighbors +=1
                                if int(i[10])>0:
                                        aromatics_w_pos_neighbors +=1
                sa_aromatics_counts.append("%s %s %s %s %s %s %s %s"%(line[0],line[2],line[4],sa_aromatics,line[11],line[13],line[15],line[17]))
                test.write("\n%s %s %s %s %s %s %s %s"%(line[0],line[2],line[4],sa_aromatics,line[11],line[13],line[15],line[17]))                
                determine_max_aromatic.append(sa_aromatics)
                determine_max_neg_neighbor.append(aromatics_w_neg_neighbors)
                determine_max_pos_neighbor.append(aromatics_w_pos_neighbors)
                aromatic_neighbor_counts.append("%s %s %s %s %s %s %s %s %s "%(line[0],line[2],line[4],aromatics_w_neg_neighbors,aromatics_w_pos_neighbors,line[11],line[13],line[15],line[17]))
                test.write("neighbor %s %s %s %s %s %s %s %s %s \n"%(line[0],line[2],line[4],aromatics_w_neg_neighbors,aromatics_w_pos_neighbors,line[11],line[13],line[15],line[17]))
       
            else:
                sa_aromatics_counts.append("%s %s %s %s %s %s %s %s"%(line[0],line[2],line[4],int('0'),line[11],line[13],line[15],line[17]))
                test.write("\n%s %s %s %s %s %s %s %s"%(line[0],line[2],line[4],int('0'),line[11],line[13],line[15],line[17]))

def mean_median_mode(q):
        vitro_prot = [0,0]
        nada = [0,0]
        vivo_vitro = [0,0]
        for a in q:
                vitro_prot = []
                nada = []
                vivo_vitro = []
                for a in q:
                    j = a.split(" ")
                    if j[5]=='y':
                        vitro_prot.append(int(j[3]))
                    if j[4]=='y' or j[5]=='y' or j[6]=='y' or j[7]=='y':
                        vivo_vitro.append(int(j[3]))
                    elif j[4]=='n' and j[5]=='n' and j[6]=='n' and j[7]=='n':
                        nada.append(int(j[3]))
                    else:
                        print('error in mean_median_mode: phase separation categories misaligned!',file=test4)
        try:
            if nada ==[] and vitro_prot==[]:
                out1.write(' aromatics   any phase separation: mean:%s stdev:%s median:%s mode:%s n:%s\n'%(statistics.mean(vivo_vitro),statistics.stdev(vivo_vitro),statistics.median(vivo_vitro),statistics.mode(vivo_vitro),len(vivo_vitro)))
            elif nada == []:
                out1.write(' aromatics   in vitro verified: mean:%s stdev:%s median:%s mode:%s n:%s\n'%(statistics.mean(vitro_prot),statistics.stdev(vitro_prot),statistics.median(vitro_prot),statistics.mode(vitro_prot),len(vitro_prot)))
                out1.write(' aromatics   any phase separation: mean:%s stdev:%s median:%s mode:%s n:%s\n'%(statistics.mean(vivo_vitro),statistics.stdev(vivo_vitro),statistics.median(vivo_vitro),statistics.mode(vivo_vitro),len(vivo_vitro)))
            else:
                out1.write(' aromatics   nollps: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(nada),statistics.stdev(nada),statistics.median(nada),statistics.mode(nada)))
                out1.write(' aromatics   in vitro verified: mean:%s stdev:%s median:%s mode:%s n:%s\n'%(statistics.mean(vitro_prot),statistics.stdev(vitro_prot),statistics.median(vitro_prot),statistics.mode(vitro_prot),len(vitro_prot)))
                out1.write(' aromatics   any phase separation: mean:%s stdev:%s median:%s mode:%s n:%s\n'%(statistics.mean(vivo_vitro),statistics.stdev(vivo_vitro),statistics.median(vivo_vitro),statistics.mode(vivo_vitro),len(vivo_vitro)))
        except:
            out1.write('aromatics statistics error \n')
def neg_mean_median_mode(q):
    vitro_prot = [0,0]
    nada = [0,0]
    vivo_vitro = [0,0]
    for a in q:
            vitro_prot = []
            nada = []
            vivo_vitro = []
            for a in q:
                j = a.split(" ")
                if j[6]=='y':
                    vitro_prot.append(int(j[3]))
                if j[5]=='y' or j[6]=='y' or j[7]=='y' or j[8]=='y':
                    vivo_vitro.append(int(j[3]))
                elif j[5]=='n' and j[6]=='n' and j[7]=='n' and j[8]=='n':
                    nada.append(int(j[3]))
                else:
                    print('error in mean_median_mode: phase separation categories misaligned!',file=test4)
    try:
        if nada==[] and vitro_prot==[]:
            out1.write(' neg_neighbor   any phase separation: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(vivo_vitro),statistics.stdev(vivo_vitro),statistics.median(vivo_vitro),statistics.mode(vivo_vitro))) 
        elif nada ==[]:
            #out1.write('nollps: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(nada),statistics.stdev(nada),statistics.median(nada),statistics.mode(nada)))
            out1.write(' neg_neighbor   in vitro verified: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(vitro_prot),statistics.stdev(vitro_prot),statistics.median(vitro_prot),statistics.mode(vitro_prot)))
            out1.write(' neg_neighbor   any phase separation: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(vivo_vitro),statistics.stdev(vivo_vitro),statistics.median(vivo_vitro),statistics.mode(vivo_vitro)))
        else:
            out1.write(' neg_neighbor   nollps: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(nada),statistics.stdev(nada),statistics.median(nada),statistics.mode(nada)))
            out1.write(' neg_neighbor   in vitro verified: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(vitro_prot),statistics.stdev(vitro_prot),statistics.median(vitro_prot),statistics.mode(vitro_prot)))
            out1.write(' neg_neighbor   any phase separation: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(vivo_vitro),statistics.stdev(vivo_vitro),statistics.median(vivo_vitro),statistics.mode(vivo_vitro)))
    except:
        out1.write('statistics error \n')
def pos_mean_median_mode(q):
    vitro_prot = [0,0]
    nada = [0,0]
    vivo_vitro = [0,0]
    for a in q:
        #if q=='sa_aromatics_counts':
        #make a list for each category (phase separated vs not, find mean,stdev,median,mode)
            vitro_prot = []
            nada = []
            vivo_vitro = []
            for a in q:
                j = a.split(" ")
                if j[6]=='y':
                    vitro_prot.append(int(j[4]))
                if j[5]=='y' or j[6]=='y' or j[7]=='y' or j[8]=='y':
                    vivo_vitro.append(int(j[4]))
                elif j[5]=='n' and j[6]=='n' and j[7]=='n' and j[8]=='n':
                    nada.append(int(j[4]))
                else:
                    print('error in mean_median_mode: phase separation categories misaligned!',file=test4)
    try:
        if nada==[] and vitro_prot==[]:
            out1.write(' pos_neighbor   any phase separation: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(vivo_vitro),statistics.stdev(vivo_vitro),statistics.median(vivo_vitro),statistics.mode(vivo_vitro)))  
        elif nada==[]:
            #out1.write('nollps: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(nada),statistics.stdev(nada),statistics.median(nada),statistics.mode(nada)))
            out1.write(' pos_neighbor   in vitro verified: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(vitro_prot),statistics.stdev(vitro_prot),statistics.median(vitro_prot),statistics.mode(vitro_prot)))
            out1.write(' pos_neighbor   any phase separation: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(vivo_vitro),statistics.stdev(vivo_vitro),statistics.median(vivo_vitro),statistics.mode(vivo_vitro)))  
        else:
            out1.write(' pos_neighbor   nollps: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(nada),statistics.stdev(nada),statistics.median(nada),statistics.mode(nada)))
            out1.write(' pos_neighbor   in vitro verified: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(vitro_prot),statistics.stdev(vitro_prot),statistics.median(vitro_prot),statistics.mode(vitro_prot)))
            out1.write(' pos_neighbor   any phase separation: mean:%s stdev:%s median:%s mode:%s \n'%(statistics.mean(vivo_vitro),statistics.stdev(vivo_vitro),statistics.median(vivo_vitro),statistics.mode(vivo_vitro)))
    except:
        out1.write('statistics error \n')
                   
def tally_aromatics(i):
        no_llps = 0
        any_condensate = 0
        in_vitro_protein = 0
        for a in sa_aromatics_counts:
                j = a.split(" ")
                #print(j[3])
                if int(j[3])== i:
                        if j[5]=='y':
                                in_vitro_protein+=1
                        if j[4]=='y' or j[5]=='y' or j[6]=='y' or j[7]=='y':
                                any_condensate+=1
                        else:
                                no_llps+=1
        sa_aromatic_histogram_no_llps.append("%s %s \n"%(i,no_llps))
        out1.write("%s %s"%(i,no_llps))
        sa_aromatic_histogram_any_condensate.append("%s %s \n"%(i,any_condensate))
        out1.write("    %s %s"%(i,any_condensate))
        sa_aromatic_histogram_in_vitro_protein.append("%s %s \n"%(i,in_vitro_protein))
        out1.write("    %s %s \n"%(i,in_vitro_protein))
def tally_neg_neighbors(i):
        no_llps = 0
        any_condensate = 0
        in_vitro_protein = 0
        for a in aromatic_neighbor_counts:
                j = a.split(" ")

                if int(j[3])== i:
                        if j[6]=='y':
                                in_vitro_protein+=1
                        if j[5]=='y' or j[6]=='y' or j[7]=='y' or j[8]=='y':
                                any_condensate+=1
                        else:
                                no_llps+=1
        sa_aromatic_neg_no_llps.append("%s %s \n"%(i,no_llps))
        out1.write("%s %s"%(i,no_llps))
        sa_aromatic_neg_any_condensate.append("%s %s \n"%(i,any_condensate))
        out1.write("    %s %s"%(i,any_condensate))
        sa_aromatic_neg_in_vitro_protein.append("%s %s \n"%(i,in_vitro_protein))
        out1.write("    %s %s \n"%(i,in_vitro_protein))
def  tally_pos_neighbors(i):
        no_llps = 0
        any_condensate = 0
        in_vitro_protein = 0
        for a in aromatic_neighbor_counts:
                j = a.split(" ")

                if int(j[4])== i:
                        if j[6]=='y':
                                in_vitro_protein+=1
                        if j[5]=='y' or j[6]=='y' or j[7]=='y' or j[8]=='y':
                                any_condensate+=1
                        else:
                                no_llps+=1
        sa_aromatic_pos_no_llps.append("%s %s \n"%(i,no_llps))
        out1.write("%s %s"%(i,no_llps))
        sa_aromatic_pos_any_condensate.append("%s %s \n"%(i,any_condensate))
        out1.write("         %s %s"%(i,any_condensate))
        sa_aromatic_pos_in_vitro_protein.append("%s %s \n"%(i,in_vitro_protein))
        out1.write("                %s %s \n"%(i,in_vitro_protein))        

for line in RRM_list:

        #when calling the function later, use an if function to choose which lines are used
        tally_up(line)
sa_aromatic_histogram_no_llps = []
sa_aromatic_histogram_any_condensate = []
sa_aromatic_histogram_in_vitro_protein = []

sa_aromatic_neg_no_llps = []
sa_aromatic_neg_any_condensate = []
sa_aromatic_neg_in_vitro_protein = []

sa_aromatic_pos_no_llps = []
sa_aromatic_pos_any_condensate = []
sa_aromatic_pos_in_vitro_protein = []
out1.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
#to call mean_median_mode: where i is either sa_aromatics_counts or aromatic_neighbor_counts, j ==either 'neg' or 'pos'
mean_median_mode(sa_aromatics_counts)
for i in range(0,(max(determine_max_aromatic)+1)):
    tally_aromatics(i)

out1.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
#to call mean_median_mode: where i is either sa_aromatics_counts or aromatic_neighbor_counts, j ==either 'neg' or 'pos'
#mean_median_mode(q = aromatic_neighbor_counts,k ='neg')
neg_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,(max(determine_max_neg_neighbor)+1)):
        tally_neg_neighbors(i)

out1.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
pos_mean_median_mode(aromatic_neighbor_counts)
#to call mean_median_mode: where i is either sa_aromatics_counts or aromatic_neighbor_counts, j ==either 'neg' or 'pos'
#mean_median_mode(q = aromatic_neighbor_counts,k ='pos')
for i in range(0,(max(determine_max_pos_neighbor)+1)):
        tally_pos_neighbors(i)

by_residue_SA = open('by_residue_surface_area_6_0.6.txt','r')
SA = []
for line in by_residue_SA:
        result = re.split('\s+',line)
        if len(result)>7:
                SA.append(result)
                test2.write(str(result))
                test2.write('\n')
test2.close()
#1=ID 2=PDB 3=domain# 4=res# 5=res_name 6=surface_exposure 7=fraction_exposed
PHE_llps = []
PHE_none = []
PHE_vitro = []

TYR_llps = []
TYR_none = []
TYR_vitro = []

HIS_llps = []
HIS_none = []
HIS_vitro = []

TRP_llps = []
TRP_none = []
TRP_vitro = []
SA_out = open('SA_out_totals.txt','w')

for i in SA:
        for j in RRM_list:
                if i[1]==j[0]:
                        if i[5]=='PHE':
                                if j[11]=='y' or j[13] =='y' or j[15]=='y' or j[17]=='y':
                                        PHE_llps.append(i[7])
                                        if j[13]=='y':
                                                PHE_vitro.append(i[7])
                                        else:
                                                PHE_none.append(i[7])
                        elif i[5]=='TYR':
                                if j[11]=='y' or j[13] =='y' or j[15]=='y' or j[17]=='y':
                                        TYR_llps.append(i[7])
                                        if j[13]=='y':
                                                TYR_vitro.append(i[7])
                                        else:
                                                TYR_none.append(i[7])
                        elif i[5]=='HIS':
                                if j[11]=='y' or j[13] =='y' or j[15]=='y' or j[17]=='y':
                                        HIS_llps.append(i[7])
                                        if j[13]=='y':
                                                HIS_vitro.append(i[7])
                                        else:
                                                HIS_none.append(i[7])
                        elif i[5]=='TRP':
                                if j[11]=='y' or j[13] =='y' or j[15]=='y' or j[17]=='y':
                                        TRP_llps.append(i[7])
                                        if j[13]=='y':
                                                TRP_vitro.append(i[7])
                                        else:
                                                TRP_none.append(i[7])

SA_out.write("PHE TYR HIS TRP\n none cond in-vit \n")
for i in range(0,len(SA)):
        
        try:
             SA_out.write(" %f3"%float(PHE_none[i]))
        except IndexError:
             SA_out.write(" -")
        try:
             SA_out.write(" %f3"%float(PHE_llps[i]))
        except IndexError:
             SA_out.write(" -")
        try:
             SA_out.write(" %f3"%float(PHE_vitro[i]))
        except IndexError:
             SA_out.write(" -")

        try:
             SA_out.write(" %f3"%float(TYR_none[i]))
        except IndexError:
             SA_out.write(" -")
        try:
             SA_out.write(" %f3"%float(TYR_llps[i]))
        except IndexError:
             SA_out.write(" -")
        try:
             SA_out.write(" %f3"%float(TYR_vitro[i]))
        except IndexError:
             SA_out.write(" -")
        try:
             SA_out.write(" %f3"%float(HIS_none[i]))
        except IndexError:
             SA_out.write(" -")
        try:
             SA_out.write(" %f3"%float(HIS_llps[i]))
        except IndexError:
             SA_out.write(" -")
        try:
             SA_out.write(" %f3"%float(HIS_vitro[i]))
        except IndexError:
             SA_out.write(" -")
        try:
             SA_out.write(" %f3"%float(TRP_none[i]))
        except IndexError:
             SA_out.write(" -")
        try:
             SA_out.write(" %f3"%float(TRP_llps[i]))
        except IndexError:
             SA_out.write(" -")
        try:
             SA_out.write(" %f3"%float(TRP_vitro[i]))
        except IndexError:
             SA_out.write(" -")
        SA_out.write('\n')
SA_out.close()

test2.close()

test3 = open('test3_2.txt','w')
#make a list of types of condensates studied remembering in RRMlist, phasepdblocation = 19, DrLLPSDBlocation = 21, then iterate through to add an additional list mimicking the one above for each type of condensate
condensate_types = []
for line in RRM_list:
        a = line[19].split(',')
        for b in a:
                if b in condensate_types:
                        pass
                else:
                     condensate_types.append(b)
                     test3.write(b)
                     test3.write('\n')
        c = line[21].split(',')
        for d in c:
                if d in condensate_types:
                        pass
                else:
                        condensate_types.append(d)
                        test3.write(d)
                        test3.write('\n')
test3.close()
test4 = open('test4_c.txt','w')

for p in condensate_types:
        if len(p)>2:
                out1.write('%s counts \n'%p)
                test.write('%s counts \n'%p)
                sa_aromatics_counts = []
                aromatic_neighbor_counts = []
                determine_max_aromatic = [0]
                determine_max_neg_neighbor = [0]
                determine_max_pos_neighbor = [0]
                for line in RRM_list:
                        a = line[19].split(',')
                        c = line[21].split(',')
                        if p in a or p in c:
       #starting with line 0, column0 = ID, column2 = high, column4 = low, 12 = phasepdb, 14 = vitro_prot, 16 = vitro_RNP, 18 = DrLLPs, phasepdblocation = 19, DrLLPSDBlocation = 21
                                tally_up(line)
                sa_aromatic_histogram_no_llps = []
                sa_aromatic_histogram_any_condensate = []
                sa_aromatic_histogram_in_vitro_protein = []

                sa_aromatic_neg_no_llps = []
                sa_aromatic_neg_any_condensate = []
                sa_aromatic_neg_in_vitro_protein = []

                sa_aromatic_pos_no_llps = []
                sa_aromatic_pos_any_condensate = []
                sa_aromatic_pos_in_vitro_protein = []
                out1.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
                mean_median_mode(sa_aromatics_counts)
                for i in range(0,max(determine_max_aromatic)+1):
                        tally_aromatics(i)

                out1.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
                neg_mean_median_mode(aromatic_neighbor_counts)
                for i in range(0,max(determine_max_neg_neighbor)+1):
                        tally_neg_neighbors(i)

                out1.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
                pos_mean_median_mode(aromatic_neighbor_counts)
                for i in range(0,max(determine_max_pos_neighbor)+1):
                        tally_pos_neighbors(i)

#mer4RS = open("4merRS_in_condensates.txt","r")
#RS4 = [re.split('\s+',line)for line in mer4RS]
#RS4ident = [line[0]for line in RS4 if len(line)==4]
#mer6RS = open("6merRS_in_condensates.txt","r")
#RS6 = [re.split('\s',line)for line in mer6RS]
#RS6ident = [line[0]for line in RS6 if len(line)==4]
#mer8RS = open("8merRS_in_condensates.txt","r")
#RS8 = [re.split('\s',line)for line in mer8RS]
#RS8ident = [line[0] for line in RS8 if len(line)==4]
out1.write("\n RS2 counts \n")
test.write("\n RS2 counts \n")
sa_aromatics_counts = []
aromatic_neighbor_counts = []
determine_max_aromatic = []
determine_max_neg_neighbor = []
determine_max_pos_neighbor = []
for line in RRM_list:
        if line[0] in RS2ident:
        #when calling the function later, use an if function to choose which lines are used
                tally_up(line)
sa_aromatic_histogram_no_llps = []
sa_aromatic_histogram_any_condensate = []
sa_aromatic_histogram_in_vitro_protein = []

sa_aromatic_neg_no_llps = []
sa_aromatic_neg_any_condensate = []
sa_aromatic_neg_in_vitro_protein = []

sa_aromatic_pos_no_llps = []
sa_aromatic_pos_any_condensate = []
sa_aromatic_pos_in_vitro_protein = []
out1.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
test.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
mean_median_mode(sa_aromatics_counts)
for i in range(0,max(determine_max_aromatic)+1):
        tally_aromatics(i)

out1.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
test.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
neg_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_neg_neighbor)+1):
        tally_neg_neighbors(i)

out1.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
test.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
pos_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_pos_neighbor)+1):
        tally_pos_neighbors(i)

out1.write("\n RS4 counts \n")
test.write("\n RS4 counts \n")
sa_aromatics_counts = []
aromatic_neighbor_counts = []
determine_max_aromatic = []
determine_max_neg_neighbor = []
determine_max_pos_neighbor = []
for line in RRM_list:
        if line[0] in RS4ident:
        #when calling the function later, use an if function to choose which lines are used
                tally_up(line)
sa_aromatic_histogram_no_llps = []
sa_aromatic_histogram_any_condensate = []
sa_aromatic_histogram_in_vitro_protein = []

sa_aromatic_neg_no_llps = []
sa_aromatic_neg_any_condensate = []
sa_aromatic_neg_in_vitro_protein = []

sa_aromatic_pos_no_llps = []
sa_aromatic_pos_any_condensate = []
sa_aromatic_pos_in_vitro_protein = []
out1.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
test.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
mean_median_mode(sa_aromatics_counts)
for i in range(0,max(determine_max_aromatic)+1):
        tally_aromatics(i)

out1.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
test.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
neg_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_neg_neighbor)+1):
        tally_neg_neighbors(i)

out1.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
test.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
pos_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_pos_neighbor)+1):
        tally_pos_neighbors(i)


out1.write("\n RS6 counts \n")
test.write("\n RS6 counts \n")
sa_aromatics_counts = []
aromatic_neighbor_counts = []
determine_max_aromatic = []
determine_max_neg_neighbor = []
determine_max_pos_neighbor = []
for line in RRM_list:
        if line[0] in RS6ident:
        #when calling the function later, use an if function to choose which lines are used
                tally_up(line)
sa_aromatic_histogram_no_llps = []
sa_aromatic_histogram_any_condensate = []
sa_aromatic_histogram_in_vitro_protein = []

sa_aromatic_neg_no_llps = []
sa_aromatic_neg_any_condensate = []
sa_aromatic_neg_in_vitro_protein = []

sa_aromatic_pos_no_llps = []
sa_aromatic_pos_any_condensate = []
sa_aromatic_pos_in_vitro_protein = []
out1.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
test.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
mean_median_mode(sa_aromatics_counts)
for i in range(0,max(determine_max_aromatic)+1):
        tally_aromatics(i)

out1.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
test.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
neg_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_neg_neighbor)+1):
        tally_neg_neighbors(i)

out1.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
test.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
pos_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_pos_neighbor)+1):
        tally_pos_neighbors(i)


out1.write("\n RS8 counts \n")
test.write("\n RS8 counts \n")
sa_aromatics_counts = []
aromatic_neighbor_counts = []
determine_max_aromatic = []
determine_max_neg_neighbor = []
determine_max_pos_neighbor = []
for line in RRM_list:
        if line[0] in RS8ident:
                tally_up(line)
sa_aromatic_histogram_no_llps = []
sa_aromatic_histogram_any_condensate = []
sa_aromatic_histogram_in_vitro_protein = []

sa_aromatic_neg_no_llps = []
sa_aromatic_neg_any_condensate = []
sa_aromatic_neg_in_vitro_protein = []

sa_aromatic_pos_no_llps = []
sa_aromatic_pos_any_condensate = []
sa_aromatic_pos_in_vitro_protein = []
out1.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
mean_median_mode(sa_aromatics_counts)
for i in range(0,max(determine_max_aromatic)+1):
        tally_aromatics(i)

out1.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
neg_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_neg_neighbor)+1):
        tally_neg_neighbors(i)

out1.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
pos_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_pos_neighbor)+1):
        tally_pos_neighbors(i)
 
        
out1.write("\n Speckles_not_RS4 \n")
test.write("\n Speckles_not_RS4 \n")
sa_aromatics_counts = []
aromatic_neighbor_counts = []
determine_max_aromatic = []
determine_max_neg_neighbor = []
determine_max_pos_neighbor = []
for line in RRM_list:
        if line[0] in RS4ident:
            pass
        else:
            a = line[19].split(',')
            c = line[21].split(',')
            if "Nuclear-speckle" in a or "Nuclear-speckle" in c:
        #when calling the function later, use an if function to choose which lines are used
                tally_up(line)
sa_aromatic_histogram_no_llps = []
sa_aromatic_histogram_any_condensate = []
sa_aromatic_histogram_in_vitro_protein = []
sa_aromatic_neg_no_llps = []
sa_aromatic_neg_any_condensate = []
sa_aromatic_neg_in_vitro_protein = []

sa_aromatic_pos_no_llps = []
sa_aromatic_pos_any_condensate = []
sa_aromatic_pos_in_vitro_protein = []
out1.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
mean_median_mode(sa_aromatics_counts)
for i in range(0,max(determine_max_aromatic)+1):
        tally_aromatics(i)

out1.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
neg_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_neg_neighbor)+1):
        tally_neg_neighbors(i)

out1.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
pos_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_pos_neighbor)+1):
        tally_pos_neighbors(i)
        
out1.write("\n Speckles_RS4 \n")
test.write("\n Speckles_RS4 \n")
sa_aromatics_counts = []
aromatic_neighbor_counts = []
determine_max_aromatic = []
determine_max_neg_neighbor = []
determine_max_pos_neighbor = []
for line in RRM_list:
        if line[0] in RS4ident:
            a = line[19].split(',')
            c = line[21].split(',')
            if "Nuclear-speckle" in a or "Nuclear-speckle" in c:
        #when calling the function later, use an if function to choose which lines are used
                tally_up(line)
sa_aromatic_histogram_no_llps = []
sa_aromatic_histogram_any_condensate = []
sa_aromatic_histogram_in_vitro_protein = []
sa_aromatic_neg_no_llps = []
sa_aromatic_neg_any_condensate = []
sa_aromatic_neg_in_vitro_protein = []

sa_aromatic_pos_no_llps = []
sa_aromatic_pos_any_condensate = []
sa_aromatic_pos_in_vitro_protein = []
out1.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
mean_median_mode(sa_aromatics_counts)
for i in range(0,max(determine_max_aromatic)+1):
        tally_aromatics(i)

out1.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
neg_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_neg_neighbor)+1):
        tally_neg_neighbors(i)

out1.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
pos_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_pos_neighbor)+1):
        tally_pos_neighbors(i)
        
        
out1.write("\n not_RS4 \n")
test.write("\n not_RS4 \n")
sa_aromatics_counts = []
aromatic_neighbor_counts = []
determine_max_aromatic = []
determine_max_neg_neighbor = []
determine_max_pos_neighbor = []
for line in RRM_list:
        if line[0] in RS4ident:
            pass
        else: 
            tally_up(line)
sa_aromatic_histogram_no_llps = []
sa_aromatic_histogram_any_condensate = []
sa_aromatic_histogram_in_vitro_protein = []
sa_aromatic_neg_no_llps = []
sa_aromatic_neg_any_condensate = []
sa_aromatic_neg_in_vitro_protein = []

sa_aromatic_pos_no_llps = []
sa_aromatic_pos_any_condensate = []
sa_aromatic_pos_in_vitro_protein = []
out1.write("Aromatics Counts \n   no_llps      any_condensate    in-vitro_verified \n")
mean_median_mode(sa_aromatics_counts)
for i in range(0,max(determine_max_aromatic)+1):
        tally_aromatics(i)

out1.write("Negative Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
neg_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_neg_neighbor)+1):
        tally_neg_neighbors(i)

out1.write("Positive Neighbors \n   no_llps    any_condensate    in-vitro_verified \n")
pos_mean_median_mode(aromatic_neighbor_counts)
for i in range(0,max(determine_max_pos_neighbor)+1):
        tally_pos_neighbors(i)
out1.close()

test.close()
