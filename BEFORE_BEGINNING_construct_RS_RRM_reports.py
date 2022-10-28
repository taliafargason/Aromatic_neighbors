from Bio import SeqIO
import re
inputfile='reviewed_only_organism9606_existence1_fragmentfalse_reviewedtrue_2022.08.17.fasta'
search_motif=[]
for a in ['S']:
    for b in ['R']:
        for c in ['S']:
            for d in ['R']:
              # for e in ['S']:
               #     for f in ['R']:
               #         for g in ['S']:
                #            for h in ['R']:      #comment these back in to look for longer repeats!   
                                    if a+b+c+d not in search_motif:#add e-h to look for longer repeats!
                                        search_motif.append(a+b+c+d)#add e-h to look for longer repeats!
                                    
                                    if b+a+d+c not in search_motif:#add e-h to look for longer repeats!
                                        search_motif.append(b+a+d+c)#add e-h to look for longer repeats!
phase_list=[]
rrmlist = open("proteins_w_rrm_domains.txt","r")
RRMlist = [re.split('\s+',line)[0] for line in rrmlist]
with open ('combo.txt') as p:
    for line in p:
        if line[:-1] not in phase_list:
            phase_list.append(line[:-1])
        

g = open('RRM_RS_report_%s.txt'%search_motif[0],'w')
seq=list(SeqIO.parse(inputfile, "fasta"))

#list of proteins containing RS8 and SR8
name_list=[]
count_list=[]
assessment1 = []
assessment2 = []
#search_motif=['ERER', 'RERE', 'DRER', 'ERDR', 'EKDR', 'DKDK','KDKD', 'ERDK']
for i in range(len(seq)):
    for j in range(len(search_motif)):
        if (search_motif[j] in seq[i].seq) and (seq[i].id.split('|')[1] not in name_list):
            #if the protein contains RS and is not already in namelist, append uniprotID to namelist
            name_list.append(seq[i].id.split('|')[1])
            temp=[]
            for k in range(len(search_motif)):
                temp.append(seq[i].seq.count(search_motif[k]))
            count=max(temp)
            count_list.append(count)
            #number times search motif appears
        elif (search_motif[0]) not in seq[i].seq and search_motif[1] not in seq[i].seq:
            ID = seq[i].id.split('|')[1]
            if ID not in assessment2:
                assessment2.append(ID)
                if ID in phase_list and ID in RRMlist:
                    assessment1.append("%s %s %s %s"%(ID,0,1,1))
                    g.write("%s %s %s %s \n"%(ID,0,1,1))
                elif ID in phase_list:
                    assessment1.append("%s %s %s %s"%(ID,0,1,0))
                    g.write("%s %s %s %s \n"%(ID,0,1,0))
                elif ID in RRMlist:
                    assessment1.append("%s %s %s %s"%(ID,0,0,1))
                    g.write("%s %s %s %s \n"%(ID,0,0,1))
                else:
                    assessment1.append("%s %s %s %s"%(ID,0,0,0))
                    g.write("%s %s %s %s \n"%(ID,0,0,0))
judge=[]
RRM_RS_list = []

for i in range(len(name_list)):
    if name_list[i] in phase_list:
        judge.append(1)
        assessment1.append(1)
        if name_list[i] in RRMlist:
            RRM_RS_list.append(1)
        else:
            RRM_RS_list.append(0)
    else:
        judge.append(0)
        if name_list[i] in RRMlist:
            RRM_RS_list.append(1)
        else:
            RRM_RS_list.append(0)



f=open('report_%s.txt' % search_motif[0] ,'w')

f.write ('protein_ID occurence in_condensate\n')
print (len(name_list))

for i in range(len(name_list)):
    print (name_list[i], count_list[i], judge[i], file=f)
    print (name_list[i], count_list[i], judge[i], RRM_RS_list[i],file=g)
    assessment1.append("%s %s %s %s"%(name_list[i],count_list[i],judge[i],RRM_RS_list[i]))

f.close()

confirm_list=list(set(name_list) & set(phase_list))

f=open('confirmed_%s.txt' % search_motif[0], 'w')
#g=open('RRMconfirmed_%s.txt' % search_motif[0], 'w')
for i in range (len(confirm_list)):
    f.write("\n%s" %confirm_list[i])
    if confirm_list[i] in RRMlist:
        f.write(" RRM")
f.close()
h = open('b_RRM_RS_report_%s.txt'%search_motif[0],'w')
assess = open("%s_assessment_out.txt"%search_motif[0],"w")
RS_RRM_LLPs = 0
RS_LLPs = []
RS_noLLPs = []
RS_RRM_noLLPs = 0
RS_noRRM_LLPs = 0
RS_noRRM_noLLPs = 0
noRS_RRM_LLPs = 0
noRS_RRM_noLLPs = 0
noRS_noRRM_LLPs = 0
noRS_noRRM_noLLPs = 0
fell_thru_cracks = []
#for i in max number of RS domains (max a[1])
for line in assessment1:
    try:
        a = re.split('\s+',line) 
        h.write("%s \n"%str(a))#a[1] = RS, a[2] = llps, a[3] = RRM
        if int(a[1])>0 and int(a[2])==1 and int(a[3])==1:
            RS_RRM_LLPs +=1
            RS_LLPs.append(a)
        elif int(a[1])>0 and int(a[2])==0 and int(a[3])==1:
            RS_RRM_noLLPs +=1
            RS_noLLPs.append(a)
        elif int(a[1])>0 and int(a[2])==1 and int(a[3])==0:
            RS_noRRM_LLPs +=1
            RS_LLPs.append(a)
        elif int(a[1])>0 and int(a[2])==0 and int(a[3])==0:
            RS_noRRM_noLLPs +=1
            RS_noLLPs.append(a)
        elif int(a[1])==0 and int(a[2])==1 and int(a[3])==1:
            noRS_RRM_LLPs +=1
        elif int(a[1])==0 and int(a[2])==0 and int(a[3])==1:
            noRS_RRM_noLLPs +=1
        elif int(a[1])==0 and int(a[2])==1 and int(a[3])==0:
            noRS_noRRM_LLPs +=1
        elif int(a[1])==0 and int(a[2])==0 and int(a[3])==0:
            noRS_noRRM_noLLPs +=1
        else:
            fell_thru_cracks.append(a)
    except TypeError:
        print(line)
assess.write("+%s+RRM\n in_condensates not_in_condensates\n %s %s\n\n"%(search_motif[0],RS_RRM_LLPs,RS_RRM_noLLPs))
assess.write("+%s-RRM\n in_condensates not_in_condensates\n %s %s\n\n"%(search_motif[0],RS_noRRM_LLPs,RS_noRRM_noLLPs))
assess.write("-%s+RRM\n in_condensates not_in_condensates\n %s %s\n\n"%(search_motif[0],noRS_RRM_LLPs,noRS_RRM_noLLPs))
assess.write("-%s-RRM\n in_condensates not_in_condensates\n %s %s\n\n"%(search_motif[0],noRS_noRRM_LLPs,noRS_noRRM_noLLPs))
g.close()
h.close()
assess.close()
