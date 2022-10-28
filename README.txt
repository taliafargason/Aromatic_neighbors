-------BEFORE BEGINNING-----
Generate lists of proteins containing RS repeats of various lengths
1) modify BEFORE_BEGINNING_construct_RS_RRM_reports.py to desired domain length
2) run this script.  Currently the default constructs a list of proteins containing 4-mer tetrapeptides.  Comments will need to be modified to look for longer/shorter repeats.
-------STEPS 0 and 1-------

This section creates a list of RRM domain sequences and preliminary analysis of the overall sequence traits.

To run:

1)	Step0_construct_domain_list.py
2)	Step1_construct_RRM_list.py
3)	Copy the output (`RRM_list.txt') to the next folder
Brief description:
1)	Creates a list of all domains contained within RRM-containing proteins
2)	Creates a list of RRM domains and uses protein sequence along with domain boundaries to 
count the number of charged residues in each RRM domain and determine whether the protein 
is phase separated
Verbose Description/Obtaining and formatting downloads:
1)	Uniprot spreadsheet can be downloaded directly from the Uniprot search bar in the form of an 
Excel spreadsheet.  This script requires the protein name, domain, and protein sequence.  
Spaces were removed from the individual columns (can use the replace tool or "cntrl H" for this) 
and copied into a text file in space-delimited format (left-most column constructed then 
copy-pasted).  In this case, we used the advanced search tool on Uniprot to search for proteins 
containing RRM domains identified by manual assertion that were not fragments, that had been 
reviewed (Swissprot database, not TREMBL), and for which there was evidence at the protein 
level.


2)	 Phase separation data was downloaded from the following websites:
	a.	http://db.phasep.pro/download/
		i.	The "PhaSepDB2.0 data" zip file (contains 4 spreadsheets; 
			"phasepdbv2_mlolt_mlolt.xlsx" and "phasepdbv2_llps.xlsx" were used
	b.	http://llps.biocuckoo.cn/download.php 
		i.	Primary download at the top of the page
	c.	http://bio-comp.org.cn/llpsdb/download.php 
		i.	Natural proteins, protein(s), and protein + RNA
3)	This was used to construct the following lists that are read in separately in this portion of the script:
"Dr_Zhang_phase_separate.txt"
"phasepdb2_MLOhtp.txt"
"combo.txt"
"llps_proteinc.txt"
"llps_protein+RNAc.txt"
"drllpsc.txt"

--------STEP2-------
This section calculates surface accessibility of atoms in the proteins.  For this section, XplorNIH needs to be installed.  
The 'CalcSA' command will be run as part of the script.

1)change the directory in the script to the file path of a directory that contains only the human alphafold structures. 
For this script to work, all pdb files need to be in both the folder containing this script and a separate folder by 
themselves (here the "alphafold_in" folder)

2)in the directory with the script, include the list of RRM domains generated in step1.

3)this will generate a series of files that should then be copied to a new separate folder for the next step. To eliminate redundancy,
for proteins with multiple structures we only copied files from the F1 models over.

Downloads:
1) here, the alphafold structures of human proteins were downloaded from https://alphafold.com/download (Proteome ID UP000005640)
this involved 23,391 protein structures.  After unzipping the file, it should be verified that all 23,391 structures are present.
All 206 proteins containing RRM domains were located.
3090 files were generated for various atoms

-------STEP3-------
This step will determine which aromatic residues are surface exposed and what residues are neighboring those residues.

to run:
1)change directory to the folder containing the SA aromatic outfiles
2)modify parameters at the begining of the the distance_aromatic script as desired and run file

-------STEP4-------
This step will analyze the data generated in step3 and put it in a more palatable format.
