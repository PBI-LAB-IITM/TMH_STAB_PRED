#!/usr/bin/python3

from Bio.PDB import *
import sys
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis 
import pandas as pd
import os
import mdtraj as md
import numpy as np
import glob
import shutil
from tabulate import tabulate

def repair(p):

	t="./foldx_20231231 --command=RepairPDB --pdb-dir=./PDB/ --pdb="+p+".pdb --output-dir=./PDB/"    
	out = os.popen(t)
	x=out.read()
	return()

def free_energy_membrane(p):
	s="./foldx_20231231 --command=Stability --pdb-dir=./PDB/ --pdb="+p+"_Repair.pdb --output-dir=./PDB/" 
	out = os.popen(s)
	x=out.read()
	#os.system(s)
	o='./PDB/'+p+'_Repair_0_ST.fxout'
	f=open(o,'r')
	lines=f.read()
	y=lines.split('\t')
	return(y)
def PDB_Read(x):   
    f=open(x,'r')
    lines=f.readlines()
    line=[]
    for i in lines:
        if i.split()[0]=='ATOM' and i.split()[2]=='CA':
           line.append(i)
    Res=[]
    for i in line:
        Res.append(split_pdb(i)[3])
    return(Res)
def split_pdb(line):
            splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
            return(splitted_line)    
def roundd(a):
  a=float(a)
  return(round(a,3)) 
def abs_contact_order(xyz, atoms_residue, cutoff_nm=.8):
    """Return the absolute contact order."""
    contact_count = 0
    seq_distance_sum = 0
    cutoff_2 = cutoff_nm*cutoff_nm
    N = len(atoms_residue)
    for i in range(N):
        for j in range(N):
            seq_dist = atoms_residue[j] - atoms_residue[i]
            if seq_dist > 0:
                d = xyz[j] - xyz[i]
                if np.dot(d, d) < cutoff_2:
                    seq_distance_sum += seq_dist 
                    contact_count += 1


    if contact_count==0.:
        print("Warning, no contacts found!")
        return 0.

    #print("contact_count: ", contact_count)
    #print("seq_distance_sum: ", seq_distance_sum)
    return seq_distance_sum/float(contact_count)   

#### Write code to copy all the files into PDB folder inside src ##################  
entries=[]
path="../input/"
for x in os.listdir(path):
    if x.endswith(".pdb"):
        entries.append(x)
if(len(entries)==0):
 print('Enter Proper PDB File or PDB File missing in input folder/n')
 print('Check Tutorial for naming and other information /n')
#print(entries)
else:
 source = '../input/'
 destination = './PDB'
 allfiles=[]
 for i in range(0,len(entries)):
    allfiles.append(glob.glob(os.path.join(source, entries[i]), recursive=True))
 #print(allfiles)
# iterate on all files to move them to destination folder
 for file_p in allfiles:
    file_path=file_p[0]
    dst_path = os.path.join(destination, os.path.basename(file_path))
    #shutil.copy(file_path, dst_path)
    shutil.move(file_path, dst_path)


#####################################################

#####################################################
 AA={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q',
      'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
      'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

 Aminoacid=['A','R','N','D','C','E','Q','G','H','I','L','K',
                             'M','F','P','S','T','W','Y','V']
 df=pd.read_csv('aaindex.csv')
         
 header=[]
 header.append('Name')
 for m in range(2,len(df.columns)):
    head=df.iloc[0,m]
    header.append(head)
 header.append('Residues')            
 df1=pd.DataFrame(columns=header)
 print("The Input PDB :  ")
 for a in entries:
    print(a)
    CP=[]
    CP.append(a)
    for m in range(2,len(df.columns)):
        head=df.iloc[0,m]
        test_keys=list(df.iloc[2:,1])
        test_values=list(df.iloc[2:,m])
        res = {test_keys[i]: test_values[i] for i in range(len(test_keys))}

        Complete=PDB_Read('./PDB/'+a)
        Res1=[]
        for i in Complete:
            if(i!='UNK'):
                Res1.append(AA[i])
        s=0
        for i in range(0,len(Res1)):
            if(Res1[i]!='UNK'):
                vv=float(res[Res1[i]])
                s=s+vv
        if(len(Res1)!=0):
    
            CP.append(s/len(Res1))
            CP.append(len(Res1))
        else:
            CP.append(0)
    df1.loc[len(df1)]=CP
 #df1.to_csv('file1.csv')
 print('\n Sequence feature calculation Done !!!!!! \n')    
##################################################################################################################################################################################
 df2=pd.DataFrame()
       
 df2['entry']=entries 
 rel=[]

 for i in entries:
        traj = md.load('./PDB/'+i)

        seq_atoms = np.array([a.residue.resSeq for a in traj.top.atoms], dtype=int)

        abs_co = abs_contact_order(traj.xyz[0], seq_atoms, cutoff_nm=0.80)

        rel.append(abs_co/traj.n_residues)
        
 df2['Relative Contact Order']=rel
 #df2.to_csv('file2.csv')
 print('\n Contact Order calculation Done !!!!!! \n')
##########################################################################################################################################################################   
 df3=pd.DataFrame(columns=['entry','Total','BackHbond','SideHbond','Energy_VdW',
    'Electro','Energy_SolvP','Energy_SolvH','Energy_vdwclash','Entropy_sidec','Entropy_mainc','water bonds','dummy','cis_bond','energy_torsion','backbone_vdwclash','helix dipole','loop_entropy','disulfide','kn electrostatic','partial covalent interactions','Energy_Ionisation','Entropy Complex','count'])
 for i in entries:
    m=i.split(".")[0] 
    l=repair(m)
 for i in entries:
    m=i.split(".")[0]
    k=free_energy_membrane(m)
    df3.loc[len(df)]=k
 #df3.to_csv('file3.csv')   
##################################################################################################################################################################


 final=pd.DataFrame()
 final['entry']=df1['Name']
 final['Residues']=df1['Residues']
 final['Feature1']=df1['WIMW960101']
 final['Feature2']=df2['Relative Contact Order']
 final['Feature3']=float(df3['Total'])
 final['Feature4']=float(df3['BackHbond'])
 final['Feature3']=(final['Feature3'])/(final['Residues'])
 final['Feature4']=(final['Feature4'])/(final['Residues'])
 
 final['Predicted ΔG H2O kcal/mol']=-98.27+(15.98*(final['Feature1']))+(-35.21*(final['Feature2']))+(7.83*(final['Feature3']))+(-36.15*(final['Feature4'])) 
 final.to_csv('final.csv')
 table=[]
 print('\n Calculation Done !!!!!! \n')
 print('\n The following output was stored in final.csv \n')
 
 out = os.popen('sh run.sh')
 for i in range(len(final)):
      table.append(list(final.iloc[i,:]))
 print(tabulate(table,headers=list(final.columns)))     






          
