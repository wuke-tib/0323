from urllib.parse import urlencode
from urllib.request import urlopen, Request
import re
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
import json

def GENENAME_2_ACC_from_uniprot(query,gene_uniprot_outfile):
    #print(' '.join(query).replace('511145.',''))
    url = 'https://legacy.uniprot.org/uploadlists/'
    params = {
        'from': 'GENENAME',
        'to': 'ACC',
        'format': 'tab',
        'query': ' '.join(query).replace('511145.',''),
        'columns':'id,entry name,protein names,genes,organism,ec,mass,sequence,go-id,database(KEGG)'
    }
    data = urlencode(params).encode()
    request = Request(url, data)
    # Please set your email address here to help us debug in case of problems.
    contact = ""
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urlopen(request)
    page = response.read()
    outFile = open(gene_uniprot_outfile,'w') 
    namesRegex = re.compile(r'yourlist:(.*)\n')
    outFile.write(namesRegex.sub('Gene ID\n',page.decode('utf-8')))
    #print(namesRegex.sub('Protein AC\t',page.decode('utf-8')))
    outFile.close()
    
def GENEID_2_ACC_from_uniprot(query,gene_uniprot_outfile):
    #print(' '.join(query).replace('511145.',''))
    url = 'https://legacy.uniprot.org/uploadlists/'
    params = {
        'from': 'ACC+ID',
        'to': 'ACC',
        'format': 'tab',
        'query': ' '.join(query).replace('511145.',''),
        'columns':'id,entry name,protein names,genes,organism,ec,mass,sequence'
    }
    data = urlencode(params).encode()
    request = Request(url, data)
    # Please set your email address here to help us debug in case of problems.
    contact = ""
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urlopen(request)
    page = response.read()
    outFile = open(gene_uniprot_outfile,'w')  
    namesRegex = re.compile(r'yourlist:(.*)\n')
    outFile.write(namesRegex.sub('Gene ID\n',page.decode('utf-8')))
    #print(namesRegex.sub('Protein AC\t',page.decode('utf-8')))
    outFile.close()
    
def gene_uniprot_select(gene_uniprot_outfile,gene_uniprot_select_outfile):
    P_ACC_SEQ2ACC_MW=pd.read_csv(gene_uniprot_outfile,sep='\t',index_col=0)
    print('Total gene in uniprot: '+str(P_ACC_SEQ2ACC_MW.shape[0]))
    org_sta=pd.DataFrame()
    for each_org in np.unique(P_ACC_SEQ2ACC_MW['Organism']):
        org_sta.loc[each_org,'sum']=P_ACC_SEQ2ACC_MW[P_ACC_SEQ2ACC_MW['Organism']==each_org].shape[0]
    org_sta=org_sta.sort_values(by = ["sum"],ascending=False)
    
    strain_select=org_sta.index[0]
    P_ACC_SEQ2ACC_MW_select=P_ACC_SEQ2ACC_MW[P_ACC_SEQ2ACC_MW['Organism']==strain_select]
    print('Total selected gene in uniprot: '+str(P_ACC_SEQ2ACC_MW_select.shape[0]))
    P_ACC_SEQ2ACC_MW_select.to_csv(gene_uniprot_select_outfile, header=True, index=True)
    return P_ACC_SEQ2ACC_MW_select

def protein_amino_composition(gene_uniprot_select_outfile,amino_acid_information_file,seq_amino_composition_file,seq_amino_composition_MW_file):
    P_ACC_SEQ2ACC_MW=pd.read_csv(gene_uniprot_select_outfile,index_col='Gene ID')
    amino_acid_information=pd.read_csv(amino_acid_information_file,index_col=0)
    seq_amino_composition=pd.DataFrame()
    seq_amino_composition_MW=pd.DataFrame()
    
    for index, row in P_ACC_SEQ2ACC_MW.iterrows():
        mass=row['Mass'].replace(',','')
        my_seq = Seq(row['Sequence'])
        seq_amino_composition.loc[index,'seq_len']=len(my_seq)
        seq_amino_composition.loc[index,'MW']=mass
        for index2, row2 in amino_acid_information.iterrows(): 
            seq_amino_composition.loc[index,row2['abbreviation']]=my_seq.count(row2['symbol'])   
            seq_amino_composition_MW.loc[index,row2['abbreviation']]=row2['mass']*my_seq.count(row2['symbol'])/float(mass) 
    
    seq_amino_composition.to_csv(seq_amino_composition_file, header=True, index=True) 
    seq_amino_composition_MW.to_csv(seq_amino_composition_MW_file, header=True, index=True) 
    return(seq_amino_composition,seq_amino_composition_MW)

def draw_amino_numbers(seq_amino_composition_file,pngname):
    protein_amino_acid_composition=pd.read_csv(seq_amino_composition_file,index_col=0) 
    protein_amino_acid_composition=protein_amino_acid_composition.iloc[:,2:22]
    protein_amino_acid_composition=protein_amino_acid_composition.sort_index(axis = 1,ascending = True)
    
    plt.figure(figsize=(25, 10)) 
    plt.tick_params(labelsize=20)
    sns.boxplot(data=protein_amino_acid_composition)
    plt.ylabel("Numbers in protein",fontsize=22)
    plt.xlabel("Amino acid", fontsize=22)
    plt.savefig(pngname,dpi =300,bbox_inches='tight')
    plt.show()

def draw_amino_mass_in_protein(seq_amino_composition_MW_file,pngname):
    protein_amino_acid_mass=pd.read_csv(seq_amino_composition_MW_file,index_col=0) 
    protein_amino_acid_mass=protein_amino_acid_mass.iloc[:,2:22]
    protein_amino_acid_mass=protein_amino_acid_mass.sort_index(axis = 1,ascending = True)
    
    plt.figure(figsize=(25, 10)) 
    plt.tick_params(labelsize=20)
    sns.boxplot(data=protein_amino_acid_mass)
    plt.ylabel("mass in protein",fontsize=22)
    plt.xlabel("Amino acid", fontsize=22)
    plt.savefig(pngname,dpi =300,bbox_inches='tight')
    plt.show()
    
#1 g蛋白中，有多少克的氨基酸
def amino_mass_norm(seq_amino_composition_MW_file,seq_amino_composition_MW_norm_file):
    seq_amino_composition_g_g=pd.read_csv(seq_amino_composition_MW_file,index_col=0)
    seq_amino_composition_g_g_norm=seq_amino_composition_g_g/seq_amino_composition_g_g.sum(axis=1)[:,None]
    seq_amino_composition_g_g_norm.to_csv(seq_amino_composition_MW_norm_file, header=True, index=True) 
    return seq_amino_composition_g_g_norm

def draw_amino_mass(seq_amino_composition_MW_norm_file,pngname):
    ecoprot_seq_amino_composition_g_g_norm=pd.read_csv(seq_amino_composition_MW_norm_file,index_col=0) 
    ecoprot_seq_amino_composition_g_g_norm=ecoprot_seq_amino_composition_g_g_norm.sort_index(axis = 1,ascending = True)
    
    plt.figure(figsize=(25, 10)) 
    plt.tick_params(labelsize=20)
    sns.boxplot(data=ecoprot_seq_amino_composition_g_g_norm)
    plt.ylabel("Mass ratio (g/g protein)",fontsize=22)
    plt.xlabel("Amino acid", fontsize=22)
    plt.savefig(pngname,dpi =300,bbox_inches='tight')
    plt.show()

def json_write(path, dictionary_data):
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path:   The path of the JSON file that shall be written
    * dictionary: The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary_data, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)

#1 g总蛋白中，有多少克的蛋白
def protein_expression_mass_norm(proteome_file,gene_uniprot_select_outfile,protein_expression_mass_norm_outfile):
    my_pro_exp=pd.read_csv(proteome_file,index_col=0) 
    P_ACC_SEQ2ACC_MW=pd.read_csv(gene_uniprot_select_outfile,index_col='Gene ID')
    protein_expression_mass_norm={}
    my_pro_exp_T=my_pro_exp.T
    for index, row in my_pro_exp_T.iterrows():
        protein_expression_mass_norm[index]={}
        for eachgene in list(my_pro_exp_T.columns): 
            if eachgene in list(P_ACC_SEQ2ACC_MW.index):
                if isinstance(P_ACC_SEQ2ACC_MW.loc[eachgene,'Mass'],str):
                    pro_MW=float(P_ACC_SEQ2ACC_MW.loc[eachgene,'Mass'].replace(',','')) 
                else:
                    pro_MW=float(P_ACC_SEQ2ACC_MW.loc[eachgene,'Mass'].values[0].replace(',','')) 
                protein_expression_mass_norm[index][eachgene]=float(row[eachgene])*pro_MW
        total = np.sum(list(protein_expression_mass_norm[index].values()))
        protein_expression_mass_norm[index] = {k: round(v / total,6) for k, v in protein_expression_mass_norm[index].items()}
        #break
    json_write(protein_expression_mass_norm_outfile+'_exp_onecell.json', protein_expression_mass_norm)
    return protein_expression_mass_norm

def top_protein_mass_ratio(protein_expression_mass_norm_df,topnum,top_filename,pngname):    
    df_mean=pd.DataFrame(np.mean(protein_expression_mass_norm_df))
    tmp=df_mean.T    
    tmp2 = tmp.rename(index={tmp.index[0]:'bj'}) 
    tmp3=tmp2.T
    tmp4=tmp3.sort_values(by='bj',ascending=False)
    gene=np.array(tmp4.index[0:topnum])
    use_prot=protein_expression_mass_norm_df.loc[:,gene]
    use_prot.to_csv(top_filename, header=True, index=True)
    print(np.sum(use_prot.mean()))
    
    plt.figure(figsize=(25, 10)) 
    plt.tick_params(labelsize=20)
    sns.boxplot(data=use_prot)
    plt.ylabel("Mass ratio (g/g total protein)",fontsize=22)
    plt.xlabel("Protein", fontsize=22) 
    plt.xticks([]) #
    plt.savefig(pngname,dpi =300,bbox_inches='tight')
    plt.show()
    return use_prot
    
def json_load(path):
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary

#1g 总蛋白中，有多少克的氨基酸
def amino_acid_expression_mass_norm(protein_expression_mass_norm_outfile,seq_amino_composition_MW_norm_file,amino_composition_norm_onecell_outfile):
    protein_expression_mass_norm=json_load(protein_expression_mass_norm_outfile) #1g蛋白/1g总蛋白
    seq_amino_composition_MW_norm=pd.read_csv(seq_amino_composition_MW_norm_file,index_col=0)#g氨基酸/1g蛋白
    amino_acid_expression_mass_norm={}
    
    for key, value in protein_expression_mass_norm.items():
        amino_acid_expression_mass_norm[key]={}
        for eachamino in list(seq_amino_composition_MW_norm.columns): 
            amino_acid_expression_mass_norm[key][eachamino]={}
            for eachgene in value.keys(): 
                amino_acid_expression_mass_norm[key][eachamino][eachgene]=round(protein_expression_mass_norm[key][eachgene]*seq_amino_composition_MW_norm.loc[eachgene,eachamino],6)
            amino_acid_expression_mass_norm[key][eachamino]['total']=np.sum(list(amino_acid_expression_mass_norm[key][eachamino].values()))
        #break 

    json_write(amino_composition_norm_onecell_outfile, amino_acid_expression_mass_norm)
    return amino_acid_expression_mass_norm

def draw_amino_mass_in_total_protein(amino_composition_norm_onecell_df,pngname):
    twoamino=amino_composition_norm_onecell_df.sort_index(axis = 1,ascending = True)
    
    plt.figure(figsize=(25, 10)) 
    plt.tick_params(labelsize=20)
    sns.boxplot(data=twoamino)
    plt.ylabel("Mass ratio (g/g total protein)",fontsize=22)
    plt.xlabel("Amino acid", fontsize=22)
    plt.savefig(pngname,dpi =300,bbox_inches='tight')
    
    plt.show()
    
def element_composition(amino_acid_information_file,amino_composition_by_condition_outfile,element_composition_by_condition_outfile,amino_composition_norm_onecell_df_outfile):
    amino_acid_information=pd.read_csv(amino_acid_information_file,index_col='abbreviation')
    amino_composition_norm_onecell_df=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)
    df_amino=pd.DataFrame()      
    for index, row in amino_composition_norm_onecell_df.iterrows():
        N_sum=[]
        H_sum=[]
        O_sum=[]
        C_sum=[] 
        S_sum=[]
        for eachanimo in amino_composition_norm_onecell_df.columns:
            amino_total=amino_acid_information.loc[eachanimo,'formula_N']*14+amino_acid_information.loc[eachanimo,'formula_H']+\
                amino_acid_information.loc[eachanimo,'formula_O']*16+amino_acid_information.loc[eachanimo,'formula_C']*12+\
                amino_acid_information.loc[eachanimo,'formula_S']*32
            N_sum.append(amino_acid_information.loc[eachanimo,'formula_N']*14/amino_total*row[eachanimo])
            H_sum.append(amino_acid_information.loc[eachanimo,'formula_H']/amino_total*row[eachanimo])
            O_sum.append(amino_acid_information.loc[eachanimo,'formula_O']*16/amino_total*row[eachanimo])
            C_sum.append(amino_acid_information.loc[eachanimo,'formula_C']*12/amino_total*row[eachanimo])
            S_sum.append(amino_acid_information.loc[eachanimo,'formula_S']*32/amino_total*row[eachanimo])
        df_amino.loc[index,'N'] = np.sum(N_sum)
        df_amino.loc[index,'H'] = np.sum(H_sum)
        df_amino.loc[index,'O'] = np.sum(O_sum)
        df_amino.loc[index,'C'] = np.sum(C_sum)        
        df_amino.loc[index,'S'] = np.sum(S_sum)   
    df_amino=df_amino.sort_index(axis = 1,ascending = True)
    df_amino.to_csv(element_composition_by_condition_outfile, header=True, index=True)
    return df_amino

def draw_element_mass_in_total_protein(element_composition_by_condition_outfile,pngname):
    df_amino=pd.read_csv(element_composition_by_condition_outfile,index_col=0)
    df_amino=df_amino.sort_index(axis = 1,ascending = True)
    
    plt.figure(figsize=(10, 10)) 
    sns.boxplot(data=df_amino)
    plt.ylabel("Mass ratio (g/g total protein)",fontsize=20)
    plt.xlabel("Element", fontsize=20) 
    plt.savefig(pngname,dpi =300,bbox_inches='tight')
    plt.show()
    
def json_write(path, dictionary):
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path:   The path of the JSON file that shall be written
    * dictionary: The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)
        
def draw_box(data,labels,colors,positions,x_position,x_position_fmt,ylabel,pngname,xlen,ylen):
    plt.figure(figsize=(xlen, ylen)) 
    for record in range(len(data)):
        usedata=data[record]
        useposition=positions[record]
        #patch_artist=True-->箱型可以更换颜色，positions=(1,1.4,1.8)-->将同一组的三个箱间隔设置为0.4，widths=0.3-->每个箱宽度为0.3 
        bplot = plt.boxplot(usedata, patch_artist=True,labels=labels,positions=useposition,widths=0.3) 
        #将三个箱分别上色
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
            
    plt.yticks([i + 0.8 / 2 for i in x_position], x_position_fmt)
    plt.xlabel(ylabel)
    plt.grid(linestyle="--", alpha=0.3)  #绘制图中虚线 透明度0.3
    plt.legend(bplot['boxes'],labels,loc='lower right')  #绘制表示框，右下角绘制
    plt.savefig(pngname,dpi =300,bbox_inches='tight') 
    #plt.gca().invert_yaxis()
    plt.show()
    
def draw_box_vert(data,labels,colors,positions,x_position,x_position_fmt,ylabel,pngname,xlen,ylen):
    plt.figure(figsize=(xlen, ylen)) 
    for record in range(len(data)):
        usedata=data[record]
        useposition=positions[record]
        #patch_artist=True-->箱型可以更换颜色，positions=(1,1.4,1.8)-->将同一组的三个箱间隔设置为0.4，widths=0.3-->每个箱宽度为0.3 
        bplot = plt.boxplot(usedata, patch_artist=True,labels=labels,positions=useposition,widths=0.3,vert=False) 
        #将三个箱分别上色
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
            
    plt.yticks([i + 0.8 / 2 for i in x_position], x_position_fmt)
    plt.xlabel(ylabel)
    plt.grid(linestyle="--", alpha=0.3)  #绘制图中虚线 透明度0.3
    plt.legend(bplot['boxes'],labels,loc='lower right')  #绘制表示框，右下角绘制
    plt.savefig(pngname,dpi =300,bbox_inches='tight') 
    #plt.gca().invert_yaxis()
    plt.show()