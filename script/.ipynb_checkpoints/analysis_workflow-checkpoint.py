# coding=utf-8
import sys
sys.path.append(r'../script/')
import pandas as pd
import argparse, os
import re
import numpy as np
import random

from protein_composition import *


def get_gene_uniprot_file(strain, result_file, tag):
    if re.search('_pro',strain):
        proteome_file='./basic_data/omics_data/%s_proteome.csv'%strain
        # proteome_file_1='./basic_data/omics_data/%s_proteome_1.csv'%strain
    elif re.search('_single_cell',strain):
        proteome_file='./basic_data/omics_data/%s_transcriptome.csv'%strain
        # proteome_file_1='./basic_data/omics_data/%s_transcriptome_1.csv'%strain
    elif re.search('_RNA_seq',strain):
        proteome_file='./basic_data/omics_data/%s_transcriptome.csv'%strain
    else:
        proteome_file='./basic_data/omics_data/%s_trans_transcriptome.csv'%strain

    if tag !="0":
        proteome_file = transcriptome_set_data(proteome_file, tag)
    
    gene_uniprot_outfile='./download_data/Gene_name2ACC_MW_%s.txt'%strain
    # return proteome_file

    my_pro_exp=pd.read_csv(proteome_file,index_col=0) 
    print(my_pro_exp.shape) 

    if re.search('Pseudomonas',strain):
        GENEID_2_ACC_from_uniprot(my_pro_exp.index,gene_uniprot_outfile)
    else:
        GENENAME_2_ACC_from_uniprot(my_pro_exp.index,gene_uniprot_outfile)

    gene_uniprot_select_outfile = result_file + 'Gene_name2ACC_MW_%s_select.txt'%strain

    P_ACC_SEQ2ACC_MW_select=gene_uniprot_select(gene_uniprot_outfile,gene_uniprot_select_outfile)

    return (proteome_file, gene_uniprot_select_outfile)

def transcriptome_set_data(proteome_file, tag):
    proteome_file_new = proteome_file.split('.csv')[0] + "_"+tag + '.csv'
    print("proteome_file_new:",proteome_file_new)
    my_pro_exp1=pd.read_csv(proteome_file,index_col=0) 

    columns_list = my_pro_exp1.columns.tolist()
    index_list = my_pro_exp1.index.tolist()

    if tag == "1":
        arr_1 = np.ones((len(index_list), len(columns_list))) # 全1矩阵
        transcriptome_df = pd.DataFrame(arr_1, index=index_list,columns=columns_list,dtype=int)

    elif tag == 'random':   # 转录组数据随机（范围：原始数据最小值-最大值）
        # 获取df中数据的最小值和最大值
        columns_max = list(my_pro_exp1.max())  # 每列最大值
        columns_min = list(my_pro_exp1.min())  # 每列最小值
        min_value = min(columns_min)
        max_value = max(columns_max)
        
        # print(round(d,2))   # 限定小数位数
        # print("随机数：",random_data)
        arr_1 = np.random.uniform(min_value, max_value, (len(index_list), len(columns_list)))  # 随机生成[min,max),index行columns列矩阵
        transcriptome_df = pd.DataFrame(arr_1, index=index_list, columns=columns_list)

    transcriptome_df.to_csv(proteome_file_new, header=True, index=True)

    return proteome_file_new


def get_amino_composition_from_protein(strain, proteome_file, gene_uniprot_select_outfile, result_file):
    # Get amino acid composiotion from protein sequence
    # gene_uniprot_select_outfile=result_file + 'Gene_name2ACC_MW_%s_select.txt'%strain
    amino_acid_information_file='./basic_data/amino_acid_information.csv'
    seq_amino_composition_file=result_file + 'seq_amino_composition_%s.csv'%strain
    seq_amino_composition_MW_file=result_file + 'seq_amino_composition_g_g_%s.csv'%strain
    [seq_amino_composition,seq_amino_composition_MW]=protein_amino_composition(gene_uniprot_select_outfile,amino_acid_information_file,seq_amino_composition_file,seq_amino_composition_MW_file)


    pngname=result_file + 'protein_amino_acid_composition_boxplot_%s.png'%strain
    draw_amino_numbers(seq_amino_composition_file,pngname)

    # Amino acid composition condsider protein sequence (normalized 1g protein)
    seq_amino_composition_MW_norm_file=result_file + 'seq_amino_composition_g_g_norm_%s.csv'%strain
    seq_amino_composition_g_g_norm=amino_mass_norm(seq_amino_composition_MW_file,seq_amino_composition_MW_norm_file)


    pngname=result_file + 'protein_amino_acid_composition_g_g_boxplot_%s.png'%strain
    draw_amino_mass(seq_amino_composition_MW_norm_file,pngname)

    # Amino acids composition of each protein (g / g protein) consider expression level under different conditions

    protein_expression_mass_norm_outfile=result_file + '%s_exp_onecell.json'%strain
    protein_expression_mass_norm_json=protein_expression_mass_norm(proteome_file,gene_uniprot_select_outfile,protein_expression_mass_norm_outfile)

    
    protein_expression_mass_norm_data=json_load(protein_expression_mass_norm_outfile) 
    protein_expression_mass_norm_df=pd.DataFrame(protein_expression_mass_norm_data).T

    #delete outliers
    if re.search('_pro',strain):
        #condition: b3916(KO) 
        protein_expression_mass_norm_df=protein_expression_mass_norm_df.drop(['b3916(KO)'])
    elif re.search('DH1',strain):
        #condition: DH1.MD018.FPP.mevalonate-pw_IN_T5223#DH1.MD018.none.mevalonate-pw_IN;FPP-pw_IN_T5222#DH1.MD018.FPP.mevalonate-pw_IN_T5215;
        protein_expression_mass_norm_df=protein_expression_mass_norm_df.drop(['DH1.MD018.FPP.mevalonate-pw_IN_T5223','DH1.MD018.none.mevalonate-pw_IN;FPP-pw_IN_T5222','DH1.MD018.FPP.mevalonate-pw_IN_T5215'])    
    elif re.search('MG1655',strain):
        #condition: MG1655.MD097.O2-starvation;CORM-3.na_WT_T6149#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6128#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6145#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6126#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6153#MG1655.MD018.Cfs;Mcn.na_WT_T0040#MG1655.MD097.O2-starvation.na_WT_T6125#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6155#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6151#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6171#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6130#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6132#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6134#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6136#MG1655.MD018.Mcn.na_WT_T0029#MG1655.MD097.O2-starvation.na_WT_T6214#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6167#MG1655.MD010.lactate-shift.na_WT_T0840#MG1655.MD003.none.na_WT_T1588#MG1655.MD003.none.na_WT_T1590#MG1655.MD003.none.na_WT_T1593#MG1655.MD097.O2-starvation.na_WT_T6184#MG1655.MD097.O2-starvation.na_WT_T6216#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6165#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6169#MG1655.MD097.O2-starvation;CORM-3.na_WT_T6143#MG1655.MD090.none.na_WT_T8309#MG1655.MD010.lactate-shift.na_WT_T0841#MG1655.MD091.none.na_WT_T8318#MG1655.MD116.none.na_WT_T1713#MG1655.MD027.none.na_WT_T1941#MG1655.MD027.none.na_WT_T1942#MG1655.MD027.none.na_WT_T1943#MG1655.MD097.O2-starvation.na_WT_T6178
        protein_expression_mass_norm_df=protein_expression_mass_norm_df.drop(['MG1655.MD097.O2-starvation;CORM-3.na_WT_T6149','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6128','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6145','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6126','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6153','MG1655.MD018.Cfs;Mcn.na_WT_T0040','MG1655.MD097.O2-starvation.na_WT_T6125','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6155','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6151','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6171','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6130','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6132','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6134','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6136','MG1655.MD018.Mcn.na_WT_T0029','MG1655.MD097.O2-starvation.na_WT_T6214','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6167','MG1655.MD010.lactate-shift.na_WT_T0840','MG1655.MD003.none.na_WT_T1588','MG1655.MD003.none.na_WT_T1590','MG1655.MD003.none.na_WT_T1593','MG1655.MD097.O2-starvation.na_WT_T6184','MG1655.MD097.O2-starvation.na_WT_T6216','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6165','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6169','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6143','MG1655.MD090.none.na_WT_T8309','MG1655.MD010.lactate-shift.na_WT_T0841','MG1655.MD091.none.na_WT_T8318','MG1655.MD116.none.na_WT_T1713','MG1655.MD027.none.na_WT_T1941','MG1655.MD027.none.na_WT_T1942','MG1655.MD027.none.na_WT_T1943','MG1655.MD097.O2-starvation.na_WT_T6178'])    
    elif re.search('BW25113',strain):
        #condition: BW25113.MD106.none.na_WT_T5490 
        protein_expression_mass_norm_df=protein_expression_mass_norm_df.drop(['BW25113.MD106.none.na_WT_T5490'])        

    topnum=250
    top_filename=result_file + ''+strain+'_'+str(topnum)+'_expdata.csv'
    pngname=result_file + 'protein_top%s_exp%s.png'%(topnum,strain)
    use_prot=top_protein_mass_ratio(protein_expression_mass_norm_df,topnum,top_filename,pngname)


    # Amino acids composition of each condaition (g / g total protein) consider expression level
    amino_composition_norm_onecell_outfile=result_file + 'amino_composition_g_g_norm_onecell_%s.json'%strain
    amino_acid_expression_mass_norm_json=amino_acid_expression_mass_norm(protein_expression_mass_norm_outfile,seq_amino_composition_MW_norm_file,amino_composition_norm_onecell_outfile)


    # draw figure
    amino_composition_norm_onecell=json_load(amino_composition_norm_onecell_outfile)
    amino_composition_norm_onecell_df=pd.DataFrame()
    for key, value in amino_composition_norm_onecell.items():
        for key2 in value.keys():
            amino_composition_norm_onecell_df.loc[key,key2]=value[key2]['total']
            
    amino_composition_norm_onecell_df_outfile=result_file + 'amino_composition_proteome_by_condition_%s.csv'%strain
    amino_composition_norm_onecell_df.to_csv(amino_composition_norm_onecell_df_outfile, header=True, index=True) 

    pngname=result_file + '%s_twenty_amino_condition_boxplot.png'%(strain)
    draw_amino_mass_in_total_protein(amino_composition_norm_onecell_df,pngname)



def analysis_main():
    parser = argparse.ArgumentParser(description='Analysis Workfolw!')
    parser.add_argument('-i', '--strain', dest='strain', nargs='+', help="list")
    parser.add_argument('-o', '--file_Folder', type=str, default='./analysis_result/')
    parser.add_argument('-tag', '--tag', type=str)


    args = parser.parse_args()
    strain_list = args.strain
    file_Folder = args.file_Folder
    tag = args.tag   # 标签"0","1","random"


    print(strain_list)
    for strain in strain_list:
        print("start："+ strain)
        result_file = file_Folder + '%s/'%strain
        if not os.path.exists(result_file):
            os.makedirs(result_file)
    
        (proteome_file,gene_uniprot_select_outfile) = get_gene_uniprot_file(strain, result_file, tag)
        get_amino_composition_from_protein(strain, proteome_file, gene_uniprot_select_outfile, result_file)
        print("*******************")

if __name__ == '__main__':
    analysis_main()
