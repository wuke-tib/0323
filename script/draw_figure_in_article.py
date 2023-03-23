import sys
sys.path.append(r'../script/')
from protein_composition import *
import pandas as pd
import numpy as np
import re
import os
import argparse
import openpyxl
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from scipy import stats

# 判断文件是否存在:
def isFileExists(fileDir):
    if not os.path.exists(fileDir):
        os.makedirs(fileDir)

'''1 各物种蛋白MW比较'''
def Comparison_of_protein_MW_by_species(strain_list1,strain_list2,fig_file,analysis_result):
    MW_total=pd.DataFrame()

    Species_Number = len(strain_list1)

    createVar = locals()
    newindex=0

    for i in range(0, Species_Number):
        createVar['protein_amino_acid_composition_'+strain_list1[i]]=pd.read_csv('%s/%s/seq_amino_composition_%s.csv'%(analysis_result,strain_list1[i],strain_list1[i]),index_col=0)
        createVar[strain_list1[i]+'_MW']=pd.DataFrame(createVar['protein_amino_acid_composition_'+strain_list1[i]]['MW']).sort_values(by='MW')

        for index, row in createVar[strain_list1[i]+'_MW'].iterrows():
            MW_total.loc[newindex,'strain']=strain_list2[i] 
            MW_total.loc[newindex,'value']=row['MW']/1000
            newindex+=1
    # print(createVar['protein_amino_acid_composition_MG1655']) # df
    # print(createVar['MG1655_MW']) # df

    MW_total.to_csv('%s/MW_total.csv'%analysis_result, sep=',', header=True, index=False) 

    pngname='%s/protein_MW_violinplot_total.png'%fig_file
    #plt.rcParams['font.style'] ='italic' 
    plt.figure(figsize=(20, 10)) 
    plt.tick_params(labelsize=25)
    sns.violinplot(x='strain', y='value',data=MW_total)
    plt.ylabel("MW (kDa)",fontsize=22)
    plt.xlabel("", fontsize=22)

    ax=plt.gca() #获得坐标轴的句柄
    labels = ax.get_xticklabels()
    [label.set_fontstyle('italic') for label in labels]

    #plt.yticks((100,200,300,400,500,600), ('10e4','20e4','30e4','40e4','50e4','60e4'))
    plt.savefig(pngname,dpi =300,bbox_inches='tight')
    plt.show()
    #plt.rcParams['font.style'] ='normal'

    # strain_list2=['E. coli','B. subtilis','S. cerevisiae','P. aeruginosa','C. glutamicum']
    mw_stat_std=pd.DataFrame()
    for each_i in strain_list2:
        data= MW_total[MW_total['strain']==each_i]['value']
        mw_stat_std.loc[each_i,'mean']=np.mean(data)
        mw_stat_std.loc[each_i,'sd']=np.std(data)
        mw_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)

    print('mw_stat_std:')
    print(mw_stat_std)
    mw_stat_std.to_excel('%s/mw_stat_std.xlsx'%fig_file,index = True)

'''2.1 蛋白组-蛋白数目'''
def Proteome_Number_of_proteins(fig_file,strain,analysis_result):
    protein_amino_acid_composition=pd.read_csv('%s/%s/seq_amino_composition_%s.csv'%(analysis_result,strain,strain),index_col=0) 
    protein_amino_acid_composition=protein_amino_acid_composition.iloc[:,2:24]
    protein_amino_acid_composition=protein_amino_acid_composition.sort_index(axis = 1,ascending = True)

    pngname='%s/protein_amino_acid_composition_boxplot_%s.png'%(fig_file,strain)

    plt.figure(figsize=(25, 10)) 
    plt.tick_params(labelsize=20)
    sns.boxplot(data=protein_amino_acid_composition)
    plt.ylabel("Numbers",fontsize=22)
    plt.xlabel("Amino acid", fontsize=22)
    plt.savefig(pngname,dpi =300,bbox_inches='tight')
    # plt.show()

    amino_list=protein_amino_acid_composition.columns.values.tolist()
    amino_num_stat_std=pd.DataFrame()
                
    for each_i in amino_list:
        data=protein_amino_acid_composition[each_i]
        amino_num_stat_std.loc[each_i,'mean']=np.mean(data)
        amino_num_stat_std.loc[each_i,'sd']=np.std(data)
        amino_num_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)
    
    print('amino_num_stat_std:')
    print(amino_num_stat_std)
    amino_num_stat_std.to_excel('%s/amino_num_stat_std.xlsx'%fig_file,index = True)

'''2.2 蛋白组-各蛋白氨基酸质量比'''
def Proteome_Amino_acid_mass_ratio_of_each_protein(fig_file,strain,analysis_result):
    ecoprot_seq_amino_composition_g_g_norm=pd.read_csv('%s/%s/seq_amino_composition_g_g_norm_%s.csv'%(analysis_result,strain,strain)) 
    ecoprot_seq_amino_composition_g_g_norm=ecoprot_seq_amino_composition_g_g_norm.iloc[:,1:21]
    ecoprot_seq_amino_composition_g_g_norm=ecoprot_seq_amino_composition_g_g_norm.sort_index(axis = 1,ascending = True)

    pngname='%s/protein_amino_acid_composition_g_g_boxplot_%s.png'%(fig_file,strain)
    plt.figure(figsize=(25, 10)) 
    plt.tick_params(labelsize=20)
    sns.boxplot(data=ecoprot_seq_amino_composition_g_g_norm)
    plt.ylabel("Mass ratio (g/g protein)",fontsize=22)
    plt.xlabel("Amino acid", fontsize=22)
    plt.savefig(pngname,dpi =300,bbox_inches='tight')

    # plt.show()

    amino_list=ecoprot_seq_amino_composition_g_g_norm.columns.values.tolist()
    amino_mass_stat_std=pd.DataFrame()
                
    for each_i in amino_list:
        data=ecoprot_seq_amino_composition_g_g_norm[each_i]
        amino_mass_stat_std.loc[each_i,'mean']=np.mean(data)
        amino_mass_stat_std.loc[each_i,'sd']=np.std(data)
        amino_mass_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)
        
    print('amino_mass_stat_std:')
    print(amino_mass_stat_std)
    amino_mass_stat_std.to_excel('%s/amino_mass_stat_std.xlsx'%fig_file,index = True)

'''2.3 蛋白组-总氨基酸的质量比'''
def Proteome_Mass_ratio_of_total_amino_acids(fig_file,strain,analysis_result):
    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_by_condition_%s.csv'%(analysis_result,strain,strain)
    amino_composition_norm_onecell_df=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)

    pngname_boxplot='%s/%s_twenty_amino_condition_boxplot.png'%(fig_file,strain)
    pngname_violinplot='%s/%s_twenty_amino_condition_violinplot.png'%(fig_file,strain)
    draw_amino_mass_in_total_protein(amino_composition_norm_onecell_df,pngname_boxplot,pngname_violinplot)

    print(np.sum(np.mean(amino_composition_norm_onecell_df)))

    print(strain)
    amino_composition_norm_onecell_df=amino_composition_norm_onecell_df.sort_index(axis = 1,ascending = True)
    amino_list=amino_composition_norm_onecell_df.columns.values.tolist()
    amino_mass_norm_stat_std=pd.DataFrame()
                
    for each_i in amino_list:
        data=amino_composition_norm_onecell_df[each_i]
        amino_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        amino_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        amino_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)
        
    print('amino_mass_norm_stat_std:')
    print(amino_mass_norm_stat_std)
    amino_mass_norm_stat_std.to_excel('%s/amino_mass_norm_stat_std_amino_composition_proteome_by_condition.xlsx'%fig_file,index = True)

'''2.4 蛋白组-topnum个蛋白质量分布情况'''
def Proteome_Protein_Distribution(topnum,fig_file,strain,analysis_result):
    protein_expression_mass_norm_outfile='%s/%s/%s_exp_onecell.json'%(analysis_result,strain,strain)
    protein_expression_mass_norm=json_load(protein_expression_mass_norm_outfile) 
    protein_expression_mass_norm_df=pd.DataFrame(protein_expression_mass_norm).T
    # protein_expression_mass_norm_df=protein_expression_mass_norm_df.drop(['b3916(KO)'])
    
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

    top_filename = os.path.join(analysis_result,strain,strain+'_'+str(topnum)+'_expdata.csv')
    pngname='%s/protein_top%s_exp%s.png'%(fig_file,topnum,strain)
    use_prot=top_protein_mass_ratio(protein_expression_mass_norm_df,topnum,top_filename,pngname)

    print(topnum/protein_expression_mass_norm_df.shape[1])

    # 蛋白质量分布标准差
    protein_list=use_prot.columns.values.tolist()
    protein_mass_norm_stat_std=pd.DataFrame()
                
    for each_i in protein_list:
        data=use_prot[each_i]
        protein_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        protein_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        protein_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)

    protein_mass_norm_stat_std.to_excel('%s/protein_mass_norm_stat_std_exp_onecell.xlsx'%fig_file,index = True)
    print('protein_mass_norm_stat_std:')
    print(protein_mass_norm_stat_std)

    return use_prot,protein_mass_norm_stat_std

'''2.5 蛋白组-topnum个蛋白的氨基酸总质量'''
def Proteome_Total_amino_acid_mass_of_protein(fig_file,strain,analysis_result,topnum,use_prot,protein_mass_norm_stat_std):
    protein_expression_mass_norm_json=use_prot.T.to_json()
    top_filename_json = os.path.join(analysis_result,strain,strain+'_'+str(topnum)+'_expdata.json')
    with open(top_filename_json, "w", encoding="utf-8") as f:
        f.write(protein_expression_mass_norm_json)
        
    seq_amino_composition_MW_norm_file='%s/%s/seq_amino_composition_g_g_norm_%s.csv'%(analysis_result,strain,strain)
    amino_composition_norm_onecell_outfile='%s/%s/amino_composition_g_g_norm_onecell_%s_top%sprotein.json'%(analysis_result,strain,strain,topnum)
    amino_acid_expression_mass_norm_json=amino_acid_expression_mass_norm(top_filename_json,seq_amino_composition_MW_norm_file,amino_composition_norm_onecell_outfile)


    # draw figure
    amino_composition_norm_onecell=json_load(amino_composition_norm_onecell_outfile)
    amino_composition_norm_onecell_df=pd.DataFrame()
    for key, value in amino_composition_norm_onecell.items():
        for key2 in value.keys():
            amino_composition_norm_onecell_df.loc[key,key2]=value[key2]['total']
            
    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_top%sprotein_%s.csv'%(analysis_result,strain,topnum,strain)
    amino_composition_norm_onecell_df.to_csv(amino_composition_norm_onecell_df_outfile, header=True, index=True) 

    pngname_boxplot='%s/%s_twenty_amino_condition_boxplot_top%s_protein.png'%(fig_file,strain,topnum)
    pngname_violinplot='%s/%s_twenty_amino_condition_violinplot_top%s_protein.png'%(fig_file,strain,topnum)
    draw_amino_mass_in_total_protein(amino_composition_norm_onecell_df,pngname_boxplot,pngname_violinplot)

    print(np.sum(np.mean(amino_composition_norm_onecell_df)))

    # 氨基酸质量标准差
    amino_composition_norm_onecell_df=amino_composition_norm_onecell_df.sort_index(axis = 1,ascending = True)

    amino_list=amino_composition_norm_onecell_df.columns.values.tolist()
    amino_mass_norm_stat_std=pd.DataFrame()
                
    for each_i in amino_list:
        data=amino_composition_norm_onecell_df[each_i]
        amino_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        amino_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        amino_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)

    amino_mass_norm_stat_std.to_excel('%s/amino_mass_norm_stat_std_amino_composition_proteome_topnum.xlsx'%fig_file,index = True)
    print('amino_mass_norm_stat_std:')
    print(amino_mass_norm_stat_std)

    # two independent-sample t-test：
    # 医学研究中常用于完全随机设计两样本均数的比较，即将受试对象完全随机分配到两个不同处理组，
    # 研究者关心的是两样本均数所代表的两总体均数是否不等。
    mw_stat=pd.DataFrame()

    data1=protein_mass_norm_stat_std['sd/mean']
    data2=amino_mass_norm_stat_std['sd/mean']
    mw_stat.loc['exp2amino','mean_in_exp']=np.mean(data1)
    mw_stat.loc['exp2amino','mean_in_AA']=np.mean(data2)

    if stats.levene(data1,data2).pvalue>0.5:#如果返回结果的p值远大于0.05，那么我们认为两总体具有方差齐性。
        print('Standard independent 2 sample test')
        mw_stat.loc['exp2amino','p_value']=stats.ttest_ind(data1,data2).pvalue#独立样本T检验
    else:
        print('Welch’s t-test')
        mw_stat.loc['exp2amino','p_value']=stats.ttest_ind(data1,data2,equal_var=False).pvalue
    mw_stat=mw_stat.sort_values(by='p_value',ascending = False)
    print(mw_stat)

'''2.6 蛋白组-topnum个的蛋白种表达变化比较大的蛋白'''
def Proteome_High_protein_expression(topnum,fig_file,strain,analysis_result,use_prot):
    protein_expression_mass_norm_outfile='%s/%s/%s_exp_onecell.json'%(analysis_result,strain,strain)
    protein_expression_mass_norm=json_load(protein_expression_mass_norm_outfile) 
    protein_expression_mass_norm_df=pd.DataFrame(protein_expression_mass_norm).T

    #变化范围大的蛋白
    # diff_max_min=np.max(protein_expression_mass_norm_df,axis = 0)-np.min(protein_expression_mass_norm_df,axis = 0)
    # diff_max_min_sort=diff_max_min.sort_values(ascending = False)
    diff_max_min_norm=(np.max(protein_expression_mass_norm_df,axis = 0)-np.min(protein_expression_mass_norm_df,axis = 0))/np.max(protein_expression_mass_norm_df,axis = 0)
    diff_max_min_sort=diff_max_min_norm.sort_values(ascending = False)

    max_gene_list=list(use_prot.columns)
    #diff_max_min_list=list(diff_max_min_sort[0:topnum].index)
    # diff_max_min_list=list(diff_max_min_sort[diff_max_min_sort>0.3].index)
    diff_max_min_list=list(diff_max_min_sort[diff_max_min_sort>0.5].index)

    intersection_list=list(set(max_gene_list).intersection(set(diff_max_min_list)))
    print(len(intersection_list))
    print(len(intersection_list)/protein_expression_mass_norm_df.shape[1])

    topnum=len(intersection_list)
    top_filename = os.path.join(analysis_result,strain,strain+'_'+str(topnum)+'_expdata.csv')
    pngname='%s/protein_top%s_exp%s.png'%(fig_file,topnum,strain)
    diff_max_min_prot=protein_expression_mass_norm_df.loc[:,intersection_list]
    use_prot_diff=top_protein_mass_ratio(diff_max_min_prot,topnum,top_filename,pngname)

    diff_max_min_prot=protein_expression_mass_norm_df.loc[:,intersection_list]
    protein_expression_mass_norm_json=diff_max_min_prot.T.to_json()
    top_filename_json = os.path.join(analysis_result,strain,strain+'_intersection'+'_expdata.json')
    with open(top_filename_json, "w", encoding="utf-8") as f:
        f.write(protein_expression_mass_norm_json)
        
    seq_amino_composition_MW_norm_file='%s/%s/seq_amino_composition_g_g_norm_%s.csv'%(analysis_result,strain,strain)
    amino_composition_norm_onecell_outfile='%s/%s/amino_composition_g_g_norm_onecell_%s_intersection_protein.json'%(analysis_result,strain,strain)
    amino_acid_expression_mass_norm_json=amino_acid_expression_mass_norm(top_filename_json,seq_amino_composition_MW_norm_file,amino_composition_norm_onecell_outfile)

    # draw figure
    amino_composition_norm_onecell=json_load(amino_composition_norm_onecell_outfile)
    amino_composition_norm_onecell_df=pd.DataFrame()
    for key, value in amino_composition_norm_onecell.items():
        for key2 in value.keys():
            amino_composition_norm_onecell_df.loc[key,key2]=value[key2]['total']
            
    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_intersection_protein_%s.csv'%(analysis_result,strain,strain)
    amino_composition_norm_onecell_df.to_csv(amino_composition_norm_onecell_df_outfile, header=True, index=True) 

    pngname_boxplot='%s/%s_twenty_amino_condition_boxplot_intersection_protein.png'%(fig_file,strain)
    pngname_violinplot='%s/%s_twenty_amino_condition_violinplot_intersection_protein.png'%(fig_file,strain)
    draw_amino_mass_in_total_protein(amino_composition_norm_onecell_df,pngname_boxplot,pngname_violinplot)
    print(np.sum(np.mean(amino_composition_norm_onecell_df)))

    protein_list=diff_max_min_prot.columns.values.tolist()
    protein_mass_norm_stat_std=pd.DataFrame()
                
    for each_i in protein_list:
        data=diff_max_min_prot[each_i]
        protein_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        protein_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        protein_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)

    protein_mass_norm_stat_std.to_excel('%s/protein_mass_norm_stat_std_exp_onecell_intersection_protein.xlsx'%fig_file,index = True)
    print('protein_mass_norm_stat_std:')
    print(protein_mass_norm_stat_std)

    amino_composition_norm_onecell_df=amino_composition_norm_onecell_df.sort_index(axis = 1,ascending = True)
    amino_list=amino_composition_norm_onecell_df.columns.values.tolist()
    amino_mass_norm_stat_std=pd.DataFrame()
            
    for each_i in amino_list:
        data=amino_composition_norm_onecell_df[each_i]
        amino_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        amino_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        amino_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)
    
    amino_mass_norm_stat_std.to_excel('%s/amino_mass_norm_stat_std_amino_composition_proteome_intersection_protein.xlsx'%fig_file,index = True)
    print('amino_mass_norm_stat_std:')
    print(amino_mass_norm_stat_std)

    mw_stat=pd.DataFrame()

    data1=protein_mass_norm_stat_std['sd/mean']
    data2=amino_mass_norm_stat_std['sd/mean']
    mw_stat.loc['exp2amino','mean_in_exp']=np.mean(data1)
    mw_stat.loc['exp2amino','mean_in_AA']=np.mean(data2)

    if stats.levene(data1,data2).pvalue>0.5:#如果返回结果的p值远大于0.05，那么我们认为两总体具有方差齐性。
        print('Standard independent 2 sample test')
        mw_stat.loc['exp2amino','p_value']=stats.ttest_ind(data1,data2).pvalue#独立样本T检验
    else:
        print('Welch’s t-test')
        mw_stat.loc['exp2amino','p_value']=stats.ttest_ind(data1,data2,equal_var=False).pvalue
    mw_stat=mw_stat.sort_values(by='p_value',ascending = False)
    print(mw_stat)

'''2.7 蛋白组-相关性分析'''
def Proteome_Correlation_Analysis(topnum,fig_file,strain,analysis_result):
    # topnum=126

    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_by_condition_%s.csv'%(analysis_result,strain,strain)
    ori_pro=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)
    ori_pro=ori_pro.sort_index(axis = 1,ascending = True)
    ori_pro_mean=np.mean(ori_pro)

    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_top%sprotein_%s.csv'%(analysis_result,strain,topnum,strain)
    top_pro=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)
    top_pro=top_pro.sort_index(axis = 1,ascending = True)
    top_pro_mean=np.mean(top_pro) 

    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_intersection_protein_%s.csv'%(analysis_result,strain,strain)
    exp_pro=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)
    exp_pro=exp_pro.sort_index(axis = 1,ascending = True)
    exp_pro_mean=np.mean(exp_pro)


    pro_mean=pd.concat([ori_pro_mean,top_pro_mean,exp_pro_mean],axis=1)
    pro_mean.columns=['Total','Mass TOP','Mass and expression']
    pro_mean.name='Amino'
    pro_mean=pro_mean.sort_index(axis = 0,ascending = True)

    species_list=['Total','Mass TOP','Mass and expression']

    Total__Mass_and_expression__mean = pro_mean[['Total','Mass TOP','Mass and expression']]

    mw_stat=pd.DataFrame()
    for each_i in range(len(species_list)-1):
        for each_j in range(each_i,len(species_list)):
            if (each_i!=each_j):
                amino_vs=species_list[each_i]+'_'+species_list[each_j]
                mw_stat.loc[amino_vs,'corr']=stats.pearsonr(pro_mean[species_list[each_i]],pro_mean[species_list[each_j]])[0]
    mw_stat=mw_stat.sort_values(by='corr',ascending = False)
    print('mw_stat:')
    print(mw_stat)

    twoamino_use_data={}
    for eachamino in list(ori_pro.columns):
        twoamino_use_data[eachamino]=[]
        twoamino_use_data[eachamino].append(list(ori_pro[eachamino]))
        twoamino_use_data[eachamino].append(list(top_pro[eachamino]))
        twoamino_use_data[eachamino].append(list(exp_pro[eachamino]))

    data=[]
    for eachkey in twoamino_use_data.keys():
        data.append(twoamino_use_data[eachkey])
    labels = ['Total','Mass TOP','Mass and expression']
    colors = ['blue', 'red', 'green']
    ori_list=[1,2,3]
    fold=3.5
    positions=[ori_list,list(np.add(ori_list,1.5*fold)),list(np.add(ori_list,3*fold)),list(np.add(ori_list,4.5*fold)),list(np.add(ori_list,6*fold)),list(np.add(ori_list,7.5*fold)),list(np.add(ori_list,9*fold)),
            list(np.add(ori_list,10.5*fold)),list(np.add(ori_list,12*fold)),list(np.add(ori_list,13.5*fold)),list(np.add(ori_list,15*fold)),list(np.add(ori_list,16.5*fold)),list(np.add(ori_list,18*fold)),
            list(np.add(ori_list,19.5*fold)),list(np.add(ori_list,21*fold)),list(np.add(ori_list,22.5*fold)),list(np.add(ori_list,24*fold)),list(np.add(ori_list,25.5*fold)),list(np.add(ori_list,27*fold)),list(np.add(ori_list,28.5*fold))]

    x_position_list=[]
    for eachp in positions:
        x_position_list.append(np.mean(eachp))
    x_position=x_position_list#np.linspace(1,np.max(list(np.add(ori_list,28.5*fold))),20)
    x_position_fmt=list(ori_pro.columns)
    ylabel='Mass ratio (g/g total protein)'
    pngname='%s/twenty_amino_boxplot_comparison_topandexp%s.png'%(fig_file,strain)
    draw_box_vert(data,labels,colors,positions,x_position,x_position_fmt,ylabel,pngname,15,25)

    columns_name=['Total','Mass TOP','Mass and expression']
    legend=['Total','50%','highly variable and mass ratio']
    pngname='%s/Total__Mass_and_expression__mean_%s.png'%(fig_file,strain)
    draw_Total_MassAndExpression(Total__Mass_and_expression__mean,columns_name,legend,pngname)

'''2 蛋白组'''
def Proteome(topnum,fig_file,strain,analysis_result):
    # 2.1蛋白数目
    Proteome_Number_of_proteins(fig_file,strain,analysis_result)
    # Proteome_Number_of_proteins(fig_file,'BW25113_pro',analysis_result)
    
    # 2.2各蛋白氨基酸质量比
    Proteome_Amino_acid_mass_ratio_of_each_protein(fig_file,strain,analysis_result)
    # Proteome_Amino_acid_mass_ratio_of_each_protein(fig_file,'BW25113_pro',analysis_result)

    # 2.3总氨基酸的质量比
    Proteome_Mass_ratio_of_total_amino_acids(fig_file,strain,analysis_result)
    # Proteome_Mass_ratio_of_total_amino_acids(fig_file,'BW25113_pro',analysis_result)

    # 2.4质量总和超过50%的蛋白质量分布情况
    (use_prot,protein_mass_norm_stat_std) = Proteome_Protein_Distribution(topnum,fig_file,strain,analysis_result)
    # (use_prot,protein_mass_norm_stat_std) = Proteome_Protein_Distribution(topnum,fig_file,'BW25113_pro',analysis_result)

    # 2.5质量总和超过50%的蛋白的氨基酸总质量
    Proteome_Total_amino_acid_mass_of_protein(fig_file,strain,analysis_result,topnum,use_prot,protein_mass_norm_stat_std)
    # Proteome_Total_amino_acid_mass_of_protein(fig_file,'BW25113_pro',analysis_result,topnum,use_prot,protein_mass_norm_stat_std)
    
    # 2.6质量总和超过50%的蛋白种表达变化比较大的蛋白
    Proteome_High_protein_expression(topnum,fig_file,strain,analysis_result,use_prot)
    # Proteome_High_protein_expression(topnum,fig_file,'BW25113_pro',analysis_result,use_prot)

    # 2.7相关性分析
    Proteome_Correlation_Analysis(topnum,fig_file,strain,analysis_result)
    # Proteome_Correlation_Analysis(topnum,fig_file,'BW25113_pro',analysis_result)

'''3.1 转录组-各蛋白氨基酸质量比'''
def Transcriptome_Amino_acid_mass_ratio_of_each_protein(fig_file,strain,analysis_result):
    ecoprot_seq_amino_composition_g_g_norm=pd.read_csv('%s/%s/seq_amino_composition_g_g_norm_%s.csv'%(analysis_result,strain,strain)) 
    ecoprot_seq_amino_composition_g_g_norm=ecoprot_seq_amino_composition_g_g_norm.iloc[:,1:21]
    ecoprot_seq_amino_composition_g_g_norm=ecoprot_seq_amino_composition_g_g_norm.sort_index(axis = 1,ascending = True)

    pngname='%s/protein_amino_acid_composition_g_g_boxplot_%s.png'%(fig_file,strain)
    plt.figure(figsize=(25, 10)) 
    plt.tick_params(labelsize=20)
    sns.boxplot(data=ecoprot_seq_amino_composition_g_g_norm)
    plt.ylabel("Mass ratio (g/g protein)",fontsize=22)
    plt.xlabel("Amino acid", fontsize=22)
    plt.savefig(pngname,dpi =300,bbox_inches='tight')

    # plt.show()

    amino_list=ecoprot_seq_amino_composition_g_g_norm.columns.values.tolist()
    amino_mass_stat_std=pd.DataFrame()
                
    for each_i in amino_list:
        data=ecoprot_seq_amino_composition_g_g_norm[each_i]
        amino_mass_stat_std.loc[each_i,'mean']=np.mean(data)
        amino_mass_stat_std.loc[each_i,'sd']=np.std(data)
        amino_mass_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)
        
    print('amino_mass_stat_std:')
    print(amino_mass_stat_std)
    amino_mass_stat_std.to_excel('%s/amino_mass_stat_std.xlsx'%fig_file,index = True)

'''3.2 转录组-topnum个蛋白的氨基酸总质量'''
def Transcriptome_Total_amino_acid_mass_of_protein(topnum,fig_file,strain,analysis_result):
    protein_expression_mass_norm_outfile='%s/%s/%s_exp_onecell.json'%(analysis_result,strain,strain)
    protein_expression_mass_norm=json_load(protein_expression_mass_norm_outfile) 
    protein_expression_mass_norm_df=pd.DataFrame(protein_expression_mass_norm).T
    # protein_expression_mass_norm_df=protein_expression_mass_norm_df.drop(['BW25113.MD106.none.na_WT_T5490'])

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

    # topnum=1010
    top_filename = os.path.join(analysis_result,strain,strain+'_'+str(topnum)+'_expdata.csv')
    pngname='%s/protein_top%s_exp%s.png'%(fig_file,topnum,strain)
    use_prot=top_protein_mass_ratio(protein_expression_mass_norm_df,topnum,top_filename,pngname)

    print(topnum/protein_expression_mass_norm_df.shape[1])

    protein_list=use_prot.columns.values.tolist()
    protein_mass_norm_stat_std=pd.DataFrame()
                
    for each_i in protein_list:
        data=use_prot[each_i]
        protein_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        protein_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        protein_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)
    
    print('protein_mass_norm_stat_std:')
    print(protein_mass_norm_stat_std)
    protein_mass_norm_stat_std.to_excel('%s/protein_mass_norm_stat_std.xlsx'%fig_file,index = True)

    protein_expression_mass_norm_json=use_prot.T.to_json()
    top_filename_json = os.path.join(analysis_result,strain,strain+'_'+str(topnum)+'_expdata.json')
    with open(top_filename_json, "w", encoding="utf-8") as f:
        f.write(protein_expression_mass_norm_json)
        
    seq_amino_composition_MW_norm_file='%s/%s/seq_amino_composition_g_g_norm_%s.csv'%(analysis_result,strain,strain)
    amino_composition_norm_onecell_outfile='%s/%s/amino_composition_g_g_norm_onecell_%s_top%sprotein.json'%(analysis_result,strain,strain,topnum)
    amino_acid_expression_mass_norm_json=amino_acid_expression_mass_norm(top_filename_json,seq_amino_composition_MW_norm_file,amino_composition_norm_onecell_outfile)

    # draw figure
    amino_composition_norm_onecell=json_load(amino_composition_norm_onecell_outfile)
    amino_composition_norm_onecell_df=pd.DataFrame()
    for key, value in amino_composition_norm_onecell.items():
        for key2 in value.keys():
            amino_composition_norm_onecell_df.loc[key,key2]=value[key2]['total']
            
    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_top%sprotein_%s.csv'%(analysis_result,strain,topnum,strain)
    amino_composition_norm_onecell_df.to_csv(amino_composition_norm_onecell_df_outfile, header=True, index=True) 

    pngname_boxplot='%s/%s_twenty_amino_condition_boxplot_top%s_protein.png'%(fig_file,strain,topnum)
    pngname_violinplot='%s/%s_twenty_amino_condition_violinplot_top%s_protein.png'%(fig_file,strain,topnum)
    draw_amino_mass_in_total_protein(amino_composition_norm_onecell_df,pngname_boxplot,pngname_violinplot)

    print(np.sum(np.mean(amino_composition_norm_onecell_df)))

    amino_list=amino_composition_norm_onecell_df.columns.values.tolist()
    amino_mass_norm_stat_std=pd.DataFrame()
    amino_composition_norm_onecell_df=amino_composition_norm_onecell_df.sort_index(axis = 1,ascending = True)
                
    for each_i in amino_list:
        data=amino_composition_norm_onecell_df[each_i]
        amino_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        amino_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        amino_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)
        
    print('amino_mass_norm_stat_std:')
    print(amino_mass_norm_stat_std)
    amino_mass_norm_stat_std.to_excel('%s/amino_mass_norm_stat_std_amino_composition_proteome_topnum.xlsx'%fig_file,index = True)

    mw_stat=pd.DataFrame()

    data1=protein_mass_norm_stat_std['sd/mean']
    data2=amino_mass_norm_stat_std['sd/mean']
    mw_stat.loc['exp2amino','mean_in_exp']=np.mean(data1)
    mw_stat.loc['exp2amino','mean_in_AA']=np.mean(data2)

    if stats.levene(data1,data2).pvalue>0.5:#如果返回结果的p值远大于0.05，那么我们认为两总体具有方差齐性。
        print('Standard independent 2 sample test')
        mw_stat.loc['exp2amino','p_value']=stats.ttest_ind(data1,data2).pvalue#独立样本T检验
    else:
        print('Welch’s t-test')
        mw_stat.loc['exp2amino','p_value']=stats.ttest_ind(data1,data2,equal_var=False).pvalue
    mw_stat=mw_stat.sort_values(by='p_value',ascending = False)
    print(mw_stat)

    return use_prot
    
'''3.3 转录组-topnum个蛋白种表达变化比较大的蛋白'''
def Transcriptome_High_protein_expression(topnum,fig_file,strain,analysis_result,use_prot):
    protein_expression_mass_norm_outfile='%s/%s/%s_exp_onecell.json'%(analysis_result,strain,strain)
    protein_expression_mass_norm=json_load(protein_expression_mass_norm_outfile) 
    protein_expression_mass_norm_df=pd.DataFrame(protein_expression_mass_norm).T
    # protein_expression_mass_norm_df=protein_expression_mass_norm_df.drop(['BW25113.MD106.none.na_WT_T5490'])

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


    #变化范围大的蛋白
    # diff_max_min=np.max(protein_expression_mass_norm_df,axis = 0)-np.min(protein_expression_mass_norm_df,axis = 0)
    # diff_max_min_sort=diff_max_min.sort_values(ascending = False)
    diff_max_min_norm=(np.max(protein_expression_mass_norm_df,axis = 0)-np.min(protein_expression_mass_norm_df,axis = 0))/np.max(protein_expression_mass_norm_df,axis = 0)
    diff_max_min_sort=diff_max_min_norm.sort_values(ascending = False)

    max_gene_list=list(use_prot.columns)
    #diff_max_min_list=list(diff_max_min_sort[0:topnum].index)
    # diff_max_min_list=list(diff_max_min_sort[diff_max_min_sort>0.3].index)
    diff_max_min_list=list(diff_max_min_sort[diff_max_min_sort>0.5].index)

    intersection_list=list(set(max_gene_list).intersection(set(diff_max_min_list)))
    print(len(intersection_list))

    print(len(intersection_list)/protein_expression_mass_norm_df.shape[1])

    topnum=len(intersection_list)
    top_filename = os.path.join(analysis_result,strain,strain+'_'+str(topnum)+'_expdata.csv')
    pngname='%s/protein_top%s_exp%s.png'%(fig_file,topnum,strain)
    diff_max_min_prot=protein_expression_mass_norm_df.loc[:,intersection_list]
    use_prot_diff=top_protein_mass_ratio(diff_max_min_prot,topnum,top_filename,pngname)

    diff_max_min_prot=protein_expression_mass_norm_df.loc[:,intersection_list]
    protein_expression_mass_norm_json=diff_max_min_prot.T.to_json()
    top_filename_json = os.path.join(analysis_result,strain,strain+'_intersection'+'_expdata.json')
    with open(top_filename_json, "w", encoding="utf-8") as f:
        f.write(protein_expression_mass_norm_json)
        
    seq_amino_composition_MW_norm_file='%s/%s/seq_amino_composition_g_g_norm_%s.csv'%(analysis_result,strain,strain)
    amino_composition_norm_onecell_outfile='%s/%s/amino_composition_g_g_norm_onecell_%s_intersection_protein.json'%(analysis_result,strain,strain)
    amino_acid_expression_mass_norm_json=amino_acid_expression_mass_norm(top_filename_json,seq_amino_composition_MW_norm_file,amino_composition_norm_onecell_outfile)

    # draw figure
    amino_composition_norm_onecell=json_load(amino_composition_norm_onecell_outfile)
    amino_composition_norm_onecell_df=pd.DataFrame()
    for key, value in amino_composition_norm_onecell.items():
        for key2 in value.keys():
            amino_composition_norm_onecell_df.loc[key,key2]=value[key2]['total']
            
    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_intersection_protein_%s.csv'%(analysis_result,strain,strain)
    amino_composition_norm_onecell_df.to_csv(amino_composition_norm_onecell_df_outfile, header=True, index=True) 

    pngname_boxplot='%s/%s_twenty_amino_condition_boxplot_intersection_protein.png'%(fig_file,strain)
    pngname_violinplot='%s/%s_twenty_amino_condition_violinplot_intersection_protein.png'%(fig_file,strain)
    draw_amino_mass_in_total_protein(amino_composition_norm_onecell_df,pngname_boxplot,pngname_violinplot)

    print(np.sum(np.mean(amino_composition_norm_onecell_df)))

    protein_list=diff_max_min_prot.columns.values.tolist()
    protein_mass_norm_stat_std=pd.DataFrame()

    for each_i in protein_list:
        data=diff_max_min_prot[each_i]
        protein_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        protein_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        protein_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)
        
    print('protein_mass_norm_stat_std:')
    print(protein_mass_norm_stat_std)
    protein_mass_norm_stat_std.to_excel('%s/protein_mass_norm_stat_std_intersection_protein.xlsx'%fig_file,index = True)


    amino_list=amino_composition_norm_onecell_df.columns.values.tolist()
    amino_mass_norm_stat_std=pd.DataFrame()
    amino_composition_norm_onecell_df=amino_composition_norm_onecell_df.sort_index(axis = 1,ascending = True)
    
    for each_i in amino_list:
        data=amino_composition_norm_onecell_df[each_i]
        amino_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        amino_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        amino_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)
        
    print('amino_mass_norm_stat_std:')
    print(amino_mass_norm_stat_std)
    amino_mass_norm_stat_std.to_excel('%s/amino_mass_norm_stat_std_amino_composition_proteome_intersection_protein.xlsx'%fig_file,index = True)

    mw_stat=pd.DataFrame()

    data1=protein_mass_norm_stat_std['sd/mean']
    data2=amino_mass_norm_stat_std['sd/mean']
    mw_stat.loc['exp2amino','mean_in_exp']=np.mean(data1)
    mw_stat.loc['exp2amino','mean_in_AA']=np.mean(data2)

    if stats.levene(data1,data2).pvalue>0.5:#如果返回结果的p值远大于0.05，那么我们认为两总体具有方差齐性。
        print('Standard independent 2 sample test')
        mw_stat.loc['exp2amino','p_value']=stats.ttest_ind(data1,data2).pvalue#独立样本T检验
    else:
        print('Welch’s t-test')
        mw_stat.loc['exp2amino','p_value']=stats.ttest_ind(data1,data2,equal_var=False).pvalue
    mw_stat=mw_stat.sort_values(by='p_value',ascending = False)
    print('mw_stat:')
    print(mw_stat)

    return topnum

'''3.4 转录组-相关性分析'''
def Transcriptome_Correlation_Analysis(fig_file,strain,analysis_result,topnum):
    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_by_condition_%s.csv'%(analysis_result,strain,strain)
    ori_pro=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)
    ori_pro=ori_pro.sort_index(axis = 1,ascending = True)
    ori_pro_mean=np.mean(ori_pro)

    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_top%sprotein_%s.csv'%(analysis_result,strain,topnum,strain)
    top_pro=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)
    top_pro=top_pro.sort_index(axis = 1,ascending = True)
    top_pro_mean=np.mean(top_pro) 

    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_intersection_protein_%s.csv'%(analysis_result,strain,strain)
    exp_pro=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)
    exp_pro=exp_pro.sort_index(axis = 1,ascending = True)
    exp_pro_mean=np.mean(exp_pro)


    pro_mean=pd.concat([ori_pro_mean,top_pro_mean,exp_pro_mean],axis=1)
    pro_mean.columns=['Total','Mass TOP','Mass and expression']
    pro_mean.name='Amino'
    pro_mean=pro_mean.sort_index(axis = 0,ascending = True)

    from scipy import stats
    species_list=['Total','Mass TOP','Mass and expression']
    mw_stat=pd.DataFrame()
    for each_i in range(len(species_list)-1):
        for each_j in range(each_i,len(species_list)):
            if (each_i!=each_j):
                amino_vs=species_list[each_i]+'_'+species_list[each_j]
                mw_stat.loc[amino_vs,'corr']=stats.pearsonr(pro_mean[species_list[each_i]],pro_mean[species_list[each_j]])[0]
    mw_stat=mw_stat.sort_values(by='corr',ascending = False)
    print(mw_stat.head())

    twoamino_use_data={}
    for eachamino in list(ori_pro.columns):
        twoamino_use_data[eachamino]=[]
        twoamino_use_data[eachamino].append(list(ori_pro[eachamino]))
        twoamino_use_data[eachamino].append(list(top_pro[eachamino]))
        twoamino_use_data[eachamino].append(list(exp_pro[eachamino]))

    data=[]
    for eachkey in twoamino_use_data.keys():
        data.append(twoamino_use_data[eachkey])
    labels = ['Total','Mass TOP','Mass and expression']
    colors = ['blue', 'red', 'green']
    ori_list=[1,2,3]
    fold=3.5
    positions=[ori_list,list(np.add(ori_list,1.5*fold)),list(np.add(ori_list,3*fold)),list(np.add(ori_list,4.5*fold)),list(np.add(ori_list,6*fold)),list(np.add(ori_list,7.5*fold)),list(np.add(ori_list,9*fold)),
            list(np.add(ori_list,10.5*fold)),list(np.add(ori_list,12*fold)),list(np.add(ori_list,13.5*fold)),list(np.add(ori_list,15*fold)),list(np.add(ori_list,16.5*fold)),list(np.add(ori_list,18*fold)),
            list(np.add(ori_list,19.5*fold)),list(np.add(ori_list,21*fold)),list(np.add(ori_list,22.5*fold)),list(np.add(ori_list,24*fold)),list(np.add(ori_list,25.5*fold)),list(np.add(ori_list,27*fold)),list(np.add(ori_list,28.5*fold))]

    x_position_list=[]
    for eachp in positions:
        x_position_list.append(np.mean(eachp))
    x_position=x_position_list#np.linspace(1,np.max(list(np.add(ori_list,28.5*fold))),20)
    x_position_fmt=list(ori_pro.columns)
    ylabel='Mass ratio (g/g total protein)'
    pngname='%s/twenty_amino_boxplot_comparison_topandexp%s.png'%(fig_file,strain)
    draw_box_vert(data,labels,colors,positions,x_position,x_position_fmt,ylabel,pngname,15,25)

'''3 转录组'''
def Transcriptome(topnum,fig_file,strain,analysis_result):
    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_by_condition_%s.csv'%(analysis_result,strain,strain)
    amino_composition_norm_onecell_df=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)

    pngname_boxplot='%s/%s_twenty_amino_condition_boxplot.png'%(fig_file,strain)
    pngname_violinplot='%s/%s_twenty_amino_condition_violinplot.png'%(fig_file,strain)
    draw_amino_mass_in_total_protein(amino_composition_norm_onecell_df,pngname_boxplot,pngname_violinplot)

    print(strain)
    amino_mass_norm_stat_std=pd.DataFrame()
    amino_composition_norm_onecell_df=amino_composition_norm_onecell_df.sort_index(axis = 1,ascending = True)
    amino_list=amino_composition_norm_onecell_df.columns.values.tolist()    

    for each_i in amino_list:
        data=amino_composition_norm_onecell_df[each_i]
        amino_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        amino_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        amino_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)

    print('amino_mass_norm_stat_std:')
    print(amino_mass_norm_stat_std)
    amino_mass_norm_stat_std.to_excel('%s/amino_mass_norm_stat_std_amino_composition_proteome_by_condition.xlsx'%fig_file,index = True)

    # 3.1 转录组-各蛋白氨基酸质量比
    Transcriptome_Amino_acid_mass_ratio_of_each_protein(fig_file,strain,analysis_result)

    # 3.2 转录组-质量总和超过50%的蛋白的氨基酸总质量
    use_prot = Transcriptome_Total_amino_acid_mass_of_protein(topnum,fig_file,strain,analysis_result)

    # 3.3 转录组-质量总和超过50%的蛋白种表达变化比较大的蛋白
    Transcriptome_High_protein_expression(topnum,fig_file,strain,analysis_result,use_prot)

    # 3.4 转录组-相关性分析
    Transcriptome_Correlation_Analysis(fig_file,strain,analysis_result,topnum)

'''4.1 mass ratio in different protein'''
def Mass_ratio_in_different_protein(fig_file,strain_list,strainName,analysis_result):
    Species_Number = len(strain_list)
    createVar = locals()

    twoamino_use_data = {}
    twoamino_columns = []
    ori_list = []
    positions = []
    fold=4.5
    for i in range(0, Species_Number):
        seq_amino_composition_MW_norm_file='%s/%s/seq_amino_composition_g_g_norm_%s.csv'%(analysis_result,strain_list[i],strain_list[i])
        seq_amino_composition_MW_norm_df=pd.read_csv(seq_amino_composition_MW_norm_file,index_col=0)
        createVar['twoamino'+str(i+1)]=seq_amino_composition_MW_norm_df.sort_index(axis = 1,ascending = True)
        ori_list.append(i+1)
        if i==0:
            twoamino_columns = list(createVar['twoamino1'].columns)

        if i==(Species_Number-1):
            positions = [ori_list]
            num=1.5
            for j in range(1, len(twoamino_columns)):
                positions.append(list(np.add(ori_list,num*fold)))
                num+=1.5
            for eachamino in twoamino_columns:
                twoamino_use_data[eachamino]=[]
                for each in range(0,Species_Number):
                    twoamino_use_data[eachamino].append(list(createVar['twoamino'+str(each+1)][eachamino]))

    data=[]
    for eachkey in twoamino_use_data.keys():
        data.append(twoamino_use_data[eachkey])

    # labels = ["BW25113_pro", "BW25113", "DH1",'MG1655','W3110']
    labels = strainName
    colors = ['blue', 'red', 'green','darkorange','violet']

    x_position_list=[]
    for eachp in positions:
        x_position_list.append(np.mean(eachp))
    x_position=x_position_list#np.linspace(1,np.max(list(np.add(ori_list,28.5*fold))),20)
    x_position_fmt=twoamino_columns
    ylabel='Mass ratio (g/g protein)'
    pngname='%s/twenty_amino_boxplot_comparison_g_g.png'%(fig_file)

    draw_box_vert(data,labels,colors,positions,x_position,x_position_fmt,ylabel,pngname,15,25)

'''4.2 total mass ratio'''
def Total_mass_ratio(fig_file,strain_list,strainName,analysis_result):

    Species_Number = len(strain_list)
    createVar = locals()

    twoamino_use_data = {}
    twoamino_columns = []
    ori_list = []
    positions = []
    fold=4.5

    for i in range(0, Species_Number):
        amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_by_condition_%s.csv'%(analysis_result,strain_list[i],strain_list[i])
        amino_composition_norm_onecell_df=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)
        createVar['twoamino'+str(i+1)]=amino_composition_norm_onecell_df.sort_index(axis = 1,ascending = True)
        ori_list.append(i+1)

        if i==0:
            twoamino_columns = list(createVar['twoamino1'].columns)
        
        if i==(Species_Number-1):
            positions = [ori_list]
            num=1.5
            for j in range(1, len(twoamino_columns)):
                positions.append(list(np.add(ori_list,num*fold)))
                num+=1.5
            for eachamino in twoamino_columns:
                twoamino_use_data[eachamino]=[]
                for each in range(0,Species_Number):
                    twoamino_use_data[eachamino].append(list(createVar['twoamino'+str(each+1)][eachamino]))

    data=[]
    for eachkey in twoamino_use_data.keys():
        data.append(twoamino_use_data[eachkey])
    labels = strainName
    colors = ['blue', 'red', 'green','darkorange','violet']
    x_position_list=[]
    for eachp in positions:
        x_position_list.append(np.mean(eachp))
    x_position=x_position_list#np.linspace(1,np.max(list(np.add(ori_list,28.5*fold))),20)
    x_position_fmt=twoamino_columns
    ylabel='Mass ratio (g/g total protein)'
    pngname='%s/twenty_amino_boxplot_comparison.png'%(fig_file)
    draw_box_vert(data,labels,colors,positions,x_position,x_position_fmt,ylabel,pngname,15,25)

'''4.3 均值比较'''
def Mean_Value_Comparison(fig_file,strain_list,strainName,analysis_result):
    Species_Number = len(strain_list)
    createVar = locals()

    newindex = 0
    pro_mean = pd.DataFrame()
    pro = pd.DataFrame()
    amino_list=[]
    for i in range(0, Species_Number):
        amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_by_condition_%s.csv'%(analysis_result,strain_list[i],strain_list[i])
        amino_composition_norm_onecell_df=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)
        amino_list=amino_composition_norm_onecell_df.columns.values.tolist()

        for index, row in amino_composition_norm_onecell_df.iterrows():
            for j in range(0,len(amino_list)):
                pro.loc[newindex,'species']=strain_list[i]
                pro.loc[newindex,'value']=row[amino_list[j]]
                pro.loc[newindex,'amino_name']=amino_list[j]
                newindex+=1
        
        createVar['twoamino'+str(i+1)]=amino_composition_norm_onecell_df.sort_index(axis = 1,ascending = True)
        createVar['twoamino'+str(i+1)+'_mean']=np.mean(createVar['twoamino'+str(i+1)])
        
        if i==0:
            pro_mean = pd.DataFrame(createVar['twoamino'+str(i+1)+'_mean'])
        else:
            pro_mean = pd.concat([pro_mean,createVar['twoamino'+str(i+1)+'_mean']],axis=1)
    pro_mean.columns=strain_list
    pro_mean.name='Amino'
    pro_mean=pro_mean.sort_index(axis = 0,ascending = True)

    pngname='%s/protein_amino_acid_composition_mean.png'%(fig_file)
    plt.figure(figsize=(26, 12)) 
    plt.tick_params(labelsize=20)
    #pro_mean.plot.bar(figsize=(25, 10))
    ind = np.arange(20)                # the x locations for the groups
    plt.xlim(-.2, 20)
    width = 0.15

    if strain_list != strainName: # legend中为拉丁名物种的改为斜体
        plt.rcParams['font.style'] ='italic' # 'normal'#

    # colors = ['blue', 'red', 'green','darkorange','violet']

    num = 0
    for i in range(0, Species_Number):
        num = i*width
        # 计算标准误差
        std_err = np.std(pro_mean[strain_list[i]], ddof=1) / np.sqrt(len(pro_mean[strain_list[i]]))
        print('std_err:',std_err)
        plt.bar(ind+num,pro_mean[strain_list[i]],width,label = strainName[i])
        # plt.errorbar(ind+num,pro_mean[strain_list[i]], yerr=std_err, color='black')
        plt.errorbar(ind+num,pro_mean[strain_list[i]], yerr=std_err, fmt='.', color='black', capsize=5)

    plt.xticks(np.arange(20) + ((Species_Number-1)/2)*width, list(pro_mean.index),rotation=0)
    #plt.bar(data=pro_mean)
    #sns.barplot(data=pro_mean.T)
    plt.ylabel("Mass ratio (g/g total protein)",fontsize=20)
    plt.xlabel("Amino acid", fontsize=22)  # 我们设置横纵坐标的标题。
    plt.legend(loc="upper right", fontsize=22)
    plt.xticks(rotation=0)
    plt.rcParams['font.style'] ='normal' # 'normal'#
    plt.savefig(pngname,dpi =300,bbox_inches='tight')

    # 5*4的画布，每个画布里表示一种氨基酸，4个物种的全部数据
    plt.figure(figsize=(25, 20))
    fig, axs = plt.subplots(nrows=5, ncols=4, figsize=(25, 20))
    fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.3, hspace=0.4)

    amino_list = sorted(amino_list)
    print(amino_list)
    amino_acids = ['Alanine', 'Arginine', 'Asparagine', 'Asparate', 'Cysteine',
               'Glutamine', 'Glutamate', 'Glycine', 'Histidine', 'Isoleucine',
               'Leucine', 'Lysine', 'Methionine', 'Phenylalanine', 'Proline',
               'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine']

    index = 0
    for i in range(5):
        for j in range(4):
            ax = axs[i, j]
            sns.violinplot(x="species", y="value", data=pro[pro['amino_name']==amino_list[index]], ax=ax)
            # sns.boxplot(x="species", y="value", data=pro[pro['amino_name']==amino_list[index]], ax=ax)
            ax.set_title(amino_acids[index], fontdict={'fontsize': 16})
            ax.set_xticklabels(strainName) # 更改 x 轴刻度

            labels = ax.get_xticklabels()
            if strain_list != strainName: # legend中为拉丁名物种的改为斜体
                plt.setp(labels, fontstyle="italic")

            ax.set_ylim([0, 0.12])
            ax.set_xlabel("strains")
            ax.set_ylabel("Mass ratio (g/g total protein)")
            index +=1

    plt.savefig('%s/twenty_AAs_total_proteins_subplots_violinplot.png'%(fig_file),dpi =300,bbox_inches='tight')
    # plt.savefig('%s/twenty_AAs_total_proteins_subplots_boxplot.png'%(fig_file),dpi =300,bbox_inches='tight')

    species_list=strain_list
    mw_stat=pd.DataFrame()
    for each_i in range(len(species_list)-1):
        for each_j in range(each_i,len(species_list)):
            if (each_i!=each_j):
                amino_vs=species_list[each_i]+'_'+species_list[each_j]
                mw_stat.loc[amino_vs,'corr']=stats.pearsonr(pro_mean[species_list[each_i]],pro_mean[species_list[each_j]])[0]
    mw_stat=mw_stat.sort_values(by='corr',ascending = False)
    print(mw_stat)

'''4.4 STD和偏离比例'''
def STD_and_deviation_ratio(fig_file,strain,analysis_result):
    print(strain)
    amino_composition_norm_onecell_df_outfile='%s/%s/amino_composition_proteome_by_condition_%s.csv'%(analysis_result,strain,strain)
    amino_composition_norm_onecell_df=pd.read_csv(amino_composition_norm_onecell_df_outfile,index_col=0)
    twoamino=amino_composition_norm_onecell_df.sort_index(axis = 1,ascending = True)

    amino_list=twoamino.columns.values.tolist()
    amino_mass_norm_stat_std=pd.DataFrame()
                
    for each_i in amino_list:
        data=twoamino[each_i]
        amino_mass_norm_stat_std.loc[each_i,'mean']=np.mean(data)
        amino_mass_norm_stat_std.loc[each_i,'sd']=np.std(data)
        amino_mass_norm_stat_std.loc[each_i,'sd/mean']=np.std(data)/np.mean(data)
        
    print('amino_mass_norm_stat_std:')
    print(amino_mass_norm_stat_std)
    amino_mass_norm_stat_std.to_excel('%s/amino_mass_norm_stat_std_amino_composition_proteome_by_condition.xlsx'%fig_file,index = True)

    print('np.max(sd):',np.max(amino_mass_norm_stat_std['sd']))

'''4.5 2*2 canvases'''
def AAs_Mass_ratio_in_different_species(fig_file,strain_list,strainName,analysis_result):
    # 2*2的画布，每个画布里表示一个物种的全部数据
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(25, 10),dpi =300)
    fig.subplots_adjust(wspace=0.15, hspace=0.4)

    index = 0
    for i in range(2):
        for j in range(2):
            ax = axs[i, j]
            
            strain = strain_list[index]

            ecoprot_seq_amino_composition_g_g_norm=pd.read_csv('%s/%s/seq_amino_composition_g_g_norm_%s.csv'%(analysis_result,strain,strain)) 
            ecoprot_seq_amino_composition_g_g_norm=ecoprot_seq_amino_composition_g_g_norm.iloc[:,1:21]
            ecoprot_seq_amino_composition_g_g_norm=ecoprot_seq_amino_composition_g_g_norm.sort_index(axis = 1,ascending = True)
            
            sns.boxplot(data=ecoprot_seq_amino_composition_g_g_norm, ax=ax)

            ax.set_title(strainName[index], fontdict={'fontsize': 26})

            if strain_list != strainName: # legend中为拉丁名物种的改为斜体
                ax.set_title(strainName[index], fontdict={'fontsize': 26}, fontstyle="italic")

            ax.set_ylabel("Mass ratio (g/g protein)",fontsize=18)
            ax.set_xlabel("Amino acid", fontsize=18)
            index +=1
    # 5*4 canvases
    # Species_Number = len(strain_list)
    # melted = pd.DataFrame()
    # for i in range(0, Species_Number):
    #     print(strain_list[i])
    #     ecoprot_seq_amino_composition_g_g_norm=pd.read_csv('%s/%s/seq_amino_composition_g_g_norm_%s.csv'%(analysis_result,strain_list[i],strain_list[i])) 
    #     ecoprot_seq_amino_composition_g_g_norm=ecoprot_seq_amino_composition_g_g_norm.iloc[:,1:21]
    #     ecoprot_seq_amino_composition_g_g_norm=ecoprot_seq_amino_composition_g_g_norm.sort_index(axis = 1,ascending = True)

    #     amino_list=ecoprot_seq_amino_composition_g_g_norm.columns.values.tolist()

    #     melted_df = ecoprot_seq_amino_composition_g_g_norm.melt(var_name='Amino acid', value_name='Mass ratio (g/g protein)')
    #     melted_df['species'] = strain_list[i]
        
    #     melted = pd.concat([melted, melted_df], axis=0)
    #     # print(ecoprot_seq_amino_composition_g_g_norm)
    # melted
    # # 5*4 canvases, each representing an amino acid, with full data for 4 species
    # plt.figure(figsize=(25, 20))
    # fig, axs = plt.subplots(nrows=5, ncols=4, figsize=(25, 20))
    # fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.3, hspace=0.4)

    # amino_list = sorted(amino_list)
    # print(amino_list)
    # amino_acids = ['Alanine', 'Arginine', 'Asparagine', 'Asparate', 'Cysteine',
    #             'Glutamine', 'Glutamate', 'Glycine', 'Histidine', 'Isoleucine',
    #             'Leucine', 'Lysine', 'Methionine', 'Phenylalanine', 'Proline',
    #             'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine']

    # index = 0
    # for i in range(5):
    #     for j in range(4):
    #         ax = axs[i, j]
    #         sns.boxplot(x="species", y="Mass ratio (g/g protein)", data=melted[melted['Amino acid']==amino_list[index]], ax=ax)
    #         ax.set_title(amino_acids[index], fontdict={'fontsize': 16})
    #         ax.set_xticklabels(strainName) # 更改 x 轴刻度
    #         labels = ax.get_xticklabels()
    #         if strain_list != strainName: # legend中为拉丁名物种的改为斜体
    #             plt.setp(labels, fontstyle="italic")

    #         ax.set_ylim([0, 0.55])
    #         ax.set_xlabel("strains")
    #         ax.set_ylabel("Mass ratio (g/g protein)")
    #         index +=1
    # # Distribution of the mass ratio of twenty AAs per unit mass of different proteins in different species

    plt.savefig('%s/twenty_AAs_different_proteins_subplots_2_2_boxplot.png'%(fig_file),dpi =300,bbox_inches='tight')

'''4.6 2*2 canvases'''
def Total_mass_ratio_top(fig_file,strain_list,strainName,analysis_result):
    # 2*2 canvases, each representing 1 species, with all data for 20 amino acids
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(25, 15),dpi =300)
    fig.subplots_adjust(wspace=0.15, hspace=0.4)

    strain_index = 0
    # topnum = 50
    for i in range(2):
        for j in range(2):
            ax = axs[i, j]
            strain = strain_list[strain_index]

            protein_expression_mass_norm_outfile='%s/%s/%s_exp_onecell.json'%(analysis_result,strain,strain)
            protein_expression_mass_norm=json_load(protein_expression_mass_norm_outfile) 
            protein_expression_mass_norm_df=pd.DataFrame(protein_expression_mass_norm).T

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
                # protein_expression_mass_norm_df=protein_expression_mass_norm_df.drop(['MG1655.MD117.O2-starvation.na_WT_T5796','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6128','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6149','MG1655.MD091.none.na_WT_T8315','MG1655.MD018.none.na_WT_T5750','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6149','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6128','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6145','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6126','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6153','MG1655.MD018.Cfs;Mcn.na_WT_T0040','MG1655.MD097.O2-starvation.na_WT_T6125','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6155','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6151','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6171','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6130','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6132','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6134','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6136','MG1655.MD018.Mcn.na_WT_T0029','MG1655.MD097.O2-starvation.na_WT_T6214','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6167','MG1655.MD010.lactate-shift.na_WT_T0840','MG1655.MD003.none.na_WT_T1588','MG1655.MD003.none.na_WT_T1590','MG1655.MD003.none.na_WT_T1593','MG1655.MD097.O2-starvation.na_WT_T6184','MG1655.MD097.O2-starvation.na_WT_T6216','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6165','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6169','MG1655.MD097.O2-starvation;CORM-3.na_WT_T6143','MG1655.MD090.none.na_WT_T8309','MG1655.MD010.lactate-shift.na_WT_T0841','MG1655.MD091.none.na_WT_T8318','MG1655.MD116.none.na_WT_T1713','MG1655.MD027.none.na_WT_T1941','MG1655.MD027.none.na_WT_T1942','MG1655.MD027.none.na_WT_T1943','MG1655.MD097.O2-starvation.na_WT_T6178'])    
            elif re.search('BW25113',strain):
                #condition: BW25113.MD106.none.na_WT_T5490 
                protein_expression_mass_norm_df=protein_expression_mass_norm_df.drop(['BW25113.MD106.none.na_WT_T5490'])     
            
            df_mean=pd.DataFrame(np.mean(protein_expression_mass_norm_df))
            tmp=df_mean.T    
            tmp2 = tmp.rename(index={tmp.index[0]:'bj'}) 
            tmp3=tmp2.T
            tmp4=tmp3.sort_values(by='bj',ascending=False)
            # gene=np.array(tmp4.index[0:topnum])
            gene = np.array(tmp4.index)
            # use_prot=protein_expression_mass_norm_df.loc[:,gene]
            use_prot = protein_expression_mass_norm_df[gene]

            # # 计算每列的最大值、最小值和均值
            # max_values = use_prot.max(axis=0)
            # min_values = use_prot.min(axis=0)
            # mean_values = use_prot.mean(axis=0)
            # ax.set_title(strainName[strain_index], fontdict={'fontsize': 26})
            # if strain_list != strainName: # legend中为拉丁名物种的改为斜体
            #     ax.set_title(strainName[strain_index], fontdict={'fontsize': 26}, fontstyle="italic")

            # # ax.plot(max_values, 'ro-', label='Max')
            # ax.plot(mean_values, 'bo-', label='Mean')
            # # ax.plot(min_values, 'go-', label='Min')
            # ax.set_xticks([]) #

            # # 设置图例和轴标签
            # # ax.legend()
            # ax.set_xlabel("Protein", fontsize=18)
            # ax.set_ylabel("Mass ratio (g/g total protein)", fontsize=18)
            # strain_index +=1


            # sns.boxplot(data=use_prot, showfliers=False, ax=ax)
            sns.boxplot(data=use_prot, ax=ax)

            ax.set_title(strainName[strain_index], fontdict={'fontsize': 26})
            if strain_list != strainName: # legend中为拉丁名物种的改为斜体
                ax.set_title(strainName[strain_index], fontdict={'fontsize': 26}, fontstyle="italic")

            # # Rectangle
            # print(ax.get_ylim()[1] * 0.95)
            # # rect = plt.Rectangle((1, 0), 50, ax.get_ylim()[1], edgecolor='red', fill=False)
            # rect = plt.Rectangle((0, 0), 50, ax.get_ylim()[1] * 0.95, linewidth=2, edgecolor='r', fill=False, zorder=4)
            # ax.add_patch(rect)

            ax.set_xticks([]) #
            ax.set_xlabel("Protein", fontsize=18)
            ax.set_ylabel("Mass ratio (g/g total protein)", fontsize=18)
            strain_index +=1


    # plt.savefig('%s/top%s_mass_distribution_subplots_2_2_boxplot.png'%(fig_file,topnum),dpi =300,bbox_inches='tight')
    plt.savefig('%s/mass_distribution_subplots_2_2_showfliers_boxplot.png'%(fig_file),dpi =300,bbox_inches='tight')

'''4 物种比较'''
def Comparison_of_other_species(fig_file,strain_list,strainName,analysis_result):
    # 氨基酸组成
    Mass_ratio_in_different_protein(fig_file,strain_list,strainName,analysis_result)

    # total mass ratio
    Total_mass_ratio(fig_file,strain_list,strainName,analysis_result)

    # 均值比较
    Mean_Value_Comparison(fig_file,strain_list,strainName,analysis_result)

    # STD和偏离比例
    for strain in strain_list:
        STD_and_deviation_ratio(fig_file,strain,analysis_result)

    # Distribution of the mass ratio of twenty AAs per unit mass of different proteins in different species
    AAs_Mass_ratio_in_different_species(fig_file,strain_list,strainName,analysis_result)

    # The mass distribution of top 50 proteins in per unit mass of total protein.
    Total_mass_ratio_top(fig_file,strain_list,strainName,analysis_result)


def main():
    parser = argparse.ArgumentParser(description='draw figure in article!')
    parser.add_argument('-fun','--function', type=str, default='MW_Comparison')
    parser.add_argument('-ff','--fig_file', type=str, default='../analysis_result/article_figures_py')
    parser.add_argument('-strains','--strain_list', nargs='+', default=['Corynebacterium_RNA_seq','BW25113','Bacillus','Yeast_single_cell'])
    parser.add_argument('-strainName','--strainName', nargs='+', default=['C. glutamicum','B. subtilis','E. coli', 'S. cerevisiae'])
    parser.add_argument('-ar','--analysis_result', type=str, default='../analysis_result')
    parser.add_argument('-top','--topnum', type=int, default=126)

    args = parser.parse_args()
    
    function = args.function
    fig_file = args.fig_file
    strain_list = args.strain_list
    strainName = args.strainName
    analysis_result = args.analysis_result
    topnum = args.topnum

    a = datetime.now()
    print("start time: ",a)

    if not os.path.exists(fig_file):
        os.makedirs(fig_file)

    # Execute separately according to function function selection
    if function == 'MW_Comparison':
        # 1. Comparison of protein MW by species
        # python draw_figure_in_article.py -fun 'MW_Comparison' -ff '../analysis_result/article_figures' -strains 'Corynebacterium_RNA_seq' 'MG1655' 'Bacillus' 'Yeast_single_cell' -strainName 'C. glutamicum' 'E. coli' 'B. subtilis' 'S. cerevisiae' -ar '../analysis_result/initial_data'
        # strainName=['E. coli','B. subtilis','S. cerevisiae','C. glutamicum']
        print('----------------------------')
        print('strain_list: ',strain_list)
        print('strainName: ',strainName)
        print('----------------------------')
        strain_dir = '_'.join(strain_list)
        fig_file = os.path.join(fig_file,'MW_Comparison',strain_dir)
        isFileExists(fig_file)
        Comparison_of_protein_MW_by_species(strain_list,strainName,fig_file,analysis_result)
    
    elif function == 'Proteome':
        # 2.Proteome
        # python draw_figure_in_article.py -fun 'Proteome' -ff '../analysis_result/article_figures' -strains 'W3110' -ar '../analysis_result/initial_data' -top 910
        print('----------------------------')
        strain = strain_list[0]
        print('strain: ',strain)
        print('----------------------------')
        fig_file = os.path.join(fig_file,'Proteome',strain)
        isFileExists(fig_file)
        Proteome(topnum,fig_file,strain,analysis_result)

    elif function == 'Transcriptome':
        # 3.Transcriptome
        # python draw_figure_in_article.py -fun 'Transcriptome' -ff '../analysis_result/article_figures' -strains 'Bacillus' -ar '../analysis_result/initial_data' -top 200
        print('----------------------------')
        strain = strain_list[0]
        print('strain: ',strain)
        print('----------------------------')
        fig_file = os.path.join(fig_file,'Transcriptome',strain)
        isFileExists(fig_file)                
        Transcriptome(topnum,fig_file,strain,analysis_result)
        # Transcriptome(topnum,fig_file,'BW25113',analysis_result)
    
    elif function == 'Species_Comparison':
        # 4. Comparison of species a. Comparison of several other strains of E. coli b. Comparison of other species
        # a.python draw_figure_in_article.py -fun 'Species_Comparison' -ff '../analysis_result/article_figures' -strains 'BW25113' 'DH1' 'MG1655' 'W3110' -strainName 'BW25113' 'DH1' 'MG1655' 'W3110' -ar '../analysis_result/initial_data'
        # b.python draw_figure_in_article.py -fun 'Species_Comparison' -ff '../analysis_result/article_figures' -strains 'Corynebacterium_RNA_seq' 'MG1655' 'Bacillus' 'Yeast_single_cell' -strainName 'C. glutamicum' 'E. coli' 'B. subtilis' 'S. cerevisiae' -ar '../analysis_result/initial_data'
        print('----------------------------')
        print('strain_list: ',strain_list)
        print('strainName: ',strainName)
        print('----------------------------')
        # strain_list = ['Corynebacterium_RNA_seq','BW25113','Bacillus','Yeast_single_cell']
        # strainName = ['C. glutamicum','B. subtilis','E. coli', 'S. cerevisiae']

        strain_dir = '_'.join(strain_list)
        fig_file = os.path.join(fig_file,'Species_Comparison',strain_dir)
        isFileExists(fig_file)
        Comparison_of_other_species(fig_file,strain_list,strainName,analysis_result)

    b = datetime.now()
    print("end time: ",b)
    minutes = ((b-a).seconds)/60
    print("draw_figure_in_article.py,run time(minutes): ",minutes)

if __name__ == '__main__':
    main()