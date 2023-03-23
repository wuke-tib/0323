import requests, re
import pandas as pd
import json

# 下载.keg 文件
def download_data(url, outfile):
    res = requests.get(url, stream=True)
    with open(outfile, "wb") as pypkg:
        for chunk in res.iter_content(chunk_size=1024):
            if chunk:
                pypkg.write(chunk)
    pypkg.close()

def json_write(path, dictionary_data):
    
    json_output = json.dumps(dictionary_data, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)

# C开头：kegg的pathway的ID所在行，D开头：属于它的kegg的所有的基因；
# A,B是kegg的分类，总共是6个大类，42个小类
def get_path_from_kegg(infile, gene_path,A_B_C_file): # pathID_def, classA_classB, classB_classC
    # pathID_def_df = pd.DataFrame()
    gene_path_df = pd.DataFrame()
    # classA_classB_df = pd.DataFrame()
    # classB_classC_df = pd.DataFrame()
    A_B_C_def = pd.DataFrame()
    
    log = ''
    A_name = ''
    B_name = ''
    indexA = ''
    indexB = ''
    indexC = ''
    
    with open(infile) as f_in:      
        while 1:           
            line = f_in.readline()
            if not line:
                break
            elif re.match('\#',line) or re.match('\+',line) or re.match('\!',line) or re.match('\%',line):
                pass
            else:
                
                if line.startswith('D') == False:
                    data_line = line.split(' ')
                    if len(data_line) > 1:  # 去除只有B的单独一行
                        # print (data_line)
                        log = data_line[0]
                        lastValue = (data_line[-1]).replace('\n','') # data_line数组最后一个元素带有换行符！

                        if data_line[0].startswith('A'):
                            indexA = data_line[0]
                            # pathID_def_df.loc[indexA, 'Pathway_ID'] = indexA
                            # 用空格连接第二个以及之后的元素
                            temp = data_line[1:-1]
                            temp.append(lastValue)
                            A_name = ' '.join(temp)
                            # pathID_def_df.loc[indexA, 'Def'] = A_name 

                        if(data_line[0] == 'B'):  # ['B', '', '09121', 'Transcription']
                            indexB = ''.join([log,data_line[2]])
                            # pathID_def_df.loc[indexB, 'Pathway_ID'] = indexB

                            temp = data_line[3:-1]
                            temp.append(lastValue)
                            B_name = ' '.join(temp)
                            # pathID_def_df.loc[indexB, 'Def'] = B_name 
                            
                            # classA_classB_df.loc[indexA+'_'+indexB, 'Source_A'] = indexA 
                            # classA_classB_df.loc[indexA+'_'+indexB, 'Target_B'] = indexB

                        if(data_line[0] == 'C'):  # ['C', '', '', '', '03040', 'Spliceosome']
                            indexC = ''.join([log,data_line[4]])
                            # pathID_def_df.loc[indexC, 'Pathway_ID'] = indexC

                            A_B_C_def.loc[indexC, 'ClassA'] = indexA
                            A_B_C_def.loc[indexC, 'ClassA_name'] = A_name
                            A_B_C_def.loc[indexC, 'ClassB'] = indexB
                            A_B_C_def.loc[indexC, 'ClassB_name'] = B_name
                            A_B_C_def.loc[indexC, 'ClassC'] = indexC

                            if '[' in data_line[-1]:  
                                # pathID_def_df.loc[indexC, '[PATH:]'] = lastValue
                                # 取下标5到-2的元素，并组合成字符串
                                # pathID_def_df.loc[indexC, 'Def'] = ' '.join(data_line[5:-1])
                                A_B_C_def.loc[indexC, 'ClassC_name'] = ' '.join(data_line[5:-1])
                            else:
                                temp = data_line[5:-1]
                                temp.append(lastValue)
                                # pathID_def_df.loc[indexC, 'Def'] = ' '.join(temp)
                                A_B_C_def.loc[indexC, 'ClassC_name'] = ' '.join(temp)

                            # classB_classC_df.loc[indexB+'_'+indexC, 'Source_B'] = indexB 
                            # classB_classC_df.loc[indexB+'_'+indexC, 'Target_C'] = indexC

                else:
                    dataD_Left = line.split('\t')[0]
                    data_line = dataD_Left.split(' ') # ['D', '', '', '', '', '', 'cg0576', 'rpoB;', 'DNA-DIRECTED RNA..'] 
                    lastValue = (data_line[-1]).replace('\n','')

                    if(data_line[0] == 'D') and len(data_line) > 1:  
                        entry = data_line[6]    # 取每一行D右边的gene
                        indexD = indexC + '_'+ entry
                        gene_path_df.loc[indexD, 'Pathway_ID'] = indexC
                        gene_path_df.loc[indexD, 'Gene_name'] = entry
                        gene_path_df.loc[indexD, 'Gene_Symbol'] = data_line[7].strip(';')
                        
                        temp = data_line[8:-1]
                        temp.append(lastValue)
                        gene_path_df.loc[indexD, '(GeneBank)Name'] = ' '.join(temp)
                        gene_path_df.loc[indexD, 'Type'] = "pathway_gene" 
               
    f_in.close()

    # pathID_def_df.to_csv(pathID_def, header=True, index=True, sep='\t', mode='w')
    gene_path_df.to_csv(gene_path, header=True, index=True, sep='\t', mode='w')
    # classA_classB_df.to_csv(classA_classB, header=True, index=True, sep='\t', mode='w')
    # classB_classC_df.to_csv(classB_classC, header=True, index=True, sep='\t', mode='w')
    A_B_C_def.to_csv(A_B_C_file, header=True, index=False, sep='\t', mode='w')

    return (A_B_C_def, gene_path_df)

# 整理出path和gene对应关系的json文件
def get_pathB_gene_json(a_list, A_B_C_def, gene_path_df,outfile): # 
    pathway_gene_kegg = {}
    for a in a_list:
        temp_df = A_B_C_def[A_B_C_def["ClassA"] == a]
        for index, col in temp_df.iterrows():
            B_path = col['ClassB']
            if B_path not in pathway_gene_kegg.keys():
                pathway_gene_kegg[B_path] = {}
                C_path = col['ClassC']
                pathway_gene_kegg[B_path]["B_pathway_name"] = col['ClassB_name']

                one_pathway = gene_path_df[gene_path_df['Pathway_ID'] == C_path]
                pathway_gene_kegg[B_path]["gene_id"] = list(one_pathway['Gene_name'])
            else:
                one_pathway = gene_path_df[gene_path_df['Pathway_ID'] == C_path]
                pathway_gene_kegg[B_path]["gene_id"] = list(set(pathway_gene_kegg[B_path]["gene_id"] + (list(one_pathway['Gene_name'])))) 
                # print(B_path)
    # print(pathway_gene_kegg)
    json_write(outfile, pathway_gene_kegg)  
    return pathway_gene_kegg


