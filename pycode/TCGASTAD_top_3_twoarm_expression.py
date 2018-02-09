import os, time

ori_path = os.getcwd()

def running(ori_path = '', N_source_path = '', T_source_path = '', Raxpy3CPU = 'Raxpy3LibmultiCPU_2013_09_15', Final_file_Name = 'TCGASTAD_NT_top3_twoarm_0to001.txt', ZeroChange = 1, Top_processing = 14):
    import os, time
    
    def mulitCUP_for_TCGA_data(ori_path, source_path, output_file_Name, Raxpy3CPU):
    
        #-----------------------------------------------------------------------------------------------
        #parameter
        #source_path = '/home/tin/Lib/20miRBase/Rax_TCGA/20130910TCGA_tumor_normal_Add_isomiR_type_delete_less_then_1rpm/STAD/normal/BCGSC__IlluminaHiSeq_miRNASeq/Level_3'
        #Raxpy3CPU = Raxpy3LibmultiCPU_2013_09_15
        ori_path = ori_path
        source_path = source_path
        output_file_Name = output_file_Name
        Raxpy3CPU = Raxpy3CPU
        #-----------------------------------------------------------------------------------------------
        os.chdir(source_path)
        type_list = [x for x in os.listdir() if '_cleaned.txt' in x]
        
        #outputfile_Name = __file__.split('.')[0]+'_Temp.txt'
        #source_path = "'''+source_path+'''"
        processing_code = '''

import os
ori_path = os.getcwd()
import pycode.RaxLib, pycode.RaxTCGALib, copy

#-----------------------------------------------------------------------------------------------
#parameter
RaxLibbasic = pycode.RaxLib
RaxLibdataTCGA = pycode.RaxTCGALib

source_path = "'''+source_path+'''"


#os.chdir(source_path)
type_list = [x for x in os.listdir() if '_cleaned.txt' in x]#@

test_type_list = type_list

outputfile_Name = __file__.split('.')[0]+'_Temp.txt'
#-----------------------------------------------------------------------------------------------

dataTCGA = RaxLibdataTCGA.dataTCGA()

#-----------------------------------------------------------------------------------------------
#input matchtable
os.chdir('/home/tin/Lib/Lib/hsa_match_table/')
TabTmt = RaxLibbasic.tabfilewithtitle()
TabTmt.open('hsaV20_miRNA_pair_table1.txt', 'precursor')

#-----------------------------------------------------------------------------------------------
#build noly one arm mir list
only_5p_mir_list = [x for x in TabTmt.main_key_list if TabTmt.read('miRNA_5p', x) != 'na' and TabTmt.read('miRNA_3p', x) == 'na']
only_3p_mir_list = [x for x in TabTmt.main_key_list if TabTmt.read('miRNA_5p', x) == 'na' and TabTmt.read('miRNA_3p', x) != 'na']

#-----------------------------------------------------------------------------------------------
#statistic item
file_mir_list_dic = {}
#file_mir_marker_dic = {}
file_mir_marker_5p_dic = {}
file_mir_marker_3p_dic = {}

#-----------------------------------------------------------------------------------------------
#loading dir
for file_name in test_type_list:
    
    os.chdir(source_path)
    #file_name = 'bcgsc.ca__IlluminaHiSeq_miRNASeq__TCGA-B7-5816-01A-21R-1602-13__isoform_quantification_1_cleaned.txt'
    
    total_isomiR_table = RaxLibbasic.tabfilewithtitle()
    total_isomiR_table.open(file_name, 'isoform_coords')
    
    #-----------------------------------------------------------------------------------------------
    #seperate the key accroding to miRNA_ID
    
    total_miRNA_ID_list = total_isomiR_table.title_box[total_isomiR_table.title_dicX['miRNA_ID']-1]
    miRNA_ID_list = list({}.fromkeys(total_miRNA_ID_list).keys())
    miRNA_ID_list.sort()
    main_key_list = total_isomiR_table.main_key_list
    
    
    key_miRNA_ID_list = []
    for miRNA in miRNA_ID_list:
        
        E = []
        for key in main_key_list:
            result = total_isomiR_table.read('miRNA_ID', key)
            
            if result == miRNA:
                E.append(key)
            
        key_miRNA_ID_list.append(E)
    
    #-----------------------------------------------------------------------------------------------
    #build table_list
    
    sub_table_list = []
    for key_list in key_miRNA_ID_list:
        
        title_dicX = total_isomiR_table.title_dicX.copy()
        
        sub_isomiR_table = RaxLibbasic.tabfilewithtitle()
        sub_isomiR_table.build(title_dicX, 'isoform_coords')
        
        for key in key_list:
            
            delete = total_isomiR_table.delete(key)
            sub_isomiR_table.append('Y', key, delete)
        
        sub_table_list.append(sub_isomiR_table)
        del sub_isomiR_table
        del title_dicX
    
    #print(len(sub_table_list))
    
    #-----------------------------------------------------------------------------------------------
    #set statistic item
    mir_list_for_file = []
    #marker_dic = {}
    marker_5p_dic = {}
    marker_3p_dic = {}
    
    #-----------------------------------------------------------------------------------------------
    
    sub_isomiR_table_list = []
    for sub_table_obj in sub_table_list:
        
        sub_table = copy.copy(sub_table_obj)
        isomiR_table = dataTCGA.isomiR()
        isomiR_table.process8(sub_table, dataTCGA)
        
        #print(file_name, isomiR_table.miRNA)
        
        #-----------------------------------------------------------------------------------------------
        #defined 3p 5p key_group
        main_key_list = isomiR_table.isomiR_tabfilewithtitle.main_key_list
        
        if isomiR_table.miRNA in only_5p_mir_list:
            CM_key_5p_E = [x for x in main_key_list if isomiR_table.isomiR_tabfilewithtitle.read('arm', x) == '5p']
            CM_key_3p_E = [x for x in main_key_list if isomiR_table.isomiR_tabfilewithtitle.read('arm', x) == 'undefined']
        elif isomiR_table.miRNA in only_3p_mir_list:
            CM_key_5p_E = [x for x in main_key_list if isomiR_table.isomiR_tabfilewithtitle.read('arm', x) == 'undefined']
            CM_key_3p_E = [x for x in main_key_list if isomiR_table.isomiR_tabfilewithtitle.read('arm', x) == '3p']
        else:
            CM_key_5p_E = [x for x in main_key_list if isomiR_table.isomiR_tabfilewithtitle.read('arm', x) == '5p']
            CM_key_3p_E = [x for x in main_key_list if isomiR_table.isomiR_tabfilewithtitle.read('arm', x) == '3p']
        
        #-----------------------------------------------------------------------------------------------
        #find top 3 key
        M_5p = {}
        for key in CM_key_5p_E:
            rpm = isomiR_table.isomiR_tabfilewithtitle.read('reads_per_million_miRNA_mapped', key)
            M_5p[key] = rpm
        
        M_3p = {}
        for key in CM_key_3p_E:
            rpm = isomiR_table.isomiR_tabfilewithtitle.read('reads_per_million_miRNA_mapped', key)
            M_3p[key] = rpm
        
        top_3_in_5p = RaxLibbasic.dicKeysortVal(M_5p, '->')[:3]
        top_3_in_3p = RaxLibbasic.dicKeysortVal(M_3p, '->')[:3]
        
        #-----------------------------------------------------------------------------------------------
        #build top 3 express
        marker_5p = 0.0
        marker_3p = 0.0
        
        i = 1
        for key in top_3_in_5p:
            rpm = isomiR_table.isomiR_tabfilewithtitle.read('reads_per_million_miRNA_mapped', key)
            marker_5p = marker_5p + rpm
            i = i + 1
        
        i = 1
        for key in top_3_in_3p:
            rpm = isomiR_table.isomiR_tabfilewithtitle.read('reads_per_million_miRNA_mapped', key)
            marker_3p = marker_3p + rpm
            i = i + 1
        
        #marker = '//'.join([str(marker_5p), str(marker_3p)])
        #marker = marker_3p + marker_5p
        
        #-----------------------------------------------------------------------------------------------
        #fix in marker_dic, mir_list_for_file
        #marker_dic[isomiR_table.miRNA] = marker
        marker_5p_dic[isomiR_table.miRNA] = marker_5p
        marker_3p_dic[isomiR_table.miRNA] = marker_3p
        mir_list_for_file.append(isomiR_table.miRNA)
        
    #-----------------------------------------------------------------------------------------------
    #fix in file_mir_list_dic, file_mir_marker_dic
    file_mir_list_dic[file_name] = mir_list_for_file
    #file_mir_marker_dic[file_name] = marker_dic
    file_mir_marker_5p_dic[file_name] = marker_5p_dic
    file_mir_marker_3p_dic[file_name] = marker_3p_dic


##-----------------------------------------------------------------------------------------------
##input matchtable
#os.chdir('/home/tin/Lib/Lib/hsa_match_table/')
#TabTmt = RaxLibbasic.tabfilewithtitle()
#TabTmt.open('hsaV20_miRNA_pair_table1.txt', 'precursor')

#-----------------------------------------------------------------------------------------------
#delete 12 miRNA exception
miRNA_exception_list = ['hsa-mir-3180-1', 'hsa-mir-3180-2', 'hsa-mir-3180-3', 'hsa-mir-522', 'hsa-mir-523', 'hsa-mir-527',
                        'hsa-mir-550a-1', 'hsa-mir-550a-2', 'hsa-mir-570', 'hsa-mir-598', 'hsa-mir-378d-1', 'hsa-mir-378d-2']
for exct in miRNA_exception_list:
    TabTmt.delete(exct)

#-----------------------------------------------------------------------------------------------
#buildmiRNA_name
mir_list = TabTmt.main_key_list
new_mir_list = [x+'-5p' for x in TabTmt.main_key_list]+[x+'-3p' for x in TabTmt.main_key_list]
#new_mir_list = TabTmt.main_key_list
new_mir_list.sort()

#-----------------------------------------------------------------------------------------------
#build a order table
TabT = RaxLibbasic.newOrderfile(['order', 'mir'], [[x] for x in new_mir_list],)
#TabT.printtable()

#-----------------------------------------------------------------------------------------------
#fix file_mir_marker in TabT
for file_name in test_type_list:
    
    _M_5p = file_mir_marker_5p_dic[file_name]
    _M_3p = file_mir_marker_3p_dic[file_name]
    _E = []
    for mir in new_mir_list:
        if '-5p' in mir:
            mir_5p = mir[:-3]
            try:
                rpm = _M_5p[mir_5p]
            except:
                rpm = 0.0
        elif '-3p' in mir:
            mir_3p = mir[:-3]
            try:
                rpm = _M_3p[mir_3p]
            except:
                rpm = 0.0
        else:
            rpm = 0.0
        
        _E.append(rpm)
    
    TabT.append('X', file_name, _E)
    
    #break

#-----------------------------------------------------------------------------------------------
#TabT modify 0to001 and delete 'order'

X = TabT.title_dicX.keys()
for X_key in X:
    j = 1
    for value in TabT.title_box[TabT.title_dicX[X_key]-1]:
        
        if value == 0.0:
            TabT.title_box[TabT.title_dicX[X_key]-1][j-1] = '''+str(ZeroChange)+'''
            
        j = j + 1

TabT.mainkeyChange('mir')
TabT.delete('order')
TabT.XYreverse()

#-----------------------------------------------------------------------------------------------
#add 'NT' and 'tissue type' annotation
NT_list = []
for key in TabT.main_key_list:
    #print(key)
    if 'normal' in source_path.split('2_TCGA_tumor_normal_Add_isomiR_type_delete_less_then_1rpm/')[-1].split('/')[1]:
        NT_list.append('N')
    elif 'tumor' in source_path.split('2_TCGA_tumor_normal_Add_isomiR_type_delete_less_then_1rpm/')[-1].split('/')[1]:
        NT_list.append('T')
    else:
        Q=Q #quit

TabT.insert('X', 'NTpair', NT_list, 1)

tissue_list = []
for key in TabT.main_key_list:
    #print(key)
    tissue_list.append(source_path.split('2_TCGA_tumor_normal_Add_isomiR_type_delete_less_then_1rpm/')[-1].split('/')[0])

TabT.insert('X', 'tissue', tissue_list, 1)


#-----------------------------------------------------------------------------------------------
#write out TabT
os.chdir(ori_path)
f = open(outputfile_Name, 'w')
f.write(TabT.report())
f.close()


        '''
        
        os.chdir(ori_path)
        Raxpy3CPU.multiprocessing(type_list, processing_code, 'type_list', Top_processing, 'N')
        
        
        #-----------------------------------------------------------------------------------------------
        #wait all processing finish. It is for the processing take more then 12 H.
        i = 1
        j = 1
        while i < 2:
            if len([x for x in os.listdir(ori_path) if '_Temp.txt' in x]) == len([x for x in os.listdir(ori_path) if 'py_end_' in x and '.py' in x]):
                break
            else:
                pass
            
            time.sleep(10)
            print('{0:.2f}'.format(j/6), '(minute)')
            j = j + 1
        
        #-----------------------------------------------------------------------------------------------
        #combined Tabfile
        
        file_list = [x for x in os.listdir(ori_path) if '_Temp.txt' in x]
        
        Title = open(file_list[0], 'r').read().split('\n')[0]
        
        oupwoud_list = []
        for file_name in file_list:
            W_ = '\n'.join(open(file_name, 'r').read().split('\n')[1:-1])
            oupwoud_list.append(W_)
        
        f = open(output_file_Name, 'w')
        f.write(Title+'\n')
        f.write('\n'.join(oupwoud_list)+'\n')
        f.close()
        
        #-----------------------------------------------------------------------------------------------
        #clean Temp
        os.chdir(ori_path)
        for file_name in [x for x in os.listdir(ori_path) if '_Temp.txt' in x]:
            os.system('rm '+file_name)
        for file_name in [x for x in os.listdir(ori_path) if 'py_' in x and '.py' in x]:
            os.system('rm '+file_name)
        for file_name in [x for x in os.listdir(ori_path) if 'job_part_' in x and '.txt' in x]:
            os.system('rm '+file_name)
        for file_name in [x for x in os.listdir(ori_path) if 'py_end_' in x and '.py' in x]:
            os.system('rm '+file_name)
    
    
    #-----------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------
    #main()
    
    #Raxpy3CPU = Raxpy3LibmultiCPU_2013_09_15
    #source_path = '/home/tin/Lib/20miRBase/Rax_TCGA/20130910TCGA_tumor_normal_Add_isomiR_type_delete_less_then_1rpm/STAD/normal/BCGSC__IlluminaHiSeq_miRNASeq/Level_3'
    source_path = N_source_path
    output_file_Name = 'TabTouT_for_N_twoarm.txt'
    mulitCUP_for_TCGA_data(ori_path, source_path, output_file_Name, Raxpy3CPU)
    
    #source_path = '/home/tin/Lib/20miRBase/Rax_TCGA/20130910TCGA_tumor_normal_Add_isomiR_type_delete_less_then_1rpm/STAD/tumor/BCGSC__IlluminaHiSeq_miRNASeq/Level_3'
    source_path = T_source_path
    output_file_Name = 'TabTouT_for_T_twoarm.txt'
    mulitCUP_for_TCGA_data(ori_path, source_path, output_file_Name, Raxpy3CPU)
    
    
    #-----------------------------------------------------------------------------------------------
    #Next processing combined Tabfile
    
    file_list = [x for x in os.listdir(ori_path) if 'TabTouT_for_' in x and 'twoarm.txt' in x]
    
    Title = open(file_list[0], 'r').read().split('\n')[0]
    
    oupwoud_list = []
    for file_name in file_list:
        W_ = '\n'.join(open(file_name, 'r').read().split('\n')[1:-1])
        oupwoud_list.append(W_)
    
    f = open(Final_file_Name, 'w')
    f.write(Title+'\n')
    f.write('\n'.join(oupwoud_list)+'\n')
    f.close()



