import os, pycode.RaxLib, pickle

ori_path = os.getcwd()

#-----------------------------------------------------------------------------------------------
#parameter
source_path = '/home/tin/Lib/20miRBase/Rax_TCGA/miRNAseq/STAD/Nm/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/'
#Raxpy3CPU = Raxpy3LibmultiCPU_2013_09_15

#-----------------------------------------------------------------------------------------------

def running(ori_path = '', cancer_type = '', source_path = '', Raxpy3CPU = 'Raxpy3LibmultiCPU_2013_09_15', Top_processing = 14):
    import os
    
    if '/Nm/' in source_path:
        tissue_type = 'normal'
    elif '/Tm/' in source_path:
        tissue_type = 'tumor'
    target_path = ori_path+'/1_TCGA_tumor_normal_Add_isomiR_type/'+cancer_type+'/'+tissue_type+'/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/'
    if os.path.exists(target_path) == False:
        os.makedirs(target_path)
    else:
        pass
    
    #-----------------------------------------------------------------------------------------------
    #open BatchID table
    os.chdir(ori_path)
    Batch_table = pycode.RaxLib.tabfilewithtitle()
    Batch_table.open('aliquotReport.txt', 'Aliquot ID')
    #f = open('aliquotReport.txt.pickletable', 'rb')
    #Batch_table = pickle.load(f)
    #f.close()
    os.chdir(ori_path)
    
    #-----------------------------------------------------------------------------------------------
    #take STADisomiR_file_list
    os.chdir(source_path)
    STADfile_list = os.listdir(source_path)
    
    if len([x for x in STADfile_list if 'hg19' in x]) != 0:
        STADisomiR_file_list = [x for x in STADfile_list if 'isoform' in x and 'hg19' in x]
    else:
        STADisomiR_file_list = [x for x in STADfile_list if 'isoform' in x]
    
    #-----------------------------------------------------------------------------------------------
    #check patient ID in BatchID table
    STADisomiR_file_list = [x for x in STADisomiR_file_list if x.split('.')[0] in Batch_table.main_key_list]
    
    
    processing_code = '''

import os
ori_path = os.getcwd()
os.chdir(ori_path+'/pycode/')
import pycode.RaxLib, pycode.RaxTCGALib, copy

#-----------------------------------------------------------------------------------------------
#parameter
RaxLibbasic = pycode.RaxLib
RaxLibdataTCGA = pycode.RaxTCGALib
source_path = "'''+source_path+'''"
targer_path = "'''+target_path+'''"

outputfile_Name = ''
#-----------------------------------------------------------------------------------------------

dataTCGA = RaxLibdataTCGA.dataTCGA()
dataTCGA.openpercursor()
dataTCGA.openhsamiRNA()
dataTCGA.openmatchtable()
dataTCGA.openlocationtable()

os.chdir(source_path)
STADfile_list = os.listdir(source_path)
STADisomiR_file_list = [x for x in STADfile_list if x.find('isoform') != -1]#@
for file_name in STADisomiR_file_list:
    
    os.chdir(source_path)
    
    #file_name = 'bcgsc.ca__IlluminaHiSeq_miRNASeq__TCGA-BR-6453-01A-11R-1802-13__isoform_quantification.txt'
    outputfile_Name = '.'.join(file_name.split('.')[0:-1])+'_1'+'.txt'
    
    miRNA_info_file = '.'.join(file_name.split('.')[0:-1])+'_1_miRNA_info'+'.txt'
    
    total_isomiR_table = RaxLibbasic.tabfilewithtitle()
    total_isomiR_table.open(file_name, 'isoform_coords')
    
    #-----------------------------------------------------------------------------------------------
    #seperate the key accroding to miRNA_ID
    
    total_miRNA_ID_list = total_isomiR_table.title_box[total_isomiR_table.title_dicX['miRNA_ID']-1]
    miRNA_ID_list = list(set(total_miRNA_ID_list))
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
            
            rowread = total_isomiR_table.readrow('Y', key)
            sub_isomiR_table.append('Y', key, rowread)
            #delete = total_isomiR_table.delete(key)
            #sub_isomiR_table.append('Y', key, delete)
        
        sub_table_list.append(sub_isomiR_table)
        del sub_isomiR_table
        del title_dicX
    #print len(sub_table_list)
    
    #-----------------------------------------------------------------------------------------------
    #extend data with ismiR info
    
    sub_isomiR_table_list = []
    #print '#-----------------------------------------------------------------------------------------------'
    
    miRNA_info = []
    percursor_hairpin_list = []
    for sub_table_obj in sub_table_list:
        
        sub_table = copy.copy(sub_table_obj)
        
        #sub_table.printtable()
        
        isomiR_table = dataTCGA.isomiR()
        isomiR_table.process1(sub_table, dataTCGA)
        isomiR_table.process5(isomiR_table.isomiR_tabfilewithtitle, dataTCGA)
        isomiR_table.process2(isomiR_table.isomiR_tabfilewithtitle, dataTCGA)
        isomiR_table.process6(isomiR_table.isomiR_tabfilewithtitle, dataTCGA)
        isomiR_table.process7(isomiR_table.isomiR_tabfilewithtitle, dataTCGA)
        
        sub_isomiR_table_list.append(isomiR_table)
        
        if isomiR_table.percursor_hairpin == 'na':
            miRNA_info.append(isomiR_table.miRNA)
        else:
            percursor_hairpin_list.append(isomiR_table.percursor_hairpin)
        
        #print '----------------------------NEXT miRNA'
        del sub_table
        del isomiR_table
    
    #-----------------------------------------------------------------------------------------------
    #rebuild table
    
    new_table = RaxLibbasic.tabfilewithtitle()
    new_table.build(sub_isomiR_table_list[0].isomiR_tabfilewithtitle.title_dicX.copy(), 'isoform_coords')
    
    for isomiR_obj in sub_isomiR_table_list:
        
        main_key_list = list(isomiR_obj.isomiR_tabfilewithtitle.title_dicY.keys()).copy()
        main_key_list.sort()
        for key in main_key_list:
            delete = isomiR_obj.isomiR_tabfilewithtitle.delete(key)
            new_table.append('Y', key, delete)
            
        del main_key_list
    
    #-----------------------------------------------------------------------------------------------
    #report new_table
    
    os.chdir(targer_path)
    
    report_W = new_table.report()
    f = open(outputfile_Name, 'w')
    f.write(report_W)
    f.close()
    
    f = open(miRNA_info_file, 'w')
    
    f.write('file_name is '+outputfile_Name+'\n')
    f.write('total_miRNA quantity is '+str(len(sub_isomiR_table_list))+'\n')
    f.write('total_precursor_hairpin quantity is '+str(len(percursor_hairpin_list))+'\n')
    
    #build total_miRNA_list
    total_miRNA_list = []
    
    for sub_table_obj in sub_isomiR_table_list:
        sub_table = copy.copy(sub_table_obj)
        total_miRNA_list.append(sub_table.miRNA)
        
        del sub_table
        del sub_table_obj
    
    f.write('------------------------------------------------------------------------'+'\n')
    f.write('total_miRNA_list is :'+'\n')
    f.write('\t'.join(total_miRNA_list)+'\n')
    f.write(str(len(total_miRNA_list))+'\n')
    f.write('------------------------------------------------------------------------'+'\n')
    f.write('not_available_miRNA_list is :'+'\n')
    f.write('\t'.join(miRNA_info)+'\n')
    f.write(str(len(miRNA_info)))
    f.close()
    
    #print '----------------------------Next file'
    
    
    '''
    
    
    os.chdir(ori_path)
    Raxpy3CPU.multiprocessing(STADisomiR_file_list, processing_code, 'STADisomiR_file_list', Top_processing, 'Y')
    return target_path




