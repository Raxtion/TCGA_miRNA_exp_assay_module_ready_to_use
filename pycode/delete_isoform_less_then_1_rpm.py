import os

ori_path = os.getcwd()

#-----------------------------------------------------------------------------------------------
#parameter
source_path = '/home/tin/Lib/20miRBase/Rax_TCGA/20130910TCGA_tumor_normal_Add_isomiR_type/STAD/normal/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/'
#Raxpy3CPU = Raxpy3LibmultiCPU_2013_09_15
#-----------------------------------------------------------------------------------------------

def running(ori_path = '', cancer_type = '', source_path = '', Raxpy3CPU = 'Raxpy3LibmultiCPU_2013_09_15', Top_processing = 14):
    import os
    
    target_path = ori_path+'/2_TCGA_tumor_normal_Add_isomiR_type_delete_less_then_1rpm/'+cancer_type+source_path.split(cancer_type)[1]
    if os.path.exists(target_path) == False:
        os.makedirs(target_path)
    else:
        pass
    
    os.chdir(source_path)
    STADfile_list = os.listdir(source_path)
    STADisomiR_file_list = [x for x in STADfile_list if '_miRNA_info.txt' not in x]
    
    processing_code = '''

import os
ori_path = os.getcwd()
os.chdir(ori_path+'/pycode/')
import pycode.RaxLib

ori_path = os.getcwd()

#-----------------------------------------------------------------------------------------------
#parameter
RaxLibbasic = pycode.RaxLib
source_path = "'''+source_path+'''"
targer_path = "'''+target_path+'''"

#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
#loading dir
file_list = os.listdir(source_path)
isomiR_file_list = [x for x in file_list if '_miRNA_info.txt' not in x]#@

for input_file in isomiR_file_list:
    #-----------------------------------------------------------------------------------------------
    #open file
    os.chdir(source_path)
    total_isomiR_table = open(input_file, 'r').read()
    Title_line = total_isomiR_table.split('\\n')[0]
    Title_list = Title_line.split('\\t')
    Y_line_list = total_isomiR_table.split('\\n')[1:-1]
    
    outputfile_Name = input_file.split('.')[0]+'_cleaned.txt'
    
    #-----------------------------------------------------------------------------------------------
    #delete less then 1 rpm
    os.chdir(targer_path)
    f = open(outputfile_Name, 'w')
    f.write(Title_line+'\\n')
    
    countread = 0.0
    i = 0
    for Y_line in Y_line_list:
        
        read = float(Y_line.split('\\t')[Title_list.index('read_count')])
        read_rpm = float(Y_line.split('\\t')[Title_list.index('reads_per_million_miRNA_mapped')])
        
        if read_rpm < 1.0:
            pass
        else:
            f.write(Y_line+'\\n')
            i = i + 1
        
        countread = countread + read
        print('build output file', i, end = '\\r')
        
    print('')
    #-----------------------------------------------------------------------------------------------
    #output total_isomiR_table
    #os.chdir(ori_path)
    #f = open(outputfile_Name, 'w')
    #f.write(total_isomiR_table.report())
    f.write(str(i))
    f.close()
    
    if countread < 1000000:
        os.system('rm '+outputfile_Name)
    else:
        pass

'''


    
    os.chdir(ori_path)
    Raxpy3CPU.multiprocessing(STADisomiR_file_list, processing_code, 'isomiR_file_list', Top_processing, 'Y')
    return target_path



