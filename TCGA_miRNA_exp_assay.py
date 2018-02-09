import os, time, atexit

t = time.time()
ori_path = os.getcwd()

#-----------------------------------------------------------------------------------------------
#parameter

source_path = '/home/tin/2014_TCGA_TEST'

#cancer_type_list = ['BLCA', 'BRCA', 'COAD', 'HNSC', 'KICH',
#                    'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC',
#                    'PRAD', 'READ', 'STAD', 'THCA', 'UCEC']

cancer_type_list = ['STAD']
Top_processing = 10


say = "It's finish !"
#-----------------------------------------------------------------------------------------------

def farewell():
    print('GOOD JOB')
    print(say)

spend_time = {}

import pycode.RaxLib as RaxLib
import pycode.RaxTCGALib as RaxTCGALib
import pycode.Raxpy3LibStatistic_2013_10_19 as RaxStat
import pycode.Raxpy3LibmultiCPU_2013_09_15 as multiCPU
import pycode.build_new_isomiR_file_NTpair as R1
import pycode.delete_isoform_less_then_1_rpm as R2
import pycode.TCGASTAD_top_3_mixarm_expression as R3
import pycode.TCGASTAD_top_3_twoarm_expression as R4
import pycode.delete_no_express_mir2 as R5
import pycode.build_mir_reverse_table as R6


for cancer_type in cancer_type_list:
    
    #cancer_type = 'STAD'
    
    tt = time.time()
    
    N_source_path_1 = source_path+'/'+cancer_type+'/Nm'+'/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/'
    T_source_path_1 = source_path+'/'+cancer_type+'/Tm'+'/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/'
    #-----------------------------------------------------------------------------------------------
    #build_new_isomiR_file_NTpair
    N_result_path_1 = R1.running(ori_path, cancer_type, N_source_path_1, multiCPU, Top_processing)
    T_result_path_1 = R1.running(ori_path, cancer_type, T_source_path_1, multiCPU, Top_processing)
    
    
    #N_result_path_1 = '/home/tin/Lib/2014_TCGASTAD_TEST/TCGA_miRNA_exp_assay_module_ready_to_use/1_TCGA_tumor_normal_Add_isomiR_type/READ/normal/BCGSC__IlluminaHiSeq_miRNASeq/Level_3'
    #T_result_path_1 = '/home/tin/Lib/2014_TCGASTAD_TEST/TCGA_miRNA_exp_assay_module_ready_to_use/1_TCGA_tumor_normal_Add_isomiR_type/READ/tumor/BCGSC__IlluminaHiSeq_miRNASeq/Level_3'
    
    #-----------------------------------------------------------------------------------------------
    #delete_isoform_less_then_1_rpm
    N_result_path_2 = R2.running(ori_path, cancer_type, N_result_path_1, multiCPU, Top_processing)
    T_result_path_2 = R2.running(ori_path, cancer_type, T_result_path_1, multiCPU, Top_processing)
    
    
    #N_result_path_2 = '/home/tin/Lib/TCGA_miRNA_exp_assay_module/20130910TCGA_tumor_normal_Add_isomiR_type_delete_less_then_1rpm/STAD/normal/BCGSC__IlluminaHiSeq_miRNASeq/Level_3'
    #T_result_path_2 = '/home/tin/Lib/TCGA_miRNA_exp_assay_module/20130910TCGA_tumor_normal_Add_isomiR_type_delete_less_then_1rpm/STAD/tumor/BCGSC__IlluminaHiSeq_miRNASeq/Level_3'
    #-----------------------------------------------------------------------------------------------
    #build dir for cancer type
    pacessing_path = ori_path+'/result/'+cancer_type
    if os.path.exists(pacessing_path) == False:
        os.makedirs(pacessing_path)
    else:
        pass
    os.chdir(pacessing_path)
    os.system('cp -a "'+ori_path+'/pycode" "'+ori_path+'/result/'+cancer_type+'/pycode"')
    
    
    #-----------------------------------------------------------------------------------------------
    #TCGASTAD_top_3_mixarm_expression
    R3.running(pacessing_path, N_result_path_2, T_result_path_2, multiCPU, 'TCGA'+cancer_type+'_NT_top3_mixarm_0to1.txt', 1, Top_processing)
    
    #-----------------------------------------------------------------------------------------------
    #TCGASTAD_top_3_twoarm_expression
    R4.running(pacessing_path, N_result_path_2, T_result_path_2, multiCPU, 'TCGA'+cancer_type+'_NT_top3_twoarm_0to1.txt', 1, Top_processing)
    
    #-----------------------------------------------------------------------------------------------
    #delete_no_express_mir2
    R5.running(pacessing_path, 'TCGA'+cancer_type, RaxLib)
    
    #-----------------------------------------------------------------------------------------------
    #TCGASTAD_top_3_twoarm_expression
    R4.running(pacessing_path, N_result_path_2, T_result_path_2, multiCPU, 'TCGA'+cancer_type+'_NT_top3_twoarm_0to0.txt', 0, Top_processing)
    
    #-----------------------------------------------------------------------------------------------
    #build_mir_reverse_table
    R6.running(pacessing_path, 'TCGA'+cancer_type, RaxLib, RaxStat)
    
    print(time.time()-tt)
    spend_time[cancer_type] = time.time()-tt
    
    

print(time.time()-t)
print(spend_time)
atexit.register(farewell)