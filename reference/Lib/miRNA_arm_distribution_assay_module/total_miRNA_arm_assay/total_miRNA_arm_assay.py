

import os
import file_seperate_by_miRNA
import file_seperate_by_arm
import list_collection_except_data
import arm_length_percentage_count
import seq_distribution


#-----------------------------------------------------------------------------------------------
# 2 parameter
assay_file_in_FasatA = 'matureV20.fa'
start_dir = '/home/tin/Lib/20miRBase/Rax_sequence_with_matureV20/Separated_by_miRNA/total_miRNA_arm_assay'+'/'
seperation_condiction = 0

#-----------------------------------------------------------------------------------------------
#Processing 1
file_name_1 = assay_file_in_FasatA
start_dir = start_dir
seperation_Address = start_dir
statistic_Address = start_dir

os.chdir(start_dir)
P1 = file_seperate_by_miRNA.processing(file_name_1, seperation_Address, statistic_Address)
seperation_Address_matureV1FastA = P1


#-----------------------------------------------------------------------------------------------
#Processing 2
file_name_1 = file_name_1
start_dir = start_dir
seperation_condiction = seperation_condiction
collection_Address = start_dir+file_name_1.split('.')[0]+'_'+'-'+str(seperation_condiction)+'/'+'1_make_-'+str(seperation_condiction)+'/'
statistic_Address = collection_Address
output_file_name = file_name_1.split('.')[0]+'-'+str(seperation_condiction)+'.'+file_name_1.split('.')[-1]

os.makedirs(collection_Address)
P2 = list_collection_except_data.processing(P1, collection_Address, seperation_condiction, output_file_name)


#-----------------------------------------------------------------------------------------------
#Processing 3
file_name_2 = output_file_name
start_dir = start_dir
seperation_Address = start_dir+file_name_1.split('.')[0]+'_'+'-'+str(seperation_condiction)+'/'+'2_S_miRNA'+'/'
statistic_Address = seperation_Address

os.makedirs(seperation_Address)
P3 = file_seperate_by_miRNA.processing(file_name_2, seperation_Address, statistic_Address)
seperation_Address_matureV1FastA = P3


#-----------------------------------------------------------------------------------------------
#Processing 4
file_name_2 = file_name_2
start_dir = start_dir
seperation_Address = start_dir+file_name_1.split('.')[0]+'_'+'-'+str(seperation_condiction)+'/'+'3_S_arm'+'/'
statistic_Address = seperation_Address

os.makedirs(seperation_Address)
P4 = file_seperate_by_arm.processing(P3, seperation_Address, statistic_Address)


#-----------------------------------------------------------------------------------------------
#Processing 5
file_name_2 = file_name_2
start_dir = start_dir
seperation_Address = start_dir+file_name_1.split('.')[0]+'_'+'-'+str(seperation_condiction)+'/'+'4_l_p'+'/'
statistic_Address = seperation_Address

os.makedirs(seperation_Address)
P5 = arm_length_percentage_count.processing(P3, seperation_Address)


#-----------------------------------------------------------------------------------------------
#Processing 6
file_name_1 = file_name_1
file_name_2 = file_name_2
start_dir = start_dir
statistic_Address = start_dir+file_name_1.split('.')[0]+'_'+'-'+str(seperation_condiction)+'/'+'5_len_distribution'+'/'
os.makedirs(statistic_Address)

os.chdir(start_dir)
P6 = seq_distribution.fileprocessing(file_name_1 , None, seperation_Address)
f1_display = []+ P6[0]
f1_display_count = P6[1]
f1_len_range_list = []+ P6[2]


os.chdir(start_dir+file_name_1.split('.')[0]+'_'+'-'+str(seperation_condiction)+'/'+'1_make_-'+str(seperation_condiction)+'/')
P7 = seq_distribution.fileprocessing(file_name_2 , None, seperation_Address)
f2_display = []+P7[0]
f2_display_count = P7[1]

P8 = seq_distribution.listprocessing(P4, '5p', None, seperation_Address)
arm_5p_display = []+P8[0]
arm_5p_display_count = P8[1]


P9 = seq_distribution.listprocessing(P4, '3p', None, seperation_Address)
arm_3p_display = []+P9[0]
arm_3p_display_count = P9[1]


P10 = seq_distribution.listprocessing(P4, '3p5p', None, seperation_Address)
arm_3p5p_display = []+P10[0]
arm_3p5p_display_count = P10[1]


f1_len_range_list.insert(0,'seq_len')
f1_len_range_list.insert(-1,'total')
f1_display.insert(0,file_name_1)
f1_display.insert(-1,f1_display_count)
f2_display.insert(0,file_name_2)
f2_display.insert(-1,f2_display_count)
arm_5p_display.insert(0,"".join(file_name_2.split('.')[0])+'-'+'5p')
arm_5p_display.insert(-1,arm_5p_display_count)
arm_3p_display.insert(0,"".join(file_name_2.split('.')[0])+'-'+'3p')
arm_3p_display.insert(-1,arm_3p_display_count)
arm_3p5p_display.insert(0,"".join(file_name_2.split('.')[0])+'-'+'3p5p')
arm_3p5p_display.insert(-1,arm_3p5p_display_count)


E = []
W= ''
range_list = []
for x in f1_len_range_list:
    range_list.append("".join(x.split('mer')[0]))

E.append(range_list)
E.append(f1_display)
E.append(f2_display)
E.append(arm_5p_display)
E.append(arm_3p_display)
E.append(arm_3p5p_display)

def listXYreverse(XY = 'list'):
    
    i = 1
    j = 1
    E = []
    YX = []
    while i < len(XY[0]) + 1:
        
        j = 1
        E = []
        while j < len(XY) + 1:
            
            E.append(XY[j-1][i-1])
            j = j + 1
            
        YX.append(E)
      
        i = i + 1

    return YX



R = listXYreverse(E)
GG = []
for x in R:
    GG.append('\t'.join(x))


os.chdir(statistic_Address)
f = open('len_distribution.txt','w')
f.write('\n'.join(GG))
f.close()

os.chdir(start_dir)





