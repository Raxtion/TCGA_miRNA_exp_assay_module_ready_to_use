import os, RaxLibbasic_2013_05_16, copy

ori_path = os.getcwd()

os.chdir('/home/tin/Lib/Lib/hsagff3_location/modify_2013-10-01_hsagff3/2_modify_table/')
#open table
table = RaxLibbasic_2013_05_16.tabfilewithtitle()
table.open('hsagff3modification.txt', 'order')

os.chdir('/home/tin/Lib/Lib/hsamiRNA/')
#open FastA
FastA = RaxLibbasic_2013_05_16.openFastAfile()
FastA.open('hsaV20.fa')

FastA_key_list = FastA.FastA_dictionary.keys()

#set Name_hsaV2_list
Name_hsaV2_list = []

accession_number_list = table.title_box[table.title_dicX['accession_number']-1]
Name_list = table.title_box[table.title_dicX['Name']-1]
i = 1
for value in accession_number_list:
    if 'MIMA' in value:
        
        E = [x for x in FastA_key_list if value in x]
        Name_hsaV2_list.append("".join(E).split(' ')[0])
    
    else:
        
        precursor = Name_list[i-1]
        Name_hsaV2_list.append(precursor)
    
    i = i + 1

#add Name_hsaV2_list
table.insert('X', 'Name_hsaV2', Name_hsaV2_list, 8)


os.chdir(ori_path)
#report table
f = open('hsagff3modifyhsaV2.txt', 'w')
f.write(table.report())
f.close()