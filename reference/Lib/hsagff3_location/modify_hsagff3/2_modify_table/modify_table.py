import os, RaxLibbasic_2013_05_16, copy

ori_path = os.getcwd()

os.chdir('/home/tin/Lib/Lib/hsagff3_location/modify_hsagff3/1_generate_tbale')
table = RaxLibbasic_2013_05_16.tabfilewithtitle()
table.open('modifyhsagff3_build_table.txt', 'order')

#delete '.' column
table.delete('A')
table.delete('B')
table.delete('C')
old_ID_list = [x.replace('\r', "") for x in table.delete('ID')]

#copy a new table
new_table = copy.copy(table)

#set new column
ID_list = []
accession_number_list = []
Name_list = []
derives_from = []

for line in old_ID_list:
    
    ID = line.split(';')[0].split('=')[-1]
    ac = line.split(';')[1].split('=')[-1]
    Na = line.split(';')[2].split('=')[-1]
    
    ID_list.append(ID)
    accession_number_list.append(ac)
    Name_list.append(Na)
    
    
    if line.split(';')[-1].find('Name') != -1:
        derives_from.append('na')
    else:
        de = line.split(';')[-1].split('=')[-1]
        derives_from.append(de)


#add four column
new_table.append('X', 'ID', ID_list)
new_table.append('X', 'accession_number', accession_number_list)
new_table.append('X', 'Name', Name_list)
new_table.append('X', 'derives_from', derives_from)

os.chdir(ori_path)
#report table
f = open('hsagff3modification.txt', 'w')
f.write(new_table.report())
f.close()

