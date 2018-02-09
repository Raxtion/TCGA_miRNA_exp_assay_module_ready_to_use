import os, RaxLibbasic_2013_05_16

ori_path = os.getcwd()

os.chdir('/home/tin/Lib/Lib/')

rowdata = open('modify_2013-10-01_hsa.gff3.txt', 'r').read()
Line_list = rowdata.split('\n')
del Line_list[-1]

E = []
i = 1
for Line in Line_list:
    E.append('\t'.join([str(i)+'_', Line]))
    
    i = i + 1

result =  '\n'.join(E)

title = '\t'.join(['order', 'chromosome', 'A', 'type', 'start', 'end', 'B', 'strand', 'C', 'ID'])

os.chdir(ori_path)

f = open('modifyhsagff3_build_table.txt', 'w')
f.write(title+'\n')
f.write(result+'\n')
f.write(str(len(E)))
f.close()
