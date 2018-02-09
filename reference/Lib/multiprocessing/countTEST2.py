import os, Raxpy3Libbasic_2013_06_10, Raxpy3LibdataTCGA_2013_06_11, copy

ori_path = os.getcwd()

#-----------------------------------------------------------------------------------------------
#parameter
RaxLibbasic = Raxpy3Libbasic_2013_06_10
RaxLibdataTCGA = Raxpy3LibdataTCGA_2013_06_11
source_path = '/home/tin/Lib/Rax_TCGA/tumor_normal/LUAD/tumor/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/'
targer_path = ori_path
outputfile_Name = ''
#-----------------------------------------------------------------------------------------------

dataTCGA = RaxLibdataTCGA.dataTCGA()
dataTCGA.openpercursor()
dataTCGA.openhsamiRNA()
dataTCGA.openmatchtable()
dataTCGA.openlocationtable()

os.chdir(source_path)
STADfile_list = os.listdir(source_path)
STADisomiR_file_list = [x for x in STADfile_list if x.find('isoform') != -1]
for file_name in STADisomiR_file_list:
    file_name = 'bcgsc.ca__IlluminaHiSeq_miRNASeq__TCGA-05-4384-01A-01T-1754-13__isoform_quantification.txt'
    
    os.chdir(source_path)
    
    #file_name = 'bcgsc.ca__IlluminaHiSeq_miRNASeq__TCGA-B7-5816-01A-21R-1602-13__isoform_quantification.txt'
    outputfile_Name = '.'.join(file_name.split('.')[0:-1])+'_1'+'.txt'
    
    miRNA_info_file = '.'.join(file_name.split('.')[0:-1])+'_1_miRNA_info'+'.txt'
    
    total_isomiR_table = RaxLibbasic.tabfilewithtitle()
    total_isomiR_table.open(file_name, 'isoform_coords')
    total_isomiR_table.delete('barcode')
    
    #-----------------------------------------------------------------------------------------------
    #seperate the key accroding to miRNA_ID
    print('part1')
    total_miRNA_ID_list = total_isomiR_table.title_box[total_isomiR_table.title_dicX['miRNA_ID']-1]
    miRNA_ID_list = {}.fromkeys(total_miRNA_ID_list).keys()
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
    #print('part2')
    #-----------------------------------------------------------------------------------------------
    #import threading, Queue
    #class ClientThread(threading.Thread):
    #    def run(self):
    #        while True:
    #        item = Pool.get()
    #        if item != None:
    #            
    #            sub_isomiR_table = RaxLibbasic.tabfilewithtitle()
    #            sub_isomiR_table.build(title_dicX, 'isoform_coords')
    #            
    #            for key in key_list:
    #                
    #                delete = total_isomiR_table.delete(key)
    #                sub_isomiR_table.append('Y', key, delete)
    #            
    #            sub_table_list.append(sub_isomiR_table)
    #            del sub_isomiR_table
    #            
    #            Pool.task_done()
    #
    #
    #Pool = Queue.Queue(0)
    #
    #for x in range(2):
    #    ClientThread.start()
    #
    #for key_list in key_miRNA_ID_list:
    #    Pool.put(key_list)
    #
    #Pool.join()
    #print(Pool)
    
    #-----------------------------------------------------------------------------------------------
    #build table_list
    print('part2')
    sub_table_list = []
    title_dicX = total_isomiR_table.title_dicX
    import Text_1
    R = Text_1.draadje(key_miRNA_ID_list)
    print(len(R))
    
    
    
    
    #for key_list in key_miRNA_ID_list:
    #    
        #sub_isomiR_table = RaxLibbasic.tabfilewithtitle()
        #sub_isomiR_table.build(title_dicX, 'isoform_coords')
        #
        #for key in key_list:
        #    
        #    delete = total_isomiR_table.delete(key)
        #    sub_isomiR_table.append('Y', key, delete)
        #
        #sub_table_list.append(sub_isomiR_table)
        #del sub_isomiR_table
    #print len(sub_table_list)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #-----------------------------------------------------------------------------------------------
    #isomiR select: delete the isomiRNA read_count<miRNA_total_count*(10%)
    print('part3')
    for sub_table in sub_table_list:
        miRNA_total_count = 0.0
        for key in sub_table.main_key_list:
            read_count = sub_table.read('read_count', key)
            miRNA_total_count = miRNA_total_count + read_count
        
        
        for key in sub_table.main_key_list:
            read_count = sub_table.read('read_count', key)
            if read_count < miRNA_total_count*(0.1):
                sub_table.delete(key)
    
    #-----------------------------------------------------------------------------------------------
    #rebuild table
    print('part4')
    new_table = RaxLibbasic.tabfilewithtitle()
    new_table.build(total_isomiR_table.title_dicX, 'isoform_coords')
    
    for sub_table in sub_table_list:
        
        main_key_list = copy.copy(sub_table.main_key_list)
        for key in main_key_list:
            delete = sub_table.delete(key)
            new_table.append('Y', key, delete)
            
        del main_key_list
    
    #-----------------------------------------------------------------------------------------------
    #report new_table
    print('part5')
    os.chdir(targer_path)
    
    report_W = new_table.report()
    f = open(outputfile_Name, 'w')
    f.write(report_W)
    f.close()
    quit()
    
    #print '----------------------------Next file'