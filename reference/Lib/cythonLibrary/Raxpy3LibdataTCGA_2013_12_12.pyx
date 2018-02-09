'''
history =  2013.05.13, 2013_05_19, 2013_05_24, 2013_06_11, 2013_06_20, 2013_06_21, 2013.07.18,
           2013.07.26, 2013.08.18
'''
#-----------------------------------------------------------------------------------------------
import os, RaxLib
#-----------------------------------------------------------------------------------------------
RaxLib = RaxLib

class dataTCGA:
    
    def __init__(self):
        
        #file list
        self.isomiR_class_list = []
        
        #genome file
        self.hg_obj_name = ['chr1.fa', 'chr2.fa', 'chr3.fa', 'chr4.fa', 'chr5.fa',
                       'chr6.fa', 'chr7.fa', 'chr8.fa', 'chr9.fa', 'chr10.fa',
                       'chr11.fa', 'chr12.fa', 'chr13.fa', 'chr14.fa', 'chr15.fa',
                       'chr16.fa', 'chr17.fa', 'chr18.fa', 'chr19.fa', 'chr20.fa',
                       'chr21.fa', 'chr22.fa', 'chrM.fa', 'chrX.fa', 'chrY.fa']
        self.hg_obj_index = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                        '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                        '21', '22', 'MT', 'X', 'Y']
        self.hg_obj_list = [[]]*25
        
        #percursor file
        self.percursor = 'obj'
        
        #hsamiRNA file
        self.mature_miRNA = 'obj'
        
        #hsamatchtable
        self.matchtable = 'obj'
        
        #hsalocationtable
        self.locationtable = 'obj'
    
    #-----------------------------------------------------------------------------------------------
    
    def openHumangenome(self):
        
        ori_path = os.getcwd()
        os.chdir('/home/tin/Lib/Lib/hg19')
        chromosome = 'obj'
        
        for genomefile in self.hg_obj_name:
            
            print(str(genomefile)+' is opening !! (' + str(len([x for x in self.hg_obj_list if x != []])) + '/25)')
            chromosome = RaxLib.openFastAfile()
            chromosome.open(genomefile)
            self.hg_obj_list[self.hg_obj_name.index(genomefile)] = chromosome
        
        os.chdir(ori_path)
        
        return
    
    #-----------------------------------------------------------------------------------------------
    
    def openhsamiRNA(self):
        
        ori_path = os.getcwd()
        #take miRNA from mature
        os.chdir('/home/tin/Lib/Lib/hsamiRNA/')
        mature_miRNA = RaxLib.openFastAfile()
        mature_miRNA.open('hsaV20.fa')
        self.mature_miRNA = mature_miRNA
        
        os.chdir(ori_path)
        
        return
    
    #-----------------------------------------------------------------------------------------------
    
    def openpercursor(self):
        
        ori_path = os.getcwd()
        #take percursor from hairpin
        os.chdir('/home/tin/Lib/Lib/miRNApercursor')
        percursor = RaxLib.openFastAfile()
        percursor.open('hairpin.fa')
        self.percursor = percursor
        
        os.chdir(ori_path)
        
        return
    
    #-----------------------------------------------------------------------------------------------
    
    def openmatchtable(self):
        
        ori_path = os.getcwd()
        #take hsaprecursor and hsamiRNA matchtable
        os.chdir('/home/tin/Lib/Lib/hsa_match_table/')
        matchtable = RaxLib.tabfilewithtitle()
        matchtable.open('group_table1.txt', 'order')
        self.matchtable = matchtable
        
        os.chdir(ori_path)
        
        return
    
    #-----------------------------------------------------------------------------------------------
    
    def openlocationtable(self):
        
        ori_path = os.getcwd()
        #take hsaprecursor and hsamiRNA locationtable
        os.chdir('/home/tin/Lib/Lib/hsagff3_location/')
        gff3table = RaxLib.tabfilewithtitle()
        gff3table.open('hsagff3modifyhsaV20.txt', 'order')
        self.locationtable = gff3table
        
        os.chdir(ori_path)
        
        return
    
    #-----------------------------------------------------------------------------------------------
    
    class isomiR:
        
        def __init__(self):
            
            #miRNA
            self.miRNA = ""
            self.situation = ""
            self.mature_3p_start = 0.0
            self.mature_3p_end = 0.0
            self.mature_5p_start = 0.0
            self.mature_5p_end = 0.0
            
            #percursor
            self.percursor_hairpin = ""
            self.percursor_sequence = ""
            self.percursor_h_sequence_ali = ""
            self.percursor_genome_left = 0
            self.percursor_genome_right = 0
            
            #genome
            self.genome_sequence = ""
            self.genome_chromosome = ""
            self.genome_version = ""
            self.genome_strand = ""
            self.genome_left = 0
            self.genome_right = 0
            
            #table
            self.isomiR_tabfilewithtitle = 'obj'
        
        #-----------------------------------------------------------------------------------------------
        
        def process1(self, tabfilewithtitle = 'table obj', genomedataTCGA = 'dataTCGA obj'): #find sequence for isomiR
            
            ori_path = os.getcwd()
            
            tabT = tabfilewithtitle
            C1 = genomedataTCGA
            
            miRNA_list = tabT.title_box[tabT.title_dicX['miRNA_ID']-1]
            
            if len(set(miRNA_list)) == 0:
                print('miRNA input Error !')
                quit()
            elif len(set(miRNA_list)) == 1:
                
                #deal with unexpected miRNA
                excepted_miRNA = [('hsa-mir-1260', 'hsa-mir-1260a'), ('hsa-mir-147', 'hsa-mir-147a'), ('hsa-mir-190', 'hsa-mir-190a'),
                                  ('hsa-mir-3150', 'hsa-mir-3150a'), ('hsa-mir-323', 'hsa-mir-323a'), ('hsa-mir-378', 'hsa-mir-378a'),
                                  ('hsa-mir-514-1', 'hsa-mir-514a-1'),('hsa-mir-514-2', 'hsa-mir-514a-2'), ('hsa-mir-514-3', 'hsa-mir-514a-3'),
                                  ('hsa-mir-544', 'hsa-mir-544a'), ('hsa-mir-549', 'hsa-mir-549a'), ('hsa-mir-644', 'hsa-mir-644a'),
                                  ('hsa-mir-663', 'hsa-mir-663a')]
                excepted_dic = {x:y for x, y in excepted_miRNA}
                miRNA = "".join(list(set(miRNA_list)))
                if miRNA in excepted_dic.keys():
                    miRNA = excepted_dic[miRNA]
                else:
                    miRNA = miRNA
                
            else:
                print('miRNA is more then one type in this block.')
                quit()
            self.miRNA = miRNA
            
            print(self.miRNA+' in processing 1 !')
            sequence_list = []
            main_key_list = tabT.main_key_list
            
            left_location_list = []
            right_location_list = []
            for key in main_key_list:
                
                W = ""
                E = key.split(':')
                
                #take the chromosome obj
                os.chdir('/home/tin/Lib/Lib/hg19')
                chromosome = 'obj'
                
                W = "".join([x for x in C1.hg_obj_index if x.find(str(E[1]).upper()) != -1 and len(x) == len(str(E[1]).upper())])    
                
                if C1.hg_obj_list[C1.hg_obj_index.index(W)] == []:
                    print(str(C1.hg_obj_name[C1.hg_obj_index.index(W)])+' is opening !! (' + str(len([x for x in C1.hg_obj_list if x != []])) + '/25)')
                    chromosome = RaxLib.openFastAfile()
                    chromosome.open(C1.hg_obj_name[C1.hg_obj_index.index(W)])
                    C1.hg_obj_list[C1.hg_obj_index.index(W)] = chromosome
                
                else:
                    chromosome = C1.hg_obj_list[C1.hg_obj_index.index(W)]
                
                os.chdir(ori_path)
                
                #take the sequence from chromosome obj
                isomiR_sequence = ""
                left_location = "".join(E[2].split('-')[0])
                left_location_list.append(left_location)
                right_location = "".join(E[2].split('-')[1])
                right_location_list.append(right_location)
                if E[-1] == '+':
                    isomiR_sequence = chromosome.read(list(chromosome.FastA_dictionary.keys())[0])[int(left_location)-1:int(right_location)]
                    isomiR_sequence = RaxLib.ntComRev(isomiR_sequence, 'nC', 'nr', 'U', 'RNA')
                elif E[-1] == '-':
                    isomiR_sequence = chromosome.read(list(chromosome.FastA_dictionary.keys())[0])[int(left_location)-1:int(right_location)]
                    isomiR_sequence = RaxLib.ntComRev(isomiR_sequence, 'C', 'r', 'U', 'RNA')
                sequence_list.append(isomiR_sequence)
            
            left_location_list.sort()
            right_location_list.sort()
            genome_left = int(left_location_list[0]) - 100
            genome_right = int(right_location_list[-1]) + 100
            
            #add the sequence_list to table obj
            tabT.append('X', 'isomiR_sequence', sequence_list)
            #tabT.printtable()
            
            self.isomiR_tabfilewithtitle = tabT
            self.genome_left = genome_left
            self.genome_right = genome_right
            self.genome_chromosome = main_key_list[0].split(':')[1]
            self.genome_version = main_key_list[0].split(':')[0]
            self.genome_strand = main_key_list[0].split(':')[-1]
            
            if self.genome_right-self.genome_left < 300:
                self.situation = 'OK'
            else:
                self.situation = 'NO'
            
            genome_sequence = chromosome.read(list(chromosome.FastA_dictionary.keys())[0])[self.genome_left-1:self.genome_right]
            if main_key_list[0].split(':')[-1] == '+':
                genome_sequence = RaxLib.ntComRev(genome_sequence, 'nC', 'nr', 'U', 'RNA')
            elif main_key_list[0].split(':')[-1] == '-':
                genome_sequence = RaxLib.ntComRev(genome_sequence, 'C', 'r', 'U', 'RNA')
            self.genome_sequence = genome_sequence
            
            Num = "".join(filter(str.isdigit, self.miRNA))
            miRNA_ID = "".join([x for x in list(C1.percursor.FastA_dictionary.keys()) if Num in "".join(filter(str.isdigit, x))])
            if miRNA_ID == "":
                self.situation = 'NO'
            else:
                self.situation = 'OK'
            
            os.chdir(ori_path)
            
            return tabT
        
        def process2(self, tabfilewithtitle = 'table obj', genomedataTCGA = 'dataTCGA obj'): #define sequence_ali for isomiR
            
            ori_path = os.getcwd()
            
            tabT = tabfilewithtitle
            C1 = genomedataTCGA
            print(self.miRNA+' in processing 2 !')
            
            if self.situation == 'OK':
                
                main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                
                genome_sequence = self.genome_sequence
                
                percursor_sequence = genome_sequence.replace('T', 'U')
                
                #make sequence_ali of isomiR in List
                seq_ali = []
                for key in main_key_list:
                    
                    Y = ['.']*len(percursor_sequence)
                    
                    X = []
                    seq_of_table_1 = tabT.read('isomiR_sequence', key).upper().replace('T', 'U')
                    if seq_of_table_1 in percursor_sequence:
                        nt_index = percursor_sequence.find(seq_of_table_1)
                        
                        for nt in seq_of_table_1:
                            X.append(nt)
                        
                        i = 1
                        while i < len(X)+1:
                            Y[nt_index+i-1] = X[i-1]
                            i = i + 1
                    
                    seq_ali.append("".join(Y))
                
                #add the sequence_list to table obj
                tabT.append('X', 'isomiR_sequence_alignment', seq_ali)
                #tabT.printtable()
                
                self.isomiR_tabfilewithtitle = tabT
                
                os.chdir(ori_path)
                
                return tabT
                
            else:
                main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                seq_ali = ['na']*len(main_key_list)
                
                #add the sequence_list to table obj
                tabT.append('X', 'isomiR_sequence_alignment', seq_ali)
                #tabT.printtable()
                
                self.isomiR_tabfilewithtitle = tabT
                
                os.chdir(ori_path)
                
                return tabT
        
        def process5(self, tabfilewithtitle = 'table obj', genomedataTCGA = 'dataTCGA obj'): #define miRNA precursor and template
            
            ori_path = os.getcwd()
            
            tabT = tabfilewithtitle
            C1 = genomedataTCGA
            print(self.miRNA+' in processing 5 !')
            
            if self.situation == 'OK':
                miRNA = self.miRNA
                #define the percursor_sequence_h
                precursor_list = C1.matchtable.title_box[C1.matchtable.title_dicX['precursor']-1]
                precursor = [x for x in precursor_list if miRNA in x and len(x.split(' ')[0]) == len(miRNA)]
                
                if len(precursor) == 1:
                    percursor_sequence = C1.percursor.read("".join(precursor))
                    
                    #take percursor_genome_left and percursor_genome_right
                    location = self.genome_sequence.find(percursor_sequence)
                    
                    if location != -1:
                        
                        if self.genome_strand == '+':
                            self.percursor_genome_left = self.genome_left + location
                            self.percursor_genome_right = self.genome_left + location + len(percursor_sequence)
                        elif self.genome_strand == '-':
                            self.percursor_genome_left = self.genome_right - location - len(percursor_sequence) + 1
                            self.percursor_genome_right = self.genome_right - location + 1
                        
                        main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                        for key in main_key_list:
                            E = key.split(':')
                            #take the chromosome obj
                            os.chdir('/home/tin/Lib/Lib/hg19')
                            chromosome = 'obj'
                            
                            W = "".join([x for x in C1.hg_obj_index if str(E[1]).upper() in x and len(x) == len(str(E[1]).upper())])    
                            
                            if C1.hg_obj_list[C1.hg_obj_index.index(W)] == []:
                                chromosome = RaxLib.openFastAfile()
                                chromosome.open(C1.hg_obj_name[C1.hg_obj_index.index(W)])
                                C1.hg_obj_list[C1.hg_obj_index.index(W)] = chromosome
                                print(str(C1.hg_obj_name[C1.hg_obj_index.index(W)])+' is opened !! (' + str(len([x for x in C1.hg_obj_list if x != []])) + '/25)')
                            
                            else:
                                chromosome = C1.hg_obj_list[C1.hg_obj_index.index(W)]
                            
                            os.chdir(ori_path)
                        
                        if main_key_list[0].split(':')[-1] == '+':
                            genome_sequence = chromosome.read(list(chromosome.FastA_dictionary.keys())[0])[self.percursor_genome_left-4-1:self.percursor_genome_left+len(percursor_sequence)+4-1]
                            genome_sequence = RaxLib.ntComRev(genome_sequence, 'nC', 'nr', 'U', 'RNA')
                        elif main_key_list[0].split(':')[-1] == '-':
                            genome_sequence = chromosome.read(list(chromosome.FastA_dictionary.keys())[0])[self.percursor_genome_left-4-1:self.percursor_genome_left+len(percursor_sequence)+4-1]
                            genome_sequence = RaxLib.ntComRev(genome_sequence, 'C', 'r', 'U', 'RNA')
                        self.genome_sequence = genome_sequence
                        
                        #add the genome_sequence to table obj
                        genome_sequence_list = [genome_sequence]*len(main_key_list)
                        tabT.append('X', 'genome_sequence', genome_sequence_list)
                        #tabT.printtable()
                        
                        #add the percursor_sequence to table obj
                        percursor_sequence_list = ['....'+percursor_sequence+'....']*len(main_key_list)
                        tabT.append('X', 'percursor_sequence_ali', percursor_sequence_list)
                        #tabT.printtable()
                        
                        self.isomiR_tabfilewithtitle = tabT
                        self.situation = 'OK'
                        self.percursor_hairpin = "".join(precursor)
                        self.percursor_sequence = percursor_sequence
                        self.percursor_h_sequence_ali = '....'+percursor_sequence+'....'
                    
                    else:
                        main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                        #add the genome_sequence to table obj
                        genome_sequence_list = ['na']*len(main_key_list)
                        tabT.append('X', 'genome_sequence', genome_sequence_list)
                        #tabT.printtable()
                        #add the percursor_sequence to table obj
                        percursor_sequence_list = ['na']*len(main_key_list)
                        tabT.append('X', 'percursor_sequence_ali', percursor_sequence_list)
                        #tabT.printtable()
                        
                        self.isomiR_tabfilewithtitle = tabT
                        self.situation = 'NO'
                        self.percursor_hairpin = "".join('na')
                        self.percursor_sequence = 'na'
                        self.percursor_h_sequence_ali = "".join('na')
                        self.percursor_genome_left = 0
                        self.percursor_genome_right = 0
                        self.genome_sequence = 'na'
                    
                else:
                    main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                    #add the genome_sequence to table obj
                    genome_sequence_list = ['na']*len(main_key_list)
                    tabT.append('X', 'genome_sequence', genome_sequence_list)
                    #tabT.printtable()
                    #add the percursor_sequence to table obj
                    percursor_sequence_list = ['na']*len(main_key_list)
                    tabT.append('X', 'percursor_sequence_ali', percursor_sequence_list)
                    #tabT.printtable()
                    
                    self.isomiR_tabfilewithtitle = tabT
                    self.situation = 'NO'
                    self.percursor_hairpin = "".join('na')
                    self.percursor_sequence = 'na'
                    self.percursor_h_sequence_ali = "".join('na')
                    self.percursor_genome_left = 0
                    self.percursor_genome_right = 0
                    self.genome_sequence = 'na'
            
            else:
                main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                #add the genome_sequence to table obj
                genome_sequence_list = ['na']*len(main_key_list)
                tabT.append('X', 'genome_sequence', genome_sequence_list)
                #tabT.printtable()
                #add the percursor_sequence to table obj
                percursor_sequence_list = ['na']*len(main_key_list)
                tabT.append('X', 'percursor_sequence_ali', percursor_sequence_list)
                #tabT.printtable()
                
                self.isomiR_tabfilewithtitle = tabT
                self.situation = 'NO'
                self.percursor_hairpin = "".join('na')
                self.percursor_sequence = 'na'
                self.percursor_h_sequence_ali = "".join('na')
                self.percursor_genome_left = 0
                self.percursor_genome_right = 0
                self.genome_sequence = 'na'
            
            return tabT
        
        def process6(self, tabfilewithtitle = 'table obj', genomedataTCGA = 'dataTCGA obj'): #define mature_start_and_end
            
            ori_path = os.getcwd()
            
            tabT = tabfilewithtitle
            C1 = genomedataTCGA
            print(self.miRNA+' in processing 6 !')
            if self.situation == 'OK':
                
                Array = C1.locationtable.find(self.percursor_hairpin.split(' ')[1]+'_1', 'derives_from')
                for locationtable_key in Array:
                    locationtable_miRNAName = C1.locationtable.read('Name_hsaV20',locationtable_key[1])
                    
                    if '-5p' in locationtable_miRNAName:
                        self.mature_5p_start = C1.locationtable.read('start',locationtable_key[1])
                        self.mature_5p_end = C1.locationtable.read('end',locationtable_key[1])
                    elif '-3p' in locationtable_miRNAName:
                        self.mature_3p_start = C1.locationtable.read('start',locationtable_key[1])
                        self.mature_3p_end = C1.locationtable.read('end',locationtable_key[1])
                    else:
                        self.mature_3p_start = C1.locationtable.read('start',locationtable_key[1])
                        self.mature_3p_end = C1.locationtable.read('end',locationtable_key[1])
                
                self.situation == 'OK'
            
            else:
                self.situation = 'NO'
                self.mature_3p_start = 0
                self.mature_3p_end = 0
                self.mature_5p_start = 0
                self.mature_5p_end = 0
            
            return tabT
        
        def process7(self, tabfilewithtitle = 'table obj', genomedataTCGA = 'dataTCGA obj'): #define isomiR type
            
            ori_path = os.getcwd()
            
            tabT = tabfilewithtitle
            C1 = genomedataTCGA
            print(self.miRNA+' in processing 7 !')
            if self.situation == 'OK':
                
                isomiR_type_dic = {}
                isomiR_arm_dic = {}
                isomiR_type_list = []
                isomiR_arm_list = []
                
                main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                for key in main_key_list:
                    
                    key_start = float(key.split(':')[2].split('-')[0])
                    key_end = float(key.split(':')[2].split('-')[1])
                    
                    if key_start >= self.mature_5p_start-7 and key_end <= self.mature_5p_end+7:
                        isomiR_type_dic[key] = (int(key_start - self.mature_5p_start),int(key_end - self.mature_5p_end))
                        isomiR_arm_dic[key] = '5p'
                    elif key_start >= self.mature_3p_start-7 and key_end <= self.mature_3p_end+7:
                        isomiR_type_dic[key] = (int(key_start - self.mature_3p_start),int(key_end - self.mature_3p_end))
                        isomiR_arm_dic[key] = '3p'
                    else:
                        isomiR_type_dic[key] = 'non'
                        
                        miRNA_region_accession = tabT.read('miRNA_region', key)
                        if miRNA_region_accession == 'precursor':
                            isomiR_arm_dic[key] = 'na'
                        elif miRNA_region_accession == 'stemloop':
                            isomiR_arm_dic[key] = 'na'
                        elif miRNA_region_accession == 'unannotated':
                            isomiR_arm_dic[key] = 'undefined'
                        else:
                            isomiR_arm_dic[key] = 'undefined'
                
                if self.genome_strand == '+':
                    isomiR_type_dic = isomiR_type_dic
                    
                elif self.genome_strand == '-':
                    modify_dic = {}
                    for x in list(isomiR_type_dic.keys()):
                        if isomiR_type_dic[x] != 'non':
                            modify_dic[x] = ((isomiR_type_dic[x][1]*(-1)),(isomiR_type_dic[x][0]*(-1)))
                        else:
                            modify_dic[x] = 'non'
                    
                    isomiR_type_dic = modify_dic
                
                for key in main_key_list:
                    isomiR_type_list.append(isomiR_type_dic[key])
                    isomiR_arm_list.append(isomiR_arm_dic[key])
                
                #add the isomiR_type_list and isomiR_arm_list to table obj
                tabT.append('X', 'isomiR_type', isomiR_type_list)
                tabT.append('X', 'arm', isomiR_arm_list)
                #tabT.printtable()
                
                self.isomiR_tabfilewithtitle = tabT
                self.situation = 'OK'
            
            else:
                main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                isomiR_type_list = ['na']*len(main_key_list)
                isomiR_arm_list = ['na']*len(main_key_list)
                
                #add the isomiR_type_list and isomiR_arm_list to table obj
                tabT.append('X', 'isomiR_type', isomiR_type_list)
                tabT.append('X', 'arm', isomiR_arm_list)
                #tabT.printtable()
                
                self.isomiR_tabfilewithtitle = tabT
                self.situation = 'NO'
            
            return tabT
        
        def process8(self, tabfilewithtitle = 'table obj', genomedataTCGA = 'dataTCGA obj'): #open full isomiR class
            
            ori_path = os.getcwd()
            
            tabT = tabfilewithtitle
            C1 = genomedataTCGA
            
            miRNA_list = tabT.title_box[tabT.title_dicX['miRNA_ID']-1]
            
            if len({}.fromkeys(miRNA_list).keys()) == 0:
                print('miRNA input Error !')
                quit()
            elif len({}.fromkeys(miRNA_list).keys()) == 1:
                
                #deal with unexpected miRNA
                excepted_miRNA = [('hsa-mir-1260', 'hsa-mir-1260a'), ('hsa-mir-147', 'hsa-mir-147a'), ('hsa-mir-190', 'hsa-mir-190a'),
                                  ('hsa-mir-3150', 'hsa-mir-3150a'), ('hsa-mir-323', 'hsa-mir-323a'), ('hsa-mir-378', 'hsa-mir-378a'),
                                  ('hsa-mir-514-1', 'hsa-mir-514a-1'),('hsa-mir-514-2', 'hsa-mir-514a-2'), ('hsa-mir-514-3', 'hsa-mir-514a-3'),
                                  ('hsa-mir-544', 'hsa-mir-544a'), ('hsa-mir-549', 'hsa-mir-549a'), ('hsa-mir-644', 'hsa-mir-644a'),
                                  ('hsa-mir-663', 'hsa-mir-663a')]
                excepted_dic = {x:y for x, y in excepted_miRNA}
                miRNA = "".join(list(set(miRNA_list)))
                if miRNA in excepted_dic.keys():
                    miRNA = excepted_dic[miRNA]
                else:
                    miRNA = miRNA
                
            else:
                print('miRNA is more then one type in this block.')
                quit()
            self.miRNA = miRNA
            
            print(self.miRNA+' in processing 8 !')
            main_key_list = tabT.main_key_list
            
            self.isomiR_tabfilewithtitle = tabT
            self.genome_chromosome = main_key_list[0].split(':')[1]
            self.genome_version = main_key_list[0].split(':')[0]
            self.genome_strand = main_key_list[0].split(':')[-1]
            
            genome_sequence = self.isomiR_tabfilewithtitle.read('genome_sequence', main_key_list[0])
            self.genome_sequence = genome_sequence
            percursor_sequence_ali = self.isomiR_tabfilewithtitle.read('percursor_sequence_ali', main_key_list[0])
            self.percursor_h_sequence_ali = percursor_sequence_ali
            self.percursor_sequence = "".join([x for x in percursor_sequence_ali if x != '.'])
            self.percursor_hairpin = miRNA
            
            if genome_sequence == 'na' or self.percursor_h_sequence_ali == 'na':
                self.situation = 'NO'
            else:
                self.situation = 'OK'
            
            os.chdir(ori_path)
            
            return tabT
        
        #-----------------------------------------------------------------------------------------------
        
        def reportisomiR(self):
            
            display_list = []
            info = ':'.join([self.genome_version, self.genome_chromosome, str(self.genome_left)+'-'+str(self.genome_right), self.genome_strand])
            One_line = '\t'.join([self.miRNA, info])
            display_list.append(One_line)
            
            Two_line = self.percursor_h_sequence_ali
            display_list.append(Two_line)
            
            three_line = self.genome_sequence
            display_list.append(three_line)
            
            
            isomiRview_main_key_list = []
            M = {}
            if self.genome_strand == '+':
                
                isomiRview_main_key_list = self.isomiR_tabfilewithtitle.main_key_list
                
                T = []
                for x in isomiRview_main_key_list:
                    T.append((x, x))
                
                M = {x:y for x,y in T}
                
            elif self.genome_strand == '-':
                
                main_key_list = self.isomiR_tabfilewithtitle.main_key_list
                
                T = []
                for key in main_key_list:
                    #ex:'hg19:3:160122398-160122416:+'
                    W_ = key.split(':')[2]
                    W_start = W_.split('-')[0]
                    W_end = W_.split('-')[1]
                    new_start = str(1000000000-int(W_start))
                    new_end = str(1000000000-int(W_end))
                    new_key = new_end + '-' + new_start
                    
                    T.append((new_key, key))
                    isomiRview_main_key_list.append(new_key)
                
                M = {x:y for x,y in T}
                
            isomiRview_main_key_list.sort()
            
            E = []
            for new_key in isomiRview_main_key_list:
                key = M[new_key]
                
                result1 = self.isomiR_tabfilewithtitle.read('isomiR_sequence_alignment', key)
                result2 = self.isomiR_tabfilewithtitle.read('read_count', key)
                result3 = self.isomiR_tabfilewithtitle.read('reads_per_million_miRNA_mapped', key)
                result4 = self.isomiR_tabfilewithtitle.read('isomiR_type', key)
                
                E.append('\t'.join([result1, str(result2), str(result3), str(result4)]))
            
            four_line = '\n'.join(E)
            display_list.append(four_line)
            
            five_line = str(self.isomiR_tabfilewithtitle.len_Y)
            display_list.append(five_line)
            
            display = '\n'.join(display_list)+'\n'
            print(display)
            
            return display