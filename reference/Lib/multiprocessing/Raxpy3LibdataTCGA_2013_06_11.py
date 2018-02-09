'''
history =  2013.05.13, 2013_05_19, 2013_05_24
'''
#-----------------------------------------------------------------------------------------------
import os, Raxpy3Libbasic_2013_06_10
#-----------------------------------------------------------------------------------------------
RaxLib = Raxpy3Libbasic_2013_06_10

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
                        '21', '22', 'M', 'X', 'Y']
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
        mature_miRNA.open('hsaV1.fa')
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
        matchtable.open('group_table8.txt', 'order')
        self.matchtable = matchtable
        
        os.chdir(ori_path)
        
        return
    
    #-----------------------------------------------------------------------------------------------
    
    def openlocationtable(self):
        
        ori_path = os.getcwd()
        #take hsaprecursor and hsamiRNA locationtable
        os.chdir('/home/tin/Lib/Lib/hsagff3_location/')
        gff3table = RaxLib.tabfilewithtitle()
        gff3table.open('hsagff3modifyhsaV2.txt', 'order')
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
            
            if len({}.fromkeys(miRNA_list).keys()) == 0:
                print('miRNA input Error !')
                quit()
            elif len({}.fromkeys(miRNA_list).keys()) == 1:
                miRNA = "".join({}.fromkeys(miRNA_list).keys())
            else:
                print('miRNA is more then one type in this block.')
                quit()
            self.miRNA = miRNA
            
            print(self.miRNA+' in processing 1 !')
            sequence_list = []
            main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
            
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
                    isomiR_sequence = chromosome.read(chromosome.FastA_dictionary.keys()[0])[int(left_location)-1:int(right_location)]
                elif E[-1] == '-':
                    isomiR_sequence = chromosome.read(chromosome.FastA_dictionary.keys()[0])[int(left_location)-1:int(right_location)]
                    isomiR_sequence = RaxLib.ntComRev(isomiR_sequence, 'Ucase', 'r')
                sequence_list.append(isomiR_sequence.replace('T', 'U'))
            
            left_location_list.sort()
            right_location_list.sort()
            genome_left = int(left_location_list[0]) - 50
            genome_right = int(right_location_list[-1]) + 50
            
            #add the sequence_list to table obj
            tabT.append('X', 'isomiR_sequence', sequence_list)
            tabT.printtable()
            
            self.isomiR_tabfilewithtitle = tabT
            self.genome_left = genome_left
            self.genome_right = genome_right
            self.genome_chromosome = main_key_list[0].split(':')[1]
            self.genome_version = main_key_list[0].split(':')[0]
            self.genome_strand = main_key_list[0].split(':')[-1]
            
            if self.genome_right-self.genome_left < 200:
                self.situation = 'OK'
            else:
                self.situation = 'NO'
            
            genome_sequence = chromosome.read(chromosome.FastA_dictionary.keys()[0])[self.genome_left-1:self.genome_right]
            if main_key_list[0].split(':')[-1] == '+':
                genome_sequence = genome_sequence
            elif main_key_list[0].split(':')[-1] == '-':
                genome_sequence = RaxLib.ntComRev(genome_sequence, 'Ucase', 'r')
            self.genome_sequence = genome_sequence.upper()
            
            Num = "".join(filter(str.isdigit, self.miRNA))
            miRNA_ID = "".join([x for x in C1.mature_miRNA.FastA_dictionary.keys() if Num in x])
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
                tabT.printtable()
                
                self.isomiR_tabfilewithtitle = tabT
                
                os.chdir(ori_path)
                
                return tabT
                
            else:
                main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                seq_ali = ['na']*len(main_key_list)
                
                #add the sequence_list to table obj
                tabT.append('X', 'isomiR_sequence_alignment', seq_ali)
                tabT.printtable()
                
                self.isomiR_tabfilewithtitle = tabT
                
                os.chdir(ori_path)
                
                return tabT
        
        def process3(self, tabfilewithtitle = 'table obj', genomedataTCGA = 'dataTCGA obj'): #It is broken
            
            ori_path = os.getcwd()
            
            tabT = tabfilewithtitle
            C1 = genomedataTCGA
            
            if self.situation == 'OK':
                
                main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                mature_miRNA_key_list = [x.lower() for x in C1.mature_miRNA.FastA_dictionary.keys()]
                
                #definite the isomiR type
                ##find 'MIMA'
                in_frame_list = []
                unannotated_list = []
                out_frame_list = []
                for key in main_key_list:
                    
                    result = tabT.read('miRNA_region', key)
                    
                    if 'unannotated' in result and [x for x in mature_miRNA_key_list if filter(str.isdigit, self.miRNA) in x] != []:
                        unannotated_list.append(key)
                    elif result.find('unannotated') == -1:
                        if 'MIMA' in result:
                            in_frame_list.append(key)
                        elif result.find('MIMA') == -1:
                            out_frame_list.append(key)
                        else:
                            print(str(miRNA)+' '+str(key)+' '+'is no value to defind isomiR type')
                            quit()
                    else:
                        print(str(miRNA)+' '+str(key)+' '+'has unexcept type')
                        quit()
                
                #print in_frame_list
                
                ##unexcept
                if unannotated_list == []:
                    
                    ##define miRNA sequence and arm
                    miRNA_5p_sequence = {}
                    miRNA_3p_sequence = {}
                    for key in in_frame_list:
                        result = tabT.read('miRNA_region', key)
                        Num = result.split(',')[1]
                        
                        miRNA_ID = "".join([x for x in C1.mature_miRNA.FastA_dictionary.keys() if Num in x])
                        
                        if '5p' in miRNA_ID:
                            miRNA_5p_sequence[miRNA_ID] = C1.mature_miRNA.FastA_dictionary[miRNA_ID].replace('T', 'U')
                        elif '3p' in miRNA_ID:
                            miRNA_3p_sequence[miRNA_ID] = C1.mature_miRNA.FastA_dictionary[miRNA_ID].replace('T', 'U')
                        else:
                            print(str(self.miRNA)+' '+"does't find the miNRA arm")
                            quit()
                    
                    ##define arm group
                    miRNA_5p_list = []
                    miRNA_3p_list = []
                    for key in in_frame_list:
                        result = tabT.read('miRNA_region', key)
                        Num = result.split(',')[1]
                        
                        if len(miRNA_5p_sequence.keys()) != 0 and len(miRNA_3p_sequence.keys()) != 0:
                            if Num in "".join(miRNA_5p_sequence.keys()[0]):
                                miRNA_5p_list.append(key)
                            elif Num in "".join(miRNA_3p_sequence.keys()[0]):
                                miRNA_3p_list.append(key)
                            else:
                                print(str(miRNA_ID)+' '+key+' '+"can't group in arm group !!")
                                quit()
                        
                        elif len(miRNA_5p_sequence.keys()) != 0 and len(miRNA_3p_sequence.keys()) == 0:
                            if Num in "".join(miRNA_5p_sequence.keys()[0]):
                                miRNA_5p_list.append(key)
                            else:
                                print(str(miRNA_ID)+' '+key+' '+"can't group in arm group !!")
                                quit()
                        
                        elif len(miRNA_5p_sequence.keys()) == 0 and len(miRNA_3p_sequence.keys()) != 0:
                            if Num in "".join(miRNA_3p_sequence.keys()[0]):
                                miRNA_3p_list.append(key)
                            else:
                                print(str(miRNA_ID)+' '+key+' '+"can't group in arm group !!")
                                quit()
                        
                        else:
                            print('The '+self.miRNA+' is fucking shit!')
                            quit()
                    
                    ##define isomiR type
                    ###find 0,0 type
                    zerozero_5p_type = ""
                    zerozero_3p_type = ""
                    zerozero_5p_start_number = 0
                    zerozero_5p_end_number = 0
                    zerozero_3p_start_number = 0
                    zerozero_3p_end_number = 0
                    ali_sequence_list = tabT.title_box[tabT.title_dicX['isomiR_sequence_alignment']-1]
                    if ali_sequence_list[0] != 'na':
                        for ali_sequence in ali_sequence_list:
                            
                            if len(miRNA_5p_sequence.keys()) != 0 and len(miRNA_3p_sequence.keys()) != 0:
                                W5p = miRNA_5p_sequence[miRNA_5p_sequence.keys()[0]]
                                W3p = miRNA_3p_sequence[miRNA_3p_sequence.keys()[0]]
                                if "".join([x for x in ali_sequence.split('.') if x != ""]) == W5p:
                                    zerozero_5p_type = ali_sequence
                                    
                                    Array = ali_sequence.split('.')
                                    seq_in_array = "".join([x for x in Array if x != ""])
                                    zerozero_5p_start_number = Array.index(seq_in_array)
                                    
                                    Array.reverse()
                                    seq_in_array = "".join([x for x in Array if x != ""])
                                    zerozero_5p_end_number = Array.index(seq_in_array)
                                
                                elif "".join([x for x in ali_sequence.split('.') if x != ""]) == W3p:
                                    zerozero_3p_type = ali_sequence
                                    
                                    Array = ali_sequence.split('.')
                                    seq_in_array = "".join([x for x in Array if x != ""])
                                    zerozero_3p_start_number = Array.index(seq_in_array)
                                    
                                    Array.reverse()
                                    seq_in_array = "".join([x for x in Array if x != ""])
                                    zerozero_3p_end_number = Array.index(seq_in_array)
                                    
                                else:
                                    zerozero_5p_type = '....'+W5p+'....'
                                    zerozero_3p_type = '....'+W3p+'....'
                                    zerozero_5p_start_number = 4
                                    zerozero_5p_end_number = 4
                                    zerozero_3p_start_number = 4
                                    zerozero_3p_end_number = 4
                                
                                ###find other type
                                isomiR_start_type_dic = {}
                                isomiR_end_type_dic = {}
                                for key in miRNA_5p_list:
                                    
                                    result = tabT.read('isomiR_sequence_alignment', key)
                                    Array = result.split('.')
                                    isomiR_start_type_dic[key] = Array.index("".join([x for x in Array if x != ""])) - zerozero_5p_start_number
                                    
                                    Array.reverse()
                                    isomiR_end_type_dic[key] = Array.index("".join([x for x in Array if x != ""])) - zerozero_5p_end_number
                                
                                for key in miRNA_3p_list:
                                    
                                    result = tabT.read('isomiR_sequence_alignment', key)
                                    Array = result.split('.')
                                    isomiR_start_type_dic[key] = Array.index("".join([x for x in Array if x != ""])) - zerozero_3p_start_number
                                    
                                    Array.reverse()
                                    isomiR_end_type_dic[key] = Array.index("".join([x for x in Array if x != ""])) - zerozero_3p_end_number
                            
                            elif len(miRNA_5p_sequence.keys()) != 0 and len(miRNA_3p_sequence.keys()) == 0:
                                W5p = miRNA_5p_sequence[miRNA_5p_sequence.keys()[0]]
                                if "".join([x for x in ali_sequence.split('.') if x != ""]) == W5p:
                                    zerozero_5p_type = ali_sequence
                                    
                                    Array = ali_sequence.split('.')
                                    seq_in_array = "".join([x for x in Array if x != ""])
                                    zerozero_5p_start_number = Array.index(seq_in_array)
                                    
                                    Array.reverse()
                                    seq_in_array = "".join([x for x in Array if x != ""])
                                    zerozero_5p_end_number = Array.index(seq_in_array)
                                else:
                                    zerozero_5p_type = '....'+W5p+'....'
                                    zerozero_5p_start_number = 4
                                    zerozero_5p_end_number = 4
                                
                                ###find other type
                                isomiR_start_type_dic = {}
                                isomiR_end_type_dic = {}
                                for key in miRNA_5p_list:
                                    
                                    result = tabT.read('isomiR_sequence_alignment', key)
                                    Array = result.split('.')
                                    isomiR_start_type_dic[key] = Array.index("".join([x for x in Array if x != ""])) - zerozero_5p_start_number
                                    
                                    Array.reverse()
                                    isomiR_end_type_dic[key] = Array.index("".join([x for x in Array if x != ""])) - zerozero_5p_end_number
                            
                            elif len(miRNA_5p_sequence.keys()) == 0 and len(miRNA_3p_sequence.keys()) != 0:
                                W3p = miRNA_3p_sequence[miRNA_3p_sequence.keys()[0]]
                                if "".join([x for x in ali_sequence.split('.') if x != ""]) == W3p:
                                    zerozero_3p_type = ali_sequence
                                    
                                    Array = ali_sequence.split('.')
                                    seq_in_array = "".join([x for x in Array if x != ""])
                                    zerozero_3p_start_number = Array.index(seq_in_array)
                                    
                                    Array.reverse()
                                    seq_in_array = "".join([x for x in Array if x != ""])
                                    zerozero_3p_end_number = Array.index(seq_in_array)
                                else:
                                    zerozero_3p_type = '....'+W3p+'....'
                                    zerozero_3p_start_number = 4
                                    zerozero_3p_end_number = 4
                                
                                ###find other type
                                isomiR_start_type_dic = {}
                                isomiR_end_type_dic = {}
                                
                                for key in miRNA_3p_list:
                                    
                                    result = tabT.read('isomiR_sequence_alignment', key)
                                    Array = result.split('.')
                                    isomiR_start_type_dic[key] = Array.index("".join([x for x in Array if x != ""])) - zerozero_3p_start_number
                                    
                                    Array.reverse()
                                    isomiR_end_type_dic[key] = Array.index("".join([x for x in Array if x != ""])) - zerozero_3p_end_number
                            
                            else:
                                print('The '+self.miRNA+' is fucking shit!')
                                quit()
                        
                        
                        
                        ###bulid new dic and list
                        isomiR_type_dic = {}
                        isomiR_type_list = []
                        for key in main_key_list:
                            
                            if isomiR_start_type_dic.has_key(key) == True:
                                isomiR_type_dic[key] = ','.join([str(isomiR_start_type_dic[key]),str(isomiR_end_type_dic[key]*(-1))])
                            elif isomiR_start_type_dic.has_key(key) == False:
                                isomiR_type_dic[key] = 'non'
                        
                        for key in main_key_list:
                            
                            isomiR_type_list.append(isomiR_type_dic[key])
                        
                        
                        
                        
                    else:
                        isomiR_type_list = ['na']*len(main_key_list)
                    
                    #add the isomiR_type_list to table obj
                    tabT.append('X', 'isomiR_type', isomiR_type_list)
                    tabT.printtable()
                
                elif unannotated_list != []:
                    main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                    isomiR_type_list = ['na']*len(main_key_list)
                    
                    #add the isomiR_type_list to table obj
                    tabT.append('X', 'isomiR_type', isomiR_type_list)
                    tabT.printtable()
                
                else:
                    print('different type!!')
                    quit()
                
                self.isomiR_tabfilewithtitle = tabT
                
                print(str(self.miRNA) + ' is done !! ')
                
                return tabT
            
            else:
                main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                isomiR_type_list = ['na']*len(main_key_list)
                
                #add the isomiR_type_list to table obj
                tabT.append('X', 'isomiR_type', isomiR_type_list)
                tabT.printtable()
                
                print(str(self.miRNA) + ' is done !! ')
                
                return tabT
        
        def process4(self, tabfilewithtitle = 'table obj', genomedataTCGA = 'dataTCGA obj'): #It is broken
            
            ori_path = os.getcwd()
            
            tabT = tabfilewithtitle
            C1 = genomedataTCGA
            print(self.miRNA+' in processing 4 !')
            
            if self.situation == 'OK':
                #define the percursor_sequence_h
                miRNA = self.miRNA
                
                percursor_ID_list = [x for x in C1.percursor.FastA_dictionary.keys() if 'hsa' in x]
                E = []
                for percursor_ID in percursor_ID_list:
                    
                    if miRNA in percursor_ID and [x for x in percursor_ID_list if filter(str.isdigit, miRNA) in x] != []:
                        E.append(percursor_ID)
                
                
                if len(E) != 1:
                    self.situation = 'NO'
                    self.percursor_hairpin = 'na'
                    self.percursor_sequence = 'na'
                    self.percursor_h_sequence_ali = "".join('na')
                    print(self.miRNA + " can't define the percursor_sequence ! 2")
                    
                    return tabT
                
                else:
                    self.percursor_hairpin = "".join(E[0])
                    percursor_sequence = C1.percursor.read("".join(E[0]))
                
                print(percursor_sequence)
                
                #take percursor_genome_left and percursor_genome_right
                location = self.genome_sequence.replace('T', 'U').find(percursor_sequence)
                if location != -1:
                    self.percursor_genome_left = self.genome_left + location
                    self.percursor_genome_right = self.genome_left + len(percursor_sequence)
                    
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
                    
                    
                    genome_sequence = chromosome.read(chromosome.FastA_dictionary.keys()[0])[self.percursor_genome_left-4-1:self.percursor_genome_left+len(percursor_sequence)+4]
                    if main_key_list[0].split(':')[-1] == '+':
                        genome_sequence = genome_sequence
                    elif main_key_list[0].split(':')[-1] == '-':
                        genome_sequence = RaxLib.ntComRev(genome_sequence, 'Ucase', 'r')
                    self.genome_sequence = genome_sequence.upper().replace('T', 'U')
                    
                    
                    self.situation = 'OK'
                    self.percursor_sequence = percursor_sequence
                    self.percursor_h_sequence_ali = '....'+percursor_sequence+'....'
                    
                    return tabT
                
                else:
                    self.situation = 'NO'
                    self.percursor_sequence = percursor_sequence
                    self.percursor_h_sequence_ali = '....'+percursor_sequence+'....'
                    print(self.miRNA + " can't define the percursor_sequence ! 3")
                    
                    return tabT
            
            else:
                self.situation = 'NO'
                self.percursor_hairpin = 'na'
                self.percursor_sequence = 'na'
                self.percursor_h_sequence_ali = "".join('na')
                self.percursor_genome_left = 0
                self.percursor_genome_right = 0
                
                return tabT
            
            
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
                    location = self.genome_sequence.replace('T', 'U').find(percursor_sequence)
                    if location != -1:
                        self.percursor_genome_left = self.genome_left + location
                        self.percursor_genome_right = self.genome_left + len(percursor_sequence)
                        
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
                        
                        
                        genome_sequence = chromosome.read(chromosome.FastA_dictionary.keys()[0])[self.percursor_genome_left-4-1:self.percursor_genome_left+len(percursor_sequence)+4]
                        if main_key_list[0].split(':')[-1] == '+':
                            genome_sequence = genome_sequence
                        elif main_key_list[0].split(':')[-1] == '-':
                            genome_sequence = RaxLib.ntComRev(genome_sequence, 'Ucase', 'r')
                        self.genome_sequence = genome_sequence.upper().replace('T', 'U')
                    
                    
                    self.situation = 'OK'
                    self.percursor_hairpin = "".join(precursor)
                    self.percursor_sequence = percursor_sequence
                    self.percursor_h_sequence_ali = '....'+percursor_sequence+'....'
                    
                else:
                    self.situation = 'NO'
                    self.percursor_hairpin = "".join('na')
                    self.percursor_sequence = 'na'
                    self.percursor_h_sequence_ali = "".join('na')
                    self.percursor_genome_left = 0
                    self.percursor_genome_right = 0
            
            else:
                self.situation = 'NO'
                self.percursor_hairpin = "".join('na')
                self.percursor_sequence = 'na'
                self.percursor_h_sequence_ali = "".join('na')
                self.percursor_genome_left = 0
                self.percursor_genome_right = 0
            
            return tabT
        
        def process6(self, tabfilewithtitle = 'table obj', genomedataTCGA = 'dataTCGA obj'): #define mature_start_and_end
            
            ori_path = os.getcwd()
            
            tabT = tabfilewithtitle
            C1 = genomedataTCGA
            print(self.miRNA+' in processing 6 !')
            if self.situation == 'OK':
                
                Array = C1.locationtable.find(self.percursor_hairpin.split(' ')[1]+'_1', 'derives_from')
                for locationtable_key in Array:
                    locationtable_miRNAName = C1.locationtable.read('Name_hsaV2',locationtable_key[1])
                    
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
                self.mature_3p_start = 0.0
                self.mature_3p_end = 0.0
                self.mature_5p_start = 0.0
                self.mature_5p_end = 0.0
            
            return tabT
        
        def process7(self, tabfilewithtitle = 'table obj', genomedataTCGA = 'dataTCGA obj'): #define isomiR type
            
            ori_path = os.getcwd()
            
            tabT = tabfilewithtitle
            C1 = genomedataTCGA
            print(self.miRNA+' in processing 7 !')
            if self.situation == 'OK':
                
                isomiR_type_dic = {}
                isomiR_type_list = []
                
                main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                for key in main_key_list:
                    
                    key_start = float(key.split(':')[2].split('-')[0])
                    key_end = float(key.split(':')[2].split('-')[1])
                    
                    if key_start >= self.mature_5p_start-7 and key_end <= self.mature_5p_end+7:
                        isomiR_type_dic[key] = ((key_start - self.mature_5p_start),(key_end - self.mature_5p_end))
                    elif key_start >= self.mature_3p_start-7 and key_end <= self.mature_3p_end+7:
                        isomiR_type_dic[key] = ((key_start - self.mature_3p_start),(key_end - self.mature_3p_end))
                    else:
                        isomiR_type_dic[key] = 'non'
                
                for key in main_key_list:
                    isomiR_type_list.append(isomiR_type_dic[key])
                
                #add the sequence_list to table obj
                tabT.append('X', 'isomiR_type', isomiR_type_list)
                tabT.printtable()
                
                self.isomiR_tabfilewithtitle = tabT
                self.situation = 'OK'
            
            else:
                main_key_list = tabT.title_box[tabT.title_dicX['isoform_coords']-1]
                isomiR_type_list = ['na']*len(main_key_list)
                
                #add the sequence_list to table obj
                tabT.append('X', 'isomiR_type', isomiR_type_list)
                tabT.printtable()
                
                self.isomiR_tabfilewithtitle = tabT
                self.situation = 'NO'
            
            return tabT
        
        #-----------------------------------------------------------------------------------------------
        
        def reportisomiR(self):
            
            display_list = []
            info = ':'.join([self.genome_version, self.genome_chromosome, str(self.genome_left)+'-'+str(self.genome_right), self.genome_strand])
            One_line = '\t'.join([self.miRNA, info])
            display_list.append(One_line)
            
            Two_line = self.percursor_hairpin_ali
            display_list.append(Two_line)
            
            three_line = self.genome_sequence
            display_list.append(three_line)
            
            main_key_list = self.isomiR_tabfilewithtitle.title_box[self.isomiR_tabfilewithtitle.title_dicX['isoform_coords']-1]
            E = []
            for key in main_key_list:
                result1 = self.isomiR_tabfilewithtitle.read('isomiR_sequence_alignment', key)
                result2 = self.isomiR_tabfilewithtitle.read('read_count', key)
                result3 = self.isomiR_tabfilewithtitle.read('reads_per_million_miRNA_mapped', key)
                result4 = self.isomiR_tabfilewithtitle.read('isomiR_type', key)
                
                E.append('\t'.join([result1, str(result2), str(result3), result4]))
            four_line = '\n'.join(E)
            display_list.append(four_line)
            
            five_line = str(self.isomiR_tabfilewithtitle.len_Y)
            display_list.append(five_line)
            
            display = '\n'.join(display_list)+'\n'
            print(display)
            
            return display