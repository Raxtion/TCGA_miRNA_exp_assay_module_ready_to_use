import os


def OpenmatureFastA(mature):
    
    #-----------------------------------------------------------------------------------------------
    #openfile and read FastA
    
    mature_key = []
    mature_number = []
    mature_species = []
    mature_miR = []
    mature_sequence = []
    mature_M = {}
    
    total = open(mature,'r').read()
    A = total.split('>')
    del A[0]
    
    i = 1
    for x in A:
            
        mature_key.append("".join("".join(x).split('\n')[0]).split(" ")[0])
        mature_number.append("".join("".join(x).split('\n')[0]).split(" ")[1])
        mature_species.append(" ".join("".join("".join(x).split('\n')[0]).split(" ")[2:4]))
        mature_miR.append("".join("".join(x).split('\n')[0]).split(" ")[4])    
        mature_sequence.append("".join("".join(x).split('\n')[1]))
        mature_M["".join("".join(x).split('\n')[0]).split(" ")[0]]=i-1
        i=i+1
    
    #-----------------------------------------------------------------------------------------------
    #there are 5 list can be used : mature_key, mature_number, mature_species, mature_miR and mature_sequence
    #there is 1 dictionary can be used : mature_M
    #-----------------------------------------------------------------------------------------------
    return mature_M, mature_key, mature_number, mature_species, mature_miR, mature_sequence, 


def lengthseperated(M = 'dictionary', key = 'list', sequence = 'list'):
	#-----------------------------------------------------------------------------------------------
	#set sequence list box
	
	seq_15_list = []
	seq_16_list = []
	seq_17_list = []
	seq_18_list = []
	seq_19_list = []
	seq_20_list = []
	seq_21_list = []
	seq_22_list = []
	seq_23_list = []
	seq_24_list = []
	seq_25_list = []
	seq_26_list = []
	seq_27_list = []
	seq_28_list = []
	seq_29_list = []
	seq_30_list = []
	seq_31_list = []
	seq_32_list = []
	
	seq_other_list = []
	
	#-----------------------------------------------------------------------------------------------
	#seperate data accroding mature_key
	
	for keys in key:
		if len(sequence[M[keys]]) == 15:
			seq_15_list.append(keys)
		elif len(sequence[M[keys]]) == 16:
			seq_16_list.append(keys)
		elif len(sequence[M[keys]]) == 17:
			seq_17_list.append(keys)
		elif len(sequence[M[keys]]) == 18:
			seq_18_list.append(keys)
		elif len(sequence[M[keys]]) == 19:
			seq_19_list.append(keys)
		elif len(sequence[M[keys]]) == 20:
			seq_20_list.append(keys)
		elif len(sequence[M[keys]]) == 21:
			seq_21_list.append(keys)
		elif len(sequence[M[keys]]) == 22:
			seq_22_list.append(keys)
		elif len(sequence[M[keys]]) == 23:
			seq_23_list.append(keys)
		elif len(sequence[M[keys]]) == 24:
			seq_24_list.append(keys)
		elif len(sequence[M[keys]]) == 25:
			seq_25_list.append(keys)
		elif len(sequence[M[keys]]) == 26:
			seq_26_list.append(keys)	
		elif len(sequence[M[keys]]) == 27:
			seq_27_list.append(keys)	
		elif len(sequence[M[keys]]) == 28:
			seq_28_list.append(keys)
		elif len(sequence[M[keys]]) == 29:
			seq_29_list.append(keys)
		elif len(sequence[M[keys]]) == 30:
			seq_30_list.append(keys)
		elif len(sequence[M[keys]]) == 31:
			seq_31_list.append(keys)
		elif len(sequence[M[keys]]) == 32:
			seq_32_list.append(keys)
		else:
			seq_other_list.append(keys)
			
			
	#-----------------------------------------------------------------------------------------------	
	#rebuilding output
	len_range_list = ['15mer'
			 ,'16mer'
			 ,'17mer'
			 ,'18mer'
			 ,'19mer'
			 ,'20mer'
			 ,'21mer'
			 ,'22mer'
			 ,'23mer'
			 ,'24mer'
			 ,'25mer'
			 ,'26mer'
			 ,'27mer'
			 ,'28mer'
			 ,'29mer'
			 ,'30mer'
			 ,'31mer'
			 ,'32mer'
			 ,'other']
	display = [str(len(seq_15_list))
		  ,str(len(seq_16_list))
		  ,str(len(seq_17_list))
		  ,str(len(seq_18_list))
		  ,str(len(seq_19_list))
		  ,str(len(seq_20_list))
		  ,str(len(seq_21_list))
		  ,str(len(seq_22_list))
		  ,str(len(seq_23_list))
		  ,str(len(seq_24_list))
		  ,str(len(seq_25_list))
		  ,str(len(seq_26_list))
		  ,str(len(seq_27_list))
		  ,str(len(seq_28_list))
		  ,str(len(seq_29_list))
		  ,str(len(seq_30_list))
		  ,str(len(seq_31_list))
		  ,str(len(seq_32_list))
		  ,str(len(seq_other_list))]
	display_count = str(len(seq_15_list)+len(seq_16_list)+len(seq_17_list)+len(seq_18_list)+len(seq_19_list)+len(seq_20_list)
			+len(seq_21_list)+len(seq_22_list)+len(seq_23_list)+len(seq_24_list)+len(seq_25_list)+len(seq_26_list)
			+len(seq_27_list)+len(seq_28_list)+len(seq_29_list)+len(seq_30_list)+len(seq_31_list)+len(seq_32_list)
			+len(seq_other_list))
	
	return display, display_count, len_range_list



def listmake(miRNA_matureFastA_file_address='path', Open_which_file_in_matureFastA='arm'):
	
	#-----------------------------------------------------------------------------------------------
	#4 parameter:
	#Open_which_file_in_matureFastA = '5p'
	#miRNA_matureFastA_file_address = '/home/tin/Rax_sequence_with_matureV1/Separated_by_miRNA/matureV1_-1/3_3por5p_determination/matureV1-1_3p5pFastA'
	#statistic_Address = '/home/tin/Rax_sequence_with_matureV1/lengh_distribution/matureV1/matureV1-1_3p5p_len/matureV1-1_5p_len_list/result'
	#Output_file_name = 'matureV1-1_5p_seq_list.txt'
	
	
	#-----------------------------------------------------------------------------------------------
	#make the key list from miRNA_matureFastA_file_address
	listdir_name = []
	listdir_name = os.listdir(miRNA_matureFastA_file_address)
	
	
	
	#-----------------------------------------------------------------------------------------------
	#set sequence list box
	
	seq_15_list = []
	seq_16_list = []
	seq_17_list = []
	seq_18_list = []
	seq_19_list = []
	seq_20_list = []
	seq_21_list = []
	seq_22_list = []
	seq_23_list = []
	seq_24_list = []
	seq_25_list = []
	seq_26_list = []
	seq_27_list = []
	seq_28_list = []
	seq_29_list = []
	seq_30_list = []
	seq_31_list = []
	seq_32_list = []
	
	seq_other_list = []
	
	
	#-----------------------------------------------------------------------------------------------
	#need pre-keep the objuct
	
	key_number = 1
	B = ''
	#-----------------------------------------------------------------------------------------------
	#name_cycle start as there
	
	while  key_number < len(listdir_name) + 1: 
		B = listdir_name[key_number-1]
		sample = "".join(B)
	
		#print key[0]
	
		
		#-----------------------------------------------------------------------------------------------
		#chang to target dir
		
		P = ''
		P = miRNA_matureFastA_file_address+'/'+sample
		#print P
		Savefiledir = P
		#if not os.path.isdir(Savefiledir):
		#	os.makedirs(Savefiledir)
		os.chdir(P)
		
		sample = sample.split('_')[0]
		
		#-----------------------------------------------------------------------------------------------
		#openfile and read FastA
		
		A = open(sample+'_'+Open_which_file_in_matureFastA+'.fa','r')                          
		total = A.read()
		T = total.replace('>','')                        
		data_ID = T.split('\n')[::2]
		del data_ID[-1]     #In the origin data, last row is not the sample just "\n"
		sequence = T.split('\n')[1::2]
	
	
	
		data_X = " ".join(data_ID)
		#let the list_data_ID join into string
		#and then split and cut
		key = data_X.split(" ")[0::5]
		number = data_X.split(" ")[1::5]
		first_species = data_X.split(" ")[2::5]
		second_species = data_X.split(" ")[3::5]
		miR = data_X.split(" ")[4::5]
	
	
		#recombin the species
		i = 1                                            
		E=[]                                             
		while  i < len(first_species) + 1 :                            
			E.append(first_species[i-1] + " " + second_species[i-1] )   
			i = i + 1                                    
		species_map = "\n".join(E)
		species = species_map.split("\n")
	
		#-----------------------------------------------------------------------------------------------
		#there are 5 list can be used : key, number, species, miR and sequence
		#-----------------------------------------------------------------------------------------------
	
	
		#-----------------------------------------------------------------------------------------------
		#let the key list fix into dictionary under their order
		i = 1
		W = ''
		M = {}
		while i < len(key) + 1:
			W = "".join(key[i-1])
			M[W]=i-1
			i=i+1
			
		#print M[key[5000]]
		#print key[5000] 
		
		
		#-----------------------------------------------------------------------------------------------
		#separate the length of sequence
		if key == ['']:
			key_number = key_number + 1
			continue
		else:
			for keys in key:
				if len(sequence[M[keys]]) == 15:
					seq_15_list.append(keys)
				elif len(sequence[M[keys]]) == 16:
					seq_16_list.append(keys)
				elif len(sequence[M[keys]]) == 17:
					seq_17_list.append(keys)
				elif len(sequence[M[keys]]) == 18:
					seq_18_list.append(keys)
				elif len(sequence[M[keys]]) == 19:
					seq_19_list.append(keys)
				elif len(sequence[M[keys]]) == 20:
					seq_20_list.append(keys)
				elif len(sequence[M[keys]]) == 21:
					seq_21_list.append(keys)
				elif len(sequence[M[keys]]) == 22:
					seq_22_list.append(keys)
				elif len(sequence[M[keys]]) == 23:
					seq_23_list.append(keys)
				elif len(sequence[M[keys]]) == 24:
					seq_24_list.append(keys)
				elif len(sequence[M[keys]]) == 25:
					seq_25_list.append(keys)
				elif len(sequence[M[keys]]) == 26:
					seq_26_list.append(keys)	
				elif len(sequence[M[keys]]) == 27:
					seq_27_list.append(keys)	
				elif len(sequence[M[keys]]) == 28:
					seq_28_list.append(keys)
				elif len(sequence[M[keys]]) == 29:
					seq_29_list.append(keys)
				elif len(sequence[M[keys]]) == 30:
					seq_30_list.append(keys)
				elif len(sequence[M[keys]]) == 31:
					seq_31_list.append(keys)
				elif len(sequence[M[keys]]) == 32:
					seq_32_list.append(keys)
				else:
					seq_other_list.append(keys)
			
			key_number = key_number + 1
		
		
		
	#-----------------------------------------------------------------------------------------------	
	#rebuilding output
	len_range_list = ['15mer'
			 ,'16mer'
			 ,'17mer'
			 ,'18mer'
			 ,'19mer'
			 ,'20mer'
			 ,'21mer'
			 ,'22mer'
			 ,'23mer'
			 ,'24mer'
			 ,'25mer'
			 ,'26mer'
			 ,'27mer'
			 ,'28mer'
			 ,'29mer'
			 ,'30mer'
			 ,'31mer'
			 ,'32mer'
			 ,'other']
	display = [str(len(seq_15_list))
		  ,str(len(seq_16_list))
		  ,str(len(seq_17_list))
		  ,str(len(seq_18_list))
		  ,str(len(seq_19_list))
		  ,str(len(seq_20_list))
		  ,str(len(seq_21_list))
		  ,str(len(seq_22_list))
		  ,str(len(seq_23_list))
		  ,str(len(seq_24_list))
		  ,str(len(seq_25_list))
		  ,str(len(seq_26_list))
		  ,str(len(seq_27_list))
		  ,str(len(seq_28_list))
		  ,str(len(seq_29_list))
		  ,str(len(seq_30_list))
		  ,str(len(seq_31_list))
		  ,str(len(seq_32_list))
		  ,str(len(seq_other_list))]
	display_count = str(len(seq_15_list)+len(seq_16_list)+len(seq_17_list)+len(seq_18_list)+len(seq_19_list)+len(seq_20_list)
			+len(seq_21_list)+len(seq_22_list)+len(seq_23_list)+len(seq_24_list)+len(seq_25_list)+len(seq_26_list)
			+len(seq_27_list)+len(seq_28_list)+len(seq_29_list)+len(seq_30_list)+len(seq_31_list)+len(seq_32_list)
			+len(seq_other_list))
	
	return display, display_count, len_range_list


def listprocessing(miRNA_matureFastA_file_address = 'path', Open_which_file_in_matureFastA = 'arm', Output_file_name = None, statistic_Address = 'path'):
	
	
	#-----------------------------------------------------------------------------------------------
	#Open_which_file_in_matureFastA = '5p'
	#miRNA_matureFastA_file_address = '/home/tin/Rax_sequence_with_matureV1/Separated_by_miRNA/matureV1_-1/3_3por5p_determination/matureV1-1_3p5pFastA'
	#statistic_Address = '/home/tin/Rax_sequence_with_matureV1/lengh_distribution/matureV1/matureV1-1_3p5p_len/matureV1-1_5p_len_list/result'
	#Output_file_name = 'matureV1-1_5p_seq_list.txt'
	
	#-----------------------------------------------------------------------------------------------
	P1 = listmake(miRNA_matureFastA_file_address, Open_which_file_in_matureFastA)
	display = P1[0]
	display_count = P1[1]
	len_range_list = P1[2]
	
	#-----------------------------------------------------------------------------------------------
	E = []
	i = 1
	while i < len(len_range_list) + 1:
		E.append(len_range_list[i-1]+':'+'\t'+str(display[i-1]))
		i = i + 1
	
		
	os.chdir(statistic_Address)
	#-----------------------------------------------------------------------------------------------
	#writeout statistic
	if Output_file_name != None:
		f = open(Output_file_name,'w')
		f.write(str('\n'.join(E))+'\n')
		f.write('total_sample:' + str(display_count))
		f.close()
	else:
		pass
	
	return display, display_count, len_range_list

def fileprocessing(mature = 'file.fa', Output_file_name = None, statistic_Address = 'path'):
	
	#-----------------------------------------------------------------------------------------------
	P1 = OpenmatureFastA(mature)
	mature_M = P1[0]
	mature_key = P1[1]	
	mature_sequence = P1[5]
	
	P2 = lengthseperated(mature_M, mature_key, mature_sequence)
	display = P2[0]
	display_count = P2[1]
	len_range_list = P2[2]
	
	
	#-----------------------------------------------------------------------------------------------
	E = []
	i = 1
	while i < len(len_range_list) + 1:
		E.append(len_range_list[i-1]+':'+'\t'+str(display[i-1]))
		i = i + 1
	
		
	os.chdir(statistic_Address)
	#-----------------------------------------------------------------------------------------------
	#writeout statistic
	if Output_file_name != None:
		f = open(Output_file_name,'w')
		f.write(str('\n'.join(E))+'\n')
		f.write('total_sample:' + str(display_count))
		f.close()
	else:
		pass
	
	
	return display, display_count, len_range_list


