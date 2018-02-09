


def processing(miRNA_matureFastA_file_address = 'path', statistic_Address = 'path'):
	import os
	
	#-----------------------------------------------------------------------------------------------
	#3 parameter:
	#miRNA_matureFastA_file_address = '/home/tin/Rax_sequence_with_matureV1/matureV1_-1/2_structure_assay_matureV1-1/matureV1-1FastA'	
	#statistic_Address = '/home/tin/Rax_sequence_with_matureV1/matureV1_-1/4_3p5plengthpercentage/result'
	
	
	#-----------------------------------------------------------------------------------------------
	#make the key list from miRNA_matureFastA_file_address
	listdir_name = []
	listdir_name = os.listdir(miRNA_matureFastA_file_address)
	
	
	#-----------------------------------------------------------------------------------------------
	#need pre-keep the objuct
	total_print = []
	check_threepandfivep_none = []
	check_threep_none = []
	check_fivep_none = []
	ALL_count = []
	
	total_string = ''
	key_number = 1
	B = ''
	c = 0                      #It is how much data be count
	
	V = 0
	ff = 0
	VV = []
	ALL_count_tital = ('miR_type'+'\t'+'5p_21nt'+'\t'+'5p_22nt'+'\t'+'5p_23nt'+'\t'
					  +'3p_21nt'+'\t'+'3p_22nt'+'\t'+'3p_23nt'+'\t'
					  +'3p5p_21nt'+'\t'+'3p5p_22nt'+'\t'+'3p5p_23nt'+'\t'
					  +'len_hsa_5p'+'\t'+'len_hsa_3p'+'\t'+'len_hsa_3p5p')
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
		
		
		
		#-----------------------------------------------------------------------------------------------
		#openfile and read FastA
		
		A = open(sample+'.fa','r')                          
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
		#select the key point in the 5 list, and throw out what you get
		#find miR-threep_list, miR-fivep_list and miR-non_list 
		i = 1
		W = ''                                         
		A = []    
		miR_threep_list = []
		miR_fivep_list = []
		miR_non_list = []
		miR_3p5p_list = []
		
		while  i < len(key) + 1:                            
			W = key[i-1]
						
			if W.find('3p') != -1:
				miR_threep_list.append(key[i-1])
			elif W.find('5p') != -1:
				miR_fivep_list.append(key[i-1])
			else:
				miR_non_list.append(key[i-1])
				
			i = i + 1                                    
	
			
		#print miR_fivep_list
			
		
		#-----------------------------------------------------------------------------------------------
		#define the compare sequence
		
		
		
		threep_compare_sequences = ''
		fivep_compare_sequences = ''
		threepfivep_compare_sequences = ''
		if len(miR_threep_list) == 0 and len(miR_fivep_list) == 0 :
			check_threepandfivep_none.append(sample + '\n')
			Z = [x for x in miR_non_list if not ''.join(x).find('hsa')]
			if len(Z) > 0:
				threepfivep_compare_sequences = sequence[M[Z[0]]][3:18]
			else:
				threepfivep_compare_sequences = sequence[M[miR_non_list[0]]][3:18]
				
			fivep_compare_sequences = ''
			threep_compare_sequences = ''
			#key_number = key_number + 1
			#continue
			
			
		elif len(miR_threep_list) == 0 and len(miR_fivep_list) != 0 :
			check_threep_none.append( sample + '\n' )
			Y = [x for x in miR_fivep_list if not ''.join(x).find('hsa')]
			if len(Y) > 0:
				fivep_compare_sequences = sequence[M[Y[0]]][3:18]			
			else :
				fivep_compare_sequences = sequence[M[miR_fivep_list[0]]][3:18]
				
			threep_compare_sequences = ''
			threepfivep_compare_sequences = ''
			
			
		elif len(miR_fivep_list) == 0 and len(miR_threep_list) != 0 :
			check_fivep_none.append( sample + '\n' )
			X = [x for x in miR_threep_list if not ''.join(x).find('hsa')]
			if len(X) > 0:
				threep_compare_sequences = sequence[M[X[0]]][3:18]			
			else :
				threep_compare_sequences = sequence[M[miR_threep_list[0]]][3:18]
				
			fivep_compare_sequences = ''
			threepfivep_compare_sequences = ''
			
			
		elif len(miR_threep_list) != 0 and len(miR_fivep_list) != 0 :
			Y = [x for x in miR_fivep_list if not ''.join(x).find('hsa')]
			X = [x for x in miR_threep_list if not ''.join(x).find('hsa')]
			if len(Y) > 0:
				fivep_compare_sequences = sequence[M[Y[0]]][3:18]			
			else :
				fivep_compare_sequences = sequence[M[miR_fivep_list[0]]][3:18]
			if len(X) > 0:
				threep_compare_sequences = sequence[M[X[0]]][3:18]			
			else :
				threep_compare_sequences = sequence[M[miR_threep_list[0]]][3:18]
			
			threepfivep_compare_sequences = ''
			
	
		else:
			check_threepandfivep_none.append(sample+'@@"' + '\n\n')
		
		
		#-----------------------------------------------------------------------------------------------
		#report the threep_compare_sequences and fivep_compare_sequences and threepfivep_compare_sequences
		
		if not [x for x in miR_threep_list if not ''.join(x).find('hsa')]:
			X_print = 'none'
		else:
			X_print = [x for x in miR_threep_list if not ''.join(x).find('hsa')][0]
		if not [x for x in miR_fivep_list if not ''.join(x).find('hsa')]:
			Y_print = 'none'
		else:
			Y_print = [x for x in miR_fivep_list if not ''.join(x).find('hsa')][0]
		if not [x for x in miR_non_list if not ''.join(x).find('hsa')]:
			Z_print = 'none'
		else:
			Z_print = [x for x in miR_non_list if not ''.join(x).find('hsa')][0]	
		
		#-----------------------------------------------------------------------------------------------
		#compare the sequencing order and sorting the miR-non_list into miR-threep_list and miR-fivep_list	
			
		i = 1
		unexpected_key = []
		W = ''
		while  i < len(miR_non_list) + 1:
			W = sequence[M[miR_non_list[i-1]]]
			
			if not threep_compare_sequences and fivep_compare_sequences :
				if W.find(fivep_compare_sequences) != -1:
					miR_fivep_list.append(miR_non_list[i-1])
				else:
					unexpected_key.append(miR_non_list[i-1])			
				
			elif not fivep_compare_sequences and threep_compare_sequences :
				if W.find(threep_compare_sequences) != -1:
					miR_threep_list.append(miR_non_list[i-1])			
				else:
					unexpected_key.append(miR_non_list[i-1])
					
			elif not threep_compare_sequences and not fivep_compare_sequences :
				if W.find(threepfivep_compare_sequences) != -1:
					miR_3p5p_list.append(miR_non_list[i-1])
					VV.append(miR_non_list[i-1]+'\n')
					
				else:
					unexpected_key.append(miR_non_list[i-1])
					
			else:
				if W.find(threep_compare_sequences) != -1:
					miR_threep_list.append(miR_non_list[i-1])
				elif W.find(fivep_compare_sequences) != -1:
					miR_fivep_list.append(miR_non_list[i-1])
				else:
					unexpected_key.append(miR_non_list[i-1])
			
			
			V = V + 1
			i = i + 1
			
	
		#print miR_fivep_list
		#print len(select_output)
	
		#-----------------------------------------------------------------------------------------------
		#there are 4 list can be used : miR_threep_list, miR_fivep_list, miR_3p5p_list and unexpected_key
		#-----------------------------------------------------------------------------------------------
		
		
		#-----------------------------------------------------------------------------------------------
		#The list cycle start here
		miR_FOUR_list_name = ['miR_threep_list', 'miR_fivep_list', 'miR_3p5p_list', 'unexpected_key']
		miR_FOUR_list = [miR_threep_list, miR_fivep_list, miR_3p5p_list, unexpected_key]
		
		#for FOUR_list in miR_FOUR_list:
		i = 1
		while i < len(miR_FOUR_list_name) + 1:
		
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
					
			if len(miR_FOUR_list[i-1]) == 0 :
				if miR_FOUR_list_name[i-1] == 'miR_threep_list':
					The_3p_21_percentage = 'N'
					The_3p_22_percentage = 'N'
					The_3p_23_percentage = 'N'
	
				elif miR_FOUR_list_name[i-1] == 'miR_fivep_list':
					The_5p_21_percentage = 'N'
					The_5p_22_percentage = 'N'
					The_5p_23_percentage = 'N'
	
				elif miR_FOUR_list_name[i-1] == 'miR_3p5p_list':
					The_3p5p_21_percentage = 'N'
					The_3p5p_22_percentage = 'N'
					The_3p5p_23_percentage = 'N'
				else:
					The_unexp_21_percentage = 'N'
					The_unexp_22_percentage = 'N'
					The_unexp_23_percentage = 'N'
				i = i + 1
				continue
			else:
				for keys in miR_FOUR_list[i-1]:
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
						
				display_count = (float(len(seq_15_list))+float(len(seq_16_list))+float(len(seq_17_list))+float(len(seq_18_list))+float(len(seq_19_list))+float(len(seq_20_list))+
						 float(len(seq_21_list))+float(len(seq_22_list))+float(len(seq_23_list))+float(len(seq_24_list))+float(len(seq_25_list))+float(len(seq_26_list))+
						 float(len(seq_27_list))+float(len(seq_28_list))+float(len(seq_29_list))+float(len(seq_30_list))+float(len(seq_31_list))+float(len(seq_32_list))+
						 float(len(seq_other_list)))
				print display_count
			
				if display_count != 0:
					if miR_FOUR_list_name[i-1] == 'miR_threep_list':
						The_3p_21_percentage = (float(len(seq_21_list))/display_count)*100
						The_3p_22_percentage = (float(len(seq_22_list))/display_count)*100
						The_3p_23_percentage = (float(len(seq_23_list))/display_count)*100
					
	
					elif miR_FOUR_list_name[i-1] == 'miR_fivep_list':
						The_5p_21_percentage = (float(len(seq_21_list))/display_count)*100
						The_5p_22_percentage = (float(len(seq_22_list))/display_count)*100
						The_5p_23_percentage = (float(len(seq_23_list))/display_count)*100

					elif miR_FOUR_list_name[i-1] == 'miR_3p5p_list':
						The_3p5p_21_percentage = (float(len(seq_21_list))/display_count)*100
						The_3p5p_22_percentage = (float(len(seq_22_list))/display_count)*100
						The_3p5p_23_percentage = (float(len(seq_23_list))/display_count)*100	
				else:
					display_count = 1
					if miR_FOUR_list_name[i-1] == 'miR_threep_list':
						The_3p_21_percentage = (float(len(seq_21_list))/display_count)*100
						The_3p_22_percentage = (float(len(seq_22_list))/display_count)*100
						The_3p_23_percentage = (float(len(seq_23_list))/display_count)*100
					
	
					elif miR_FOUR_list_name[i-1] == 'miR_fivep_list':
						The_5p_21_percentage = (float(len(seq_21_list))/display_count)*100
						The_5p_22_percentage = (float(len(seq_22_list))/display_count)*100
						The_5p_23_percentage = (float(len(seq_23_list))/display_count)*100

					elif miR_FOUR_list_name[i-1] == 'miR_3p5p_list':
						The_3p5p_21_percentage = (float(len(seq_21_list))/display_count)*100
						The_3p5p_22_percentage = (float(len(seq_22_list))/display_count)*100
						The_3p5p_23_percentage = (float(len(seq_23_list))/display_count)*100
						
						
						
		
			i = i + 1
		
		
		#print The_3p_21_percentage
		#print The_5p_21_percentage
		#print The_3p5p_21_percentage
		#print The_3p_22_percentage
		#print The_5p_22_percentage
		#print The_3p5p_22_percentage
		#print The_3p_23_percentage
		#print The_5p_23_percentage
		#print The_3p5p_23_percentage
		
		
	
		#-----------------------------------------------------------------------------------------------
		#re-print the FastA file
	
		i = 1
		W = ''                                            
		miR_threep = []
		miR_fivep = []
		unexpected = []
		miR_3p5p = []
	
		while i < len(miR_threep_list) + 1:
			W = "".join(miR_threep_list[i-1])
			miR_threep.append('>'+ key[M[W]] + " " + number[M[W]] + " " + species[M[W]] + " " + miR[M[W]] + "\n" + sequence[M[W]] ) 
			#miR_threep.append('>'+ '%-10s' % key[M[W]]  + "\t" + '%-25s' % sequence[M[W]] + str(len(sequence[M[W]])) + "\t" + '%s' % species[M[W]] ) 
			
			i = i + 1  
		i = 1
		while i < len(miR_fivep_list) + 1:
			W = "".join(miR_fivep_list[i-1])
			miR_fivep.append('>'+ key[M[W]] + " " + number[M[W]] + " " + species[M[W]] + " " + miR[M[W]] + "\n" + sequence[M[W]] ) 
			#miR_fivep.append('>'+ '%-10s' % key[M[W]]  + "\t" + '%-25s' % sequence[M[W]] + str(len(sequence[M[W]])) + "\t" + '%s' % species[M[W]] )
			
			i = i + 1 
		i = 1	
		while i < len(miR_3p5p_list) + 1:
			W = "".join(miR_3p5p_list[i-1])
			miR_3p5p.append('>'+ key[M[W]] + " " + number[M[W]] + " " + species[M[W]] + " " + miR[M[W]] + "\n" + sequence[M[W]] ) 
			#miR_3p5p.append('>'+ '%-10s' % key[M[W]]  + "\t" + '%-25s' % sequence[M[W]] + str(len(sequence[M[W]])) + "\t" + '%s' % species[M[W]] )
			
			i = i + 1	
		i = 1
		while i < len(unexpected_key) + 1:
			W = "".join(unexpected_key[i-1])
			unexpected.append('>'+ key[M[W]] + " " + number[M[W]] + " " + species[M[W]] + " " + miR[M[W]] + "\n" + sequence[M[W]] ) 
			#unexpected.append('>'+ '%-10s' % key[M[W]]  + "\t" + '%-25s' % sequence[M[W]] + str(len(sequence[M[W]])) + "\t" + '%s' % species[M[W]] )
			
			i = i + 1			
			
		output_miR_threep = "\n".join(miR_threep)
		output_miR_fivep = "\n".join(miR_fivep) 
		output_miR_3p5p = "\n".join(miR_3p5p)
		output_unexpected = "\n".join(unexpected)
		output_countthreep = str(len(miR_threep))
		output_countfivep = str(len(miR_fivep))
		output_count3p5p = str(len(miR_3p5p))
		output_countunexpected = str(len(unexpected))
		output_totalcount = str(len(miR_threep)+len(miR_fivep)+len(miR_3p5p)+len(unexpected))
		
		
		All = ("miR-3p:" + "\n" + output_miR_threep + "\n\n" 
			  + "miR-5p:" + "\n" + output_miR_fivep + "\n\n"
			  + "miR-3p5p:" + "\n" + output_miR_3p5p + "\n\n"
			  + "unexpected:" + "\n" + output_unexpected)
		
		print All      
		print output_totalcount
		
		
		
		#-----------------------------------------------------------------------------------------------
		#find the hsa length
		
		i = 1
		W = ''
		hsa_3p = []
		while i < len(miR_threep_list) + 1:
			W = "".join(miR_threep_list[i-1])
			if W.find("hsa") != -1:
				hsa_3p.append(miR_threep_list[i-1])
				
			i = i +1
		
		i = 1
		W = ''
		hsa_5p = []
		while i < len(miR_fivep_list) + 1:
			W = "".join(miR_fivep_list[i-1])
			if W.find("hsa") != -1:
				hsa_5p.append(miR_fivep_list[i-1])
						
			i = i +1
			
		i = 1
		W = ''
		hsa_3p5p = []
		while i < len(miR_3p5p_list) + 1:
			W = "".join(miR_3p5p_list[i-1])
			if W.find("hsa") != -1:
				hsa_3p5p.append(miR_3p5p_list[i-1])
						
			i = i +1
		
		
		len_hsa_3p = ''
		len_hsa_5p = ''
		len_hsa_3p5p = ''
		
		if not hsa_3p:
			len_hsa_3p = 'N'
		else:
			len_hsa_3p = str(len(sequence[M["".join(hsa_3p[0])]])-1)
			
		if not hsa_5p:
			len_hsa_5p = 'N'
		else:
			len_hsa_5p = str(len(sequence[M["".join(hsa_5p[0])]])-1)
			
		if not hsa_3p5p:
			len_hsa_3p5p = 'N'
		else:
			len_hsa_3p5p = str(len(sequence[M["".join(hsa_3p5p[0])]])-1)
		
		ALL_count.append(sample +'\t'+ str(The_5p_21_percentage) +'\t'+ str(The_5p_22_percentage) +'\t'+ str(The_5p_23_percentage) +'\t'
					     + str(The_3p_21_percentage) +'\t'+ str(The_3p_22_percentage) +'\t'+ str(The_3p_23_percentage) +'\t'
					     + str(The_3p5p_21_percentage) +'\t'+ str(The_3p5p_22_percentage) +'\t'+ str(The_3p5p_23_percentage) +'\t'
					     + len_hsa_3p +'\t'+ len_hsa_5p +'\t'+ len_hsa_3p5p)
					 
		
	
		
		total_print.append( All + '\n\n' + output_totalcount )
		
		c = c + len(output_totalcount)
		ff = ff + len(miR_3p5p)
	
		key_number = key_number + 1
	
		
	
	
	os.chdir(statistic_Address)
	#-----------------------------------------------------------------------------------------------
	#writeout statistic
	total_string = "\n".join(total_print)
	f = open('total_log.txt','w')
	f.write(total_string+'\n')
	f.write('total_dir:' + str(len(ALL_count))+'\n')
	f.write('totaldata_infile:' + str(c)+'\n')
	f.close()
	
	
	
	#-----------------------------------------------------------------------------------------------
	#writeout datalist
	H = ''
	H = "\n".join(ALL_count)
	f = open('total_datalist.txt','w')
	f.write(ALL_count_tital+'\n')
	f.write(H)
	f.close()
	
	
	#-----------------------------------------------------------------------------------------------
	#writeout 3pand5p none
	none_string_3p = ''.join(check_threep_none)
	none_string_5p = ''.join(check_fivep_none)
	CHECK =''.join(VV)
	none_string_3pand5p = ''.join(check_threepandfivep_none)
	H = ('none_3P' + '\n' + none_string_3p + '\n' + str(len(check_threep_none)) + '\n'  
	 +'none_5P' + '\n' + none_string_5p + '\n' + str(len(check_fivep_none))+ '\n' 
	 +'none_3pand5p' + '\n' + none_string_3pand5p + '\n' + str(len(check_threepandfivep_none))+ '\n'
	 +'fix in miR_3p5p' + '\n' + CHECK + '\n' + str(V) + '\n'+ str(ff))
	
	f = open('3p5pnone.txt','w')
	f.write(H)
	f.close()
	
	return
	

