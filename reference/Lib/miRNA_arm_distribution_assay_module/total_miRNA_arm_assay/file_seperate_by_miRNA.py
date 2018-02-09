




def processing(mature, seperation_Address, statistic_Address):
	import os
	#4 parameter
	#open_file_name = mature	
	seperation_Address = seperation_Address+"".join(mature.split('.')[0])+'FastA'
	os.makedirs(seperation_Address)
	#statistic_Address = '/home/tin/Rax_sequence_with_matureV1/result'
	
	
	A = open(mature,'r')
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
	M_ = {}
	while i < len(key) + 1:
		W = ("".join(key[i-1]), "".join(number[i-1]))
		M_[W]=i-1
		i=i+1
		
	#print M_[key[5000]]
	#print key[5000] 
	
	
	#-----------------------------------------------------------------------------------------------
	#select two key_type in to new list ONE and TWO
	ONE = []                #group for 'hsa-miR-2299'
	TWO = []                #group for 'hsa-miR2299'
	THREE = []              #group for 'hsa-miR-2299-3p'
	FOUR = []               #group for unexpected
	bantam = []             #group for 'bantam' miRNA
	i = 1
	W = ''                                         
	
	while  i < len(key) + 1:                            
		W = (key[i-1], number[i-1])
		Z = ''.join([x for x in W[0].split("-")[1:3][0] if x.isdigit()])
		
		if W[0].find('bantam') != -1:
			bantam.append((key[i-1], number[i-1]))
			i = i + 1
			continue
		else:
			if len(W[0].split("-")) <= 3:
				if len(Z) == 0 :
					ONE.append((key[i-1], number[i-1]))		
				elif len(Z) != 0 :
					TWO.append((key[i-1], number[i-1]))
			
			elif len(W[0].split("-")) > 3:
				X = W[0].split("-")[3].find('p')
				if X == -1 :
					THREE.append((key[i-1], number[i-1]))
				else:
					if len(Z) == 0 :
						ONE.append((key[i-1], number[i-1]))		
					elif len(Z) != 0 :
						TWO.append((key[i-1], number[i-1]))
			
			else:
				FOUR.append((key[i-1], number[i-1]))
			i = i + 1                                    
	
			
			
	i = 1
	W = ''
	while i < len(bantam) + 1: 
		W = bantam[i-1]
		
		if len(W[0].split("-")) <= 2:
			TWO.append((key[M_[W]], number[M_[W]]))
				
		elif len(W[0].split("-")) > 2:
			X = W[0].split("-")[2].find('p')
			if X != -1 :
				TWO.append((key[M_[W]], number[M_[W]]))
			else:
				ONE.append((key[M_[W]], number[M_[W]]))		
				
		
		else:
			FOUR.append((key[M_[W]], number[M_[W]]))
		i = i + 1                       
	
		
		
	
	#print ONE
	#print TWO
	#print THREE
	#print FOUR
	#print bantam
	
	
	
	#-----------------------------------------------------------------------------------------------
	#process in ONE
	total_print_ONE = []
	ALL_count_ONE = []
	total_string = ''
	ONE_number = 1
	B = ''
	c = 0
	while  ONE_number < len(ONE) + 1: 
		B = ONE[ONE_number-1]
		A = B[0].split("-")[1:3]
		sample = "-".join(A)
	
		#print ONE[0]
		
		#-----------------------------------------------------------------------------------------------
		#check is the sample has been read
		dir_list = os.listdir(seperation_Address)
		join_dir_list = '@'.join(dir_list)
		X = join_dir_list.find(sample + '@')
		if X == -1 :
			
		
		#-----------------------------------------------------------------------------------------------
		#select the ONE point in the 5 list, and throw out what you get 
			i = 1
			W = ''                                         
			select_output = []
			A = []         
													  #Q = raw_input("please input the ID:'\n'")
			while  i < len(ONE) + 1: 
				W = ONE[i-1]
				X = W[0].find(sample)
			
				if X != -1:
					A = W[0].split("-")[1:3]                 #start from the second(1) to in front the fourth(3)
					W = "-".join(A)                       #join them then confirm the length of you want
				
					if len(W) == len(sample):
						select_output.append(ONE[i-1])
					
				i = i + 1                                    
	
		#print select_output
		#print len(select_output)
	
	
	
	
		#-----------------------------------------------------------------------------------------------
		#re-print the FastA file
	
			select_output.sort()
			i = 1
			W = ''                                            
			E = []
			while i < len(select_output) + 1:
				W = select_output[i-1]
				E.append('>'+ key[M_[W]] + " " + number[M_[W]] + " " + species[M_[W]] + " " + miR[M_[W]] + "\n" + sequence[M_[W]] ) 
				#E.append( miR[M_[W]] + "\t\t" + sequence[M_[W]]+ "\t\t" + species[M_[W]]  )   
			
				i = i + 1                                    
			All = "\n".join(E)                               
			print All                                        
			print str(len(E))
			
			
			ALL_count_ONE.append( sample + '\t' + 'has :' + '\t' + str(len(select_output))) 
			
			
			#-----------------------------------------------------------------------------------------------
			#writeout to a .fa file
				
			P = ''
			P = seperation_Address+'/'+sample
			Savefiledir = P
			print P
			if not os.path.isdir(Savefiledir):
				os.makedirs(Savefiledir)
			os.chdir(P)
			
			
			
			W = ''
			W = sample + ".fa"
			f = open(W,'w')     
			f.write(All+'\n')	
			f.write(str(len(E)))
			f.close()
			
			
			
			total_print_ONE.append( All + "\n" + str(len(E)) )
			
	
			#-----------------------------------------------------------------------------------------------
			#count totaldata be write in file
			c = c + len(E)
			
				
			ONE_number = ONE_number + 1
		
		else:
			ONE_number = ONE_number + 1
			continue
		
	d = c 
		
		
	#-----------------------------------------------------------------------------------------------
	#process in TWO
	total_print_TWO = []
	ALL_count_TWO = []
	total_string = ''
	TWO_number = 1
	B = ''
	c = 0
	while  TWO_number < len(TWO) + 1: 
		B = TWO[TWO_number-1]
		A = B[0].split("-")[1:2]
		sample = "-".join(A)
	
		#print TWO[0]
			
		#-----------------------------------------------------------------------------------------------
		#check is the sample has been read
		dir_list = os.listdir(seperation_Address)
		join_dir_list = '@'.join(dir_list)
		X = join_dir_list.find(sample + '@')
		if X == -1 :
		
		#-----------------------------------------------------------------------------------------------
		#select the TWO point in the 5 list, and throw out what you get 
			i = 1
			W = ''                                         
			select_output = []
			A = []         
													  #Q = raw_input("please input the ID:'\n'")
			while  i < len(TWO) + 1:                            
				W = TWO[i-1]
				X = W[0].find(sample)
			
				if X != -1:
					A = W[0].split("-")[1:2]                 #start from the second(1) to in front the fourth(3)
					W = "-".join(A)                       #join them then confirm the length of you want
				
					if len(W) == len(sample):
						select_output.append(TWO[i-1])
					
				i = i + 1                                    
	
		#print select_output
		#print len(select_output)
	
	
	
	
		#-----------------------------------------------------------------------------------------------
		#re-print the FastA file
	
			select_output.sort()
			i = 1
			W = ''                                            
			E = []
			while i < len(select_output) + 1:
				W = select_output[i-1]
				E.append('>'+ key[M_[W]] + " " + number[M_[W]] + " " + species[M_[W]] + " " + miR[M_[W]] + "\n" + sequence[M_[W]] ) 
				#E.append( miR[M_[W]] + "\t\t" + sequence[M_[W]]+ "\t\t" + species[M_[W]]  )   
			
				i = i + 1                                    
			All = "\n".join(E)                               
			print All                                        
			print str(len(E))
			
			
			ALL_count_TWO.append( sample + '\t' + 'has :' + '\t' + str(len(select_output)) ) 
	
			
			#-----------------------------------------------------------------------------------------------
			#writeout to a .fa file
				
			P = ''
			P = seperation_Address+'/'+sample
			Savefiledir = P
			print P
			if not os.path.isdir(Savefiledir):
				os.makedirs(Savefiledir)
			os.chdir(P)
			
			
			
			W = ''
			W = sample + ".fa"
			f = open(W,'w')     
			f.write(All+'\n')	
			f.write(str(len(E)))
			f.close()
			
			
			total_print_TWO.append( All + "\n" + str(len(E))  )
			
	
			#-----------------------------------------------------------------------------------------------
			#count totaldata be write in file
			c = c + len(E)
			
				
			TWO_number = TWO_number + 1
			
		
		else:
			TWO_number = TWO_number + 1
			continue
	ee = d + c
	
	
	#-----------------------------------------------------------------------------------------------
	#process in THREE
	total_print_THREE = []
	ALL_count_THREE = []
	total_string = ''
	THREE_number = 1
	B = ''
	c = 0
	while  THREE_number < len(THREE) + 1: 
		B = THREE[THREE_number-1]
		A = B[0].split("-")[1:4]
		sample = "-".join(A)
	
		#print THREE[0]
		
		#-----------------------------------------------------------------------------------------------
		#check is the sample has been read
		dir_list = os.listdir(seperation_Address)
		join_dir_list = '@'.join(dir_list)
		X = join_dir_list.find(sample + '@')
		if X == -1 :
		
		#-----------------------------------------------------------------------------------------------
		#select the THREE point in the 5 list, and throw out what you get 
			i = 1
			W = ''                                         
			select_output = []
			A = []         
													  #Q = raw_input("please input the ID:'\n'")
			while  i < len(THREE) + 1:                            
				W = THREE[i-1]
				X = W[0].find(sample)
			
				if X != -1:
					A = W[0].split("-")[1:4]                 #start from the second(1) to in front the fourth(3)
					W = "-".join(A)                       #join them then confirm the length of you want
				
					if len(W) == len(sample):
						select_output.append(THREE[i-1])
					
				i = i + 1                                    
	
		#print select_output
		#print len(select_output)
	
	
	
	
		#-----------------------------------------------------------------------------------------------
		#re-print the FastA file
	
			select_output.sort()
			i = 1
			W = ''                                            
			E = []
			while i < len(select_output) + 1:
				W = select_output[i-1]
				E.append('>'+ key[M_[W]] + " " + number[M_[W]] + " " + species[M_[W]] + " " + miR[M_[W]] + "\n" + sequence[M_[W]] ) 
				#E.append( miR[M_[W]] + "\t\t" + sequence[M_[W]]+ "\t\t" + species[M_[W]]  )   
			
				i = i + 1                                    
			All = "\n".join(E)                               
			print All                                        
			print str(len(E))
			
			
			ALL_count_THREE.append( sample + '\t' + 'has :' + '\t' + str(len(select_output)) ) 
	
			
			#-----------------------------------------------------------------------------------------------
			#writeout to a .fa file
				
			P = ''
			P = seperation_Address+'/'+sample
			Savefiledir = P
			print P
			if not os.path.isdir(Savefiledir):
				os.makedirs(Savefiledir)
			os.chdir(P)
			
			
			
			W = ''
			W = sample + ".fa"
			f = open(W,'w')     
			f.write(All+'\n')	
			f.write(str(len(E)))
			f.close()
			
			
			total_print_THREE.append( All + "\n" + str(len(E))  )
			
	
			#-----------------------------------------------------------------------------------------------
			#count totaldata be write in file
			c = c + len(E)
			
				
			THREE_number = THREE_number + 1
			
		else:
			THREE_number = THREE_number + 1
			continue
	ddd = ee + c
	
	
	
	
	
	
	total_print = []
	ALL_count = []
	total_print = total_print_ONE + total_print_TWO + total_print_THREE
	ALL_count = ALL_count_ONE + ALL_count_TWO + ALL_count_THREE
	
	
	#-----------------------------------------------------------------------------------------------
	#delete the repeat item	
	ALL_count = {}.fromkeys(ALL_count).keys()
	
	
	
	
	
	os.chdir(statistic_Address)
	#-----------------------------------------------------------------------------------------------
	#writeout statistic
	total_string = "\n".join(total_print)
	f = open('total_log.txt','w')
	f.write(total_string+'\n')
	f.write(str(len(E))+'\n')
	f.write('total_dir:' + str(len(ALL_count))+'\n')
	f.write('totaldata_infile:' + str(ddd)+'\n')
	f.close()
	
	#-----------------------------------------------------------------------------------------------
	#writeout datalist
	H = ''
	H = "\n".join(ALL_count)
	f = open('total_datalist.txt','w')
	f.write(H)
	f.close()
	
	#-----------------------------------------------------------------------------------------------
	#writeout FOUR unexpected data
	H = ''
	H = "\n".join(FOUR)
	
	f = open('unexpected.txt','w')
	f.write(H)
	f.write(str(len(FOUR)))
	f.close()
	
	return seperation_Address