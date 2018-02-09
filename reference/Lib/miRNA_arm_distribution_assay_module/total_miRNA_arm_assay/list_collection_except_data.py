

def processing(miRNA_matureFastA_file_address = 'path', statistic_Address = 'path', selection_condiction = 1, Output_file_name = 'file.fa'):

	## the '/home/tin/Rax/mature/matureFastA' is a parameter
	## the P is parameter
	## the selected condition is a parameter
	
	import os
	
	#4 parameter:	
	#miRNA_matureFastA_file_address = '/home/tin/Rax_sequence_with_matureV1/matureV1FastA'
	#statistic_Address = '/home/tin/Rax_sequence_with_matureV1/matureV1_-1/1_make_matureV1-1/result'
	#selection_condiction = 1                 #keep more then 'selection_condiction'
	#Output_file_name = 'matureV1-1.fa'
	
	#-----------------------------------------------------------------------------------------------
	#make the key list from miRNA_matureFastA_file_address
	listdir_name = []
	listdir_name = os.listdir(miRNA_matureFastA_file_address)
	
	#-----------------------------------------------------------------------------------------------
	#need pre-keep the objuct
	total_print = []
	All_for_FastAdata = []
	ALL_count = []
	total_string = ''
	key_number = 1
	B = ''
	c = 0
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
		print P
		Savefiledir = P
		
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
		#re-print the FastA file
	
		i = 1
		W = ''                                            
		E = []
		while i < len(key) + 1:
			W = "".join(key[i-1])
			E.append('>'+ key[M[W]] + " " + number[M[W]] + " " + species[M[W]] + " " + miR[M[W]] + "\n" + sequence[M[W]] )#FastA
			#E.append( miR[M[W]] + "\t\t" + sequence[M[W]]+ "\t\t" + species[M[W]]  )   
				
			i = i + 1                                    
		
		print str(len(E))
		
		if len(E) > selection_condiction:                                    #the selected condition
			All_cut = "\n".join(E)  
			All_for_FastAdata.append(All_cut)
			c = c + len(E)
			
			
		All = "\n".join(E)                               
			
		print All                                        
		
		ALL_count.append( sample + '\t' +'has:'+ '\t' + str(len(key))) 
	
		
		
		
		total_print.append( All + "\n" + str(len(E))+ "\n" + P  )
		
		
		
		
		
	
		key_number = key_number + 1
		
		
	#-----------------------------------------------------------------------------------------------
	#writeout statistic
	#total_log.txt included the content, amount, amount of created dir, amount of content
	
	
	P = ''
	P = statistic_Address
	Savefiledir = P
	os.chdir(Savefiledir)
	
	
	total_string = "\n".join(total_print)
	All_for_FastAdata_string = "\n".join(All_for_FastAdata)
	
	f = open('total_log.txt','w')
	f.write(total_string+'\n')
	f.write(str(len(E))+'\n')
	f.write('total_dir:' + str(len(ALL_count))+'\n')
	f.write('total_count_data:' + str(c)+'\n')
	f.close()
	
	#-----------------------------------------------------------------------------------------------
	#writeout datalist
	
	f = open(Output_file_name,'w')
	f.write(All_for_FastAdata_string+'\n')
	f.write(str(c))
	f.close()	

	return