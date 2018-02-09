import os

ori_path = os.getcwd()

def running(ori_path = '', _input = 'TCGASTAD', RaxLib = 'Raxpy3Libbasic_2013_09_14'):
    import os, pickle
    
    #-----------------------------------------------------------------------------------------------
    #parameter
    RaxLib = RaxLib
    
    
    
    #-----------------------------------------------------------------------------------------------
    
    
    #-----------------------------------------------------------------------------------------------
    #openfile
    TabT_mixarm = RaxLib.tabfilewithtitle()
    TabT_mixarm.open(_input+'_NT_top3_mixarm_0to1.txt', 'mir')
    f = open(_input+'_NT_top3_mixarm_0to1.txt.pickletable', 'wb')
    pickle.dump(TabT_mixarm, f, )
    f.close()
    
    TabT_mixarm.XYreverse()
    f = open(_input+'_mixarm_0to1_Xdatasat_Ymir.txt', 'w')
    f.write(TabT_mixarm.report())
    f.close()
    f = open(_input+'_mixarm_0to1_Xdatasat_Ymir.txt.pickletable', 'wb')
    pickle.dump(TabT_mixarm, f, )
    f.close()
    
    TabT_twoarm = RaxLib.tabfilewithtitle()
    TabT_twoarm.open(_input+'_NT_top3_twoarm_0to1.txt', 'mir')
    f = open(_input+'_NT_top3_twoarm_0to1.txt.pickletable', 'wb')
    pickle.dump(TabT_twoarm, f, )
    f.close()
    
    TabT_twoarm.XYreverse()
    f = open(_input+'_twoarm_0to1_Xdatasat_Ymir.txt', 'w')
    f.write(TabT_twoarm.report())
    f.close()
    f = open(_input+'_twoarm_0to1_Xdatasat_Ymir.txt.pickletable', 'wb')
    pickle.dump(TabT_twoarm, f, )
    f.close()
    
    '''
    #quit()
    #-----------------------------------------------------------------------------------------------
    #openfile
    #TabT_mixarm = RaxLib.tabfilewithtitle()
    #TabT_mixarm.open(_input+'_mixarm_Xdatasat_Ymir.txt', 'mir')
    f = open(_input+'_mixarm_Xdatasat_Ymir.txt.pickletable', 'rb')
    TabT_mixarm = pickle.load(f)
    f.close()
    
    #TabT_twoarm = RaxLib.tabfilewithtitle()
    #TabT_twoarm.open(_input+'_twoarm_Xdatasat_Ymir.txt', 'mir')
    f = open(_input+'_twoarm_Xdatasat_Ymir.txt.pickletable', 'rb')
    TabT_twoarm = pickle.load(f)
    f.close()
    '''
    
    #-----------------------------------------------------------------------------------------------
    #delete mir that express less then 20% of total dataset
    tissue_TabT_mixarm = TabT_mixarm.delete('tissue')
    tissue_TabT_twoarm = TabT_twoarm.delete('tissue')
    
    N_item_list = []
    T_item_list = []
    for item in [x for x in TabT_mixarm.title_dicX.keys() if x != TabT_mixarm.main_key]:
        tissue = TabT_mixarm.read(item, 'NTpair')
        if 'N' == tissue:
            N_item_list.append(item)
        elif 'T' == tissue:
            T_item_list.append(item)
    
    
    delete_mir_list_mixarm = []
    for key in TabT_mixarm.main_key_list:
        
        N_value = []
        for N_item in N_item_list:
            read = TabT_mixarm.read(N_item, key)
            N_value.append(read)
        
        T_value = []
        for T_item in T_item_list:
            read = TabT_mixarm.read(T_item, key)
            T_value.append(read)
        
        if len([x for x in N_value if x == 1]) > len(N_item_list)*50/100 or len([x for x in T_value if x == 1]) > len(T_item_list)*50/100:
            delete_mir_list_mixarm.append(key)
            TabT_mixarm.delete(key)
        else:
            pass
    
    
    #-----------------------------------------------------------------------------------------------
    #new as mixarm and twoarm is not-paired data
    
    N_item_list = []
    T_item_list = []
    for item in [x for x in TabT_twoarm.title_dicX.keys() if x != TabT_twoarm.main_key]:
        tissue = TabT_twoarm.read(item, 'NTpair')
        if 'N' == tissue:
            N_item_list.append(item)
        elif 'T' == tissue:
            T_item_list.append(item)
    
    
    delete_mir_list_twoarm = []
    for key in TabT_twoarm.main_key_list:
        
        N_value = []
        for N_item in N_item_list:
            read = TabT_twoarm.read(N_item, key)
            N_value.append(read)
        
        T_value = []
        for T_item in T_item_list:
            read = TabT_twoarm.read(T_item, key)
            T_value.append(read)
        
        if len([x for x in N_value if x == 1]) > len(N_item_list)*50/100 or len([x for x in T_value if x == 1]) > len(T_item_list)*50/100:
            delete_mir_list_twoarm.append(key)
            TabT_twoarm.delete(key)
        else:
            pass
    
    '''
    #-----------------------------------------------------------------------------------------------
    #old as mixarm and twoarm is paired data (delete twoarm mir accroding to the mir from mixarm)
    for mir in delete_mir_list:
        
        for key in TabT_twoarm.main_key_list:
            
            if mir == key[:-3]:
                TabT_twoarm.delete(key)
    '''
    
    TabT_mixarm.insert('Y', 'tissue', tissue_TabT_mixarm, 1)
    TabT_twoarm.insert('Y', 'tissue', tissue_TabT_twoarm, 1)
    
    #-----------------------------------------------------------------------------------------------
    #open BatchID table
    os.chdir('/'.join(ori_path.split('/')[:-2]))
    Batch_table = RaxLib.tabfilewithtitle()
    Batch_table.open('aliquotReport.txt', 'Aliquot ID')
    os.chdir(ori_path)
    
    #-----------------------------------------------------------------------------------------------
    #add BatchID
    mixarm_datasets_E = RaxLib.dicKeysortVal(TabT_mixarm.title_dicX, '<-')
    mixarm_datasets_Batch_value_list = []
    for dataset in mixarm_datasets_E[1:]:
        value = Batch_table.read('BCR Batch', dataset.split('_')[0])
        mixarm_datasets_Batch_value_list.append(value)
    TabT_mixarm.insert('Y', 'BatchID', ['BatchID']+mixarm_datasets_Batch_value_list, 1)
    
    twoarm_datasets_E = RaxLib.dicKeysortVal(TabT_twoarm.title_dicX, '<-')
    twoarm_datasets_Batch_value_list = []
    for dataset in twoarm_datasets_E[1:]:
        value = Batch_table.read('BCR Batch', dataset.split('_')[0])
        twoarm_datasets_Batch_value_list.append(value)
    TabT_twoarm.insert('Y', 'BatchID', ['BatchID']+twoarm_datasets_Batch_value_list, 1)
    
    #-----------------------------------------------------------------------------------------------
    #output
    TabT_mixarm.XYreverse()
    f = open(_input+'_NT_top3_mixarm_0to1_Del.txt', 'w')
    f.write(TabT_mixarm.report())
    f.close()
    f = open(_input+'_NT_top3_mixarm_0to1_Del.txt.pickletable', 'wb')
    pickle.dump(TabT_mixarm, f, )
    f.close()
    
    TabT_twoarm.XYreverse()
    f = open(_input+'_NT_top3_twoarm_0to1_Del.txt', 'w')
    f.write(TabT_twoarm.report())
    f.close()
    f = open(_input+'_NT_top3_twoarm_0to1_Del.txt.pickletable', 'wb')
    pickle.dump(TabT_twoarm, f, )
    f.close()
    
    #-----------------------------------------------------------------------------------------------
    #output delete_mir_list
    #TabT_mixarm = RaxLib.tabfilewithtitle()
    #TabT_mixarm.open(_input+'_mixarm_Xdatasat_Ymir.txt', 'mir') 
    f = open(_input+'_mixarm_0to1_Xdatasat_Ymir.txt.pickletable', 'rb')
    TabT_mixarm = pickle.load(f)
    f.close()
    
    TabT = RaxLib.tabfilewithtitle()
    TabT.build({x:y for x ,y in [('total_mir', 1)]}, 'total_mir')
    
    for key in [x for x in TabT_mixarm.main_key_list if x != 'tissue' and x != 'NTpair']:
        TabT.append('Y', key, [key])
    
    E_ = RaxLib.multilistextend([[x for x in TabT_mixarm.main_key_list if x != 'tissue' and x != 'NTpair'], delete_mir_list_mixarm])
    TabT.append('X', 'delete_mir_mixarm', E_[1])
    
    f = open(_input+'_delete_mir_list_mixarm.txt', 'w')
    f.write(TabT.report())
    f.close()
    
    #-----------------------------------------------------------------------------------------------
    #output delete_mir_list
    #TabT_twoarm = RaxLib.tabfilewithtitle()
    #TabT_twoarm.open(_input+'_twoarm_Xdatasat_Ymir.txt', 'mir')
    f = open(_input+'_twoarm_0to1_Xdatasat_Ymir.txt.pickletable', 'rb')
    TabT_twoarm = pickle.load(f)
    f.close()
    
    TabT = RaxLib.tabfilewithtitle()
    TabT.build({x:y for x ,y in [('total_mir', 1)]}, 'total_mir')
    
    for key in [x for x in TabT_twoarm.main_key_list if x != 'tissue' and x != 'NTpair']:
        TabT.append('Y', key, [key])
    
    E_ = RaxLib.multilistextend([[x for x in TabT_twoarm.main_key_list if x != 'tissue' and x != 'NTpair'], delete_mir_list_twoarm])
    TabT.append('X', 'delete_mir_twoarm', E_[1])
    
    f = open(_input+'_delete_mir_list_twoarm.txt', 'w')
    f.write(TabT.report())
    f.close()