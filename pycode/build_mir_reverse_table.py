import os, math, pickle
#import Raxpy3Libbasic_2013_09_14, Raxpy3LibStatistic_2013_08_04

ori_path = os.getcwd()

#-----------------------------------------------------------------------------------------------
#parameter
#RaxLib = Raxpy3Libbasic_2013_09_14
#RaxStat = Raxpy3LibStatistic_2013_08_04
input_file = 'TCGASTAD_NT_top3_twoarm_nomodify_0to001_Del.txt'

outputfile_Name = 'NGS_critical_value_table.txt'

#-----------------------------------------------------------------------------------------------

def running(ori_path = '', _input = 'TCGASTAD_NT_top3_twoarm_nomodify_0to001_Del.txt', RaxLib = 'Raxpy3Libbasic_2013_09_14', RaxStat = 'Raxpy3LibStatistic_2013_08_04'):
    
    #-----------------------------------------------------------------------------------------------
    #open file
    #TabTori = RaxLib.tabfilewithtitle()
    #TabTori.open(input_file, 'mir')
    #TabTori.XYreverse()
    
    #f = open(_input, 'rb')
    #TabTori = pickle.load(f)
    #f.close()
    
    
    TabTori = RaxLib.tabfilewithtitle()
    TabTori.open(_input+'_NT_top3_twoarm_0to0.txt', 'mir')
    f = open(_input+'_NT_top3_twoarm_0to0.txt.pickletable', 'wb')
    pickle.dump(TabTori, f, )
    f.close()
    
    TabTori.XYreverse()
    f = open(_input+'_twoarm_0to0_Xdatasat_Ymir.txt', 'w')
    f.write(TabTori.report())
    f.close()
    f = open(_input+'_twoarm_0to0_Xdatasat_Ymir.txt.pickletable', 'wb')
    pickle.dump(TabTori, f, )
    f.close()
    
    '''
    NT_list = TabTori.delete('tissue')
    
    for item in [x for x in TabTori.title_dicX.keys() if x != TabTori.main_key]:
        for key in TabTori.main_key_list:
            read = TabTori.read(item, key)
            TabTori.title_box[TabTori.title_dicX[item]-1][TabTori.main_key_list.index(key)] = math.log(read, 2)
    
    TabTori.insert('Y', 'tissue', NT_list)
    
    f = open('Look.txt.pickletable', 'wb')
    pickle.dump(TabTori, f, )
    #f.write(TabTori.report())
    f.close()
    
    #-----------------------------------------------------------------------------------------------
    #open file
    #TabTori = RaxLib.tabfilewithtitle()
    #TabTori.open('Look.txt', 'mir')
    f = open('Look.txt.pickletable', 'rb')
    TabTori = pickle.load(f)
    f.close()
    #NT_list = TabTori.delete('tissue')
    '''
    
    #-----------------------------------------------------------------------------------------------
    #new TabToup
    TabToup = RaxLib.tabfilewithtitle()
    TabToup.build({x:y for x, y in [('mir', 1)]}, 'mir')
    
    new_main_key = [x[:-3] for x in TabTori.main_key_list if x != 'tissue' and x != 'NTpair']
    new_main_key = list(set(new_main_key))
    new_main_key.sort()
    for key in new_main_key:
        TabToup.append('Y', key, [key])
    
    print(len(TabToup.main_key_list))
    #TabToup.printtable()
    
    #-----------------------------------------------------------------------------------------------
    #
    N_item_list = []
    T_item_list = []
    for item in [x for x in TabTori.title_dicX.keys() if x != TabTori.main_key]:
        tissue = TabTori.read(item, 'NTpair')
        if 'N' == tissue:
            N_item_list.append(item)
        elif 'T' == tissue:
            T_item_list.append(item)
    
    factor_list = []
    critical_value_list = []
    row_CV_list = []
    N_5p_mean_list = []
    N_3p_mean_list = []
    T_5p_mean_list = []
    T_3p_mean_list = []
    i = 1
    '''
    print(TabToup.main_key_list[-1])
    print(len(TabToup.main_key_list))
    print(len([x for x in TabTori.main_key_list if x != 'tissue']))
    '''
    for oup_key in TabToup.main_key_list:
        key_3p = [x for x in TabTori.main_key_list if x != 'tissue' and x != 'NTpair'][i-1]
        key_5p = [x for x in TabTori.main_key_list if x != 'tissue' and x != 'NTpair'][i+1-1]
        '''
        print(oup_key)
        print(key_3p)
        print(key_5p)
        '''
        N_3p_value = []
        for N_item in N_item_list:
            read = TabTori.read(N_item, key_3p)
            N_3p_value.append(read)
        N_3p_mean = RaxStat.meanSD(N_3p_value)[0]
        N_3p_mean_list.append(N_3p_mean)
        
        T_3p_value = []
        for T_item in T_item_list:
            read = TabTori.read(T_item, key_3p)
            T_3p_value.append(read)
        T_3p_mean = RaxStat.meanSD(T_3p_value)[0]
        T_3p_mean_list.append(T_3p_mean)
        
        N_5p_value = []
        for N_item in N_item_list:
            read = TabTori.read(N_item, key_5p)
            N_5p_value.append(read)
        N_5p_mean = RaxStat.meanSD(N_5p_value)[0]
        N_5p_mean_list.append(N_5p_mean)
        
        T_5p_value = []
        for T_item in T_item_list:
            read = TabTori.read(T_item, key_5p)
            T_5p_value.append(read)
        T_5p_mean = RaxStat.meanSD(T_5p_value)[0]
        T_5p_mean_list.append(T_5p_mean)
        
        
        critical_value = RaxStat.armexchange(T_5p_mean, T_3p_mean, N_5p_mean, N_3p_mean)
        row_CV_list.append(critical_value)
        
        if type(critical_value) == type(''):
            factor = 'non'
            critical_value_list.append(critical_value)
            factor_list.append(factor)
        elif critical_value > 0.0:
            factor = '+'
            critical_value_list.append(math.fabs(critical_value))
            factor_list.append(factor)
        elif critical_value < 0.0:
            factor = '-'
            critical_value_list.append(math.fabs(critical_value))
            factor_list.append(factor)
        else:
            if '-' in str(critical_value):
                critical_value_list.append(math.fabs(critical_value))
                factor_list.append('-')
            else:
                critical_value_list.append(math.fabs(critical_value))
                factor_list.append('+')
        
        i = i + 2
    print(len(factor_list))
    TabToup.append('X', 'factor', factor_list)
    TabToup.append('X', 'critical_value', critical_value_list)
    TabToup.append('X', 'CV', row_CV_list)
    TabToup.append('X', 'T_5p_mean', T_5p_mean_list)
    TabToup.append('X', 'T_3p_mean', T_3p_mean_list)
    TabToup.append('X', 'N_5p_mean', N_5p_mean_list)
    TabToup.append('X', 'N_3p_mean', N_3p_mean_list)
    
    
    #-----------------------------------------------------------------------------------------------
    #output
    
    f = open(outputfile_Name, 'w')
    f.write(TabToup.report())
    f.close()



