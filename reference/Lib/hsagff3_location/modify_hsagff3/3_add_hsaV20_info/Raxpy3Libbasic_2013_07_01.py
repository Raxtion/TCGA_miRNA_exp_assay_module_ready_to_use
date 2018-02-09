'''
history = 2013.03.28, 2013.03.29, 2013.03.31, 2013.04.02, 2013.04.16, 2013.04.17, 2013.04.18, 2013.04.24, 2013.04.26,
          2013.04.28, 2013.04.30, 2013.05.09, 2013.05.13, 2013.05.15, 2013.05.16, 2013.05.24, 2013.05.30, 2013.06.02,
          2013.06.05, 2013.06.10, 2013.06.14
'''

#-----------------------------------------------------------------------------------------------

class openFastAfile:    
    
    def __init__(self):
        self.data_number = 0
        self.FastA_dictionary = {} 
    
    def open(self, file_name = '.fa'):        
        
        FastA_dictionary = {}
        total = open(file_name,'r').read()
        A = total.split('>')
        A = multilistshorten([A])[0]
        
        for x in A:
            ID = "".join("".join(x).split('\n')[0])
            sequence = "".join([x.replace('\n', "") for x in "".join(x).split('\n')[1:]])
            FastA_dictionary[ID] = sequence
        
        self.data_number = len(FastA_dictionary.keys())
        self.FastA_dictionary = FastA_dictionary
        
        return FastA_dictionary
    
    def build(self, dictionary = {}):
        self.data_number = len(dictionary.keys())
        self.FastA_dictionary = dictionary
    
    def read(self, ID = ''):        
        
        result = self.FastA_dictionary[ID]
        
        return result
    
    def report(self, ID_list = []):
        
        if ID_list == []:
            ID_list = list(self.FastA_dictionary.keys())
        else:
            ID_list = ID_list
        
        ID_list.sort()
        E = []
        for x in ID_list:
            E.append('>' + x + '\n' + self.read(x))
        
        display = '\n'.join(E)+'\n'+str(len(E))
        
        return display

#-----------------------------------------------------------------------------------------------

def listXYreverse(XY = [[]]):
    
    i = 1
    j = 1
    E = []
    YX = []
    while i < len(XY[0]) + 1:
        
        j = 1
        E = []
        while j < len(XY) + 1:
            E.append(XY[j-1][i-1])
            j = j + 1
            
        YX.append(E)
      
        i = i + 1

    return YX

#-----------------------------------------------------------------------------------------------

def listNSplit(E = [], N = 0):
    
    result = [E[i:i+N] for i in range(0, len(E), N)]

    return result

#-----------------------------------------------------------------------------------------------

def listPartSeparate(E = [], Part = 0):
    
    N = int((len(E) / Part)) + 1
    result = [E[i:i+N] for i in range(0, len(E), N)]
    
    return result

#-----------------------------------------------------------------------------------------------

def multilistextend(inp = [[]]):
    
    oup = []
    E = [int(len(x)) for x in inp]
    E.sort()
    
    for x in inp:
        if int(len(x)) < E[-1]: 
            x.extend(' '*(E[-1]-int(len(x))))
            oup.append(x)
        else:
            oup.append(x)

    return oup

#-----------------------------------------------------------------------------------------------

def multilistshorten(inp = [[]]):
    
    oup = []
    for E in inp:
        oup.append([x for x in E if x != ' ' and x != ""])
    
    return oup

#-----------------------------------------------------------------------------------------------

def foldchange(A_data_read = 0.0, B_data_read = 0.0):
    
    fold_change_between_A_B, direction, fold = '', '', ''
    if float(A_data_read) > float(B_data_read):
        fold_change_between_A_B = '+'+str(float(A_data_read)/float(B_data_read))
        direction = '+'
        fold = str(float(A_data_read)/float(B_data_read))
        
    elif float(A_data_read) < float(B_data_read):
        fold_change_between_A_B = '-'+str(float(B_data_read)/float(A_data_read))
        direction = '-'
        fold = str(float(B_data_read)/float(A_data_read))
    
    return fold_change_between_A_B, direction, fold

#-----------------------------------------------------------------------------------------------

def armexchange(A_5p = 0.0, A_3p = 0.0, B_5p = 0.0, B_3p = 0.0):
    critical_value = 0.0
    
    a = float(A_5p)/(float(A_5p) + float(A_3p))
    b = float(A_3p)/(float(A_5p) + float(A_3p))
    c = float(B_5p)/(float(B_5p) + float(B_3p))
    d = float(B_3p)/(float(B_5p) + float(B_3p))
    
    e = a-b
    f = c-d
    if f != 0.0:
        critical_value = e/f
    else:
        critical_value = 'Denominator = 0'
    
    return critical_value

#-----------------------------------------------------------------------------------------------

def meanSD(E = [0.0]):
    
    mean = 0.0
    SD = 0.0
    i = 0.0
    j = 0.0
    
    if len(E) == 0:
        
        mean = 'na'
        SD = 'na'
        
    elif len(E) > 0:
        
        for x in E:
            i = i + float(x)
        mean = float(i/float(len(E)))
        
        for x in E:
            j = j + (float(x)-mean)**2
        SD = float(j/float(len(E)))**(0.5)
    
    return mean, SD

#-----------------------------------------------------------------------------------------------

def listIntersection(E = [[]]):
    
    if len(E) < 2:
        print('input = ')
        print(E)
        print('input is too small to ListIntersection')
        quit()
    elif len(E) == 2:
        R = set(E[0]) & set(E[1])
    else:
        R = set(E[0]) & set(E[1])
        i = 3
        while i < len(E) + 1:
            R = set(R) & set(E[i-1])
            i = i + 1
    
    T = list(R)
    T.sort()
    
    return T
    
#-----------------------------------------------------------------------------------------------

def listUnion(E = [[]]):
    
    if len(E) < 2:
        print('input = ')
        print(E)
        print('input is too small to ListIntersection')
        quit()
    elif len(E) == 2:
        R = set(E[0]) | set(E[1])
    else:
        R = set(E[0]) | set(E[1])
        i = 3
        while i < len(E) + 1:
            R = set(R) | set(E[i-1])
            i = i + 1
    
    T = list(R)
    T.sort()
    
    return T

#-----------------------------------------------------------------------------------------------

class tabfilewithtitle:
    
    def __init__(self):
        self.title_tuple = [()]
        self.title_dicX = {}
        self.title_dicY = {}
        self.title_box = []
        self.main_key = ""
        self.main_key_list = []
        self.len_X = 0
        self.len_Y = 0
    
    def build(self, title_dicX = {}, main_key = ''):
        self.main_key = main_key
        self.title_dicX = title_dicX
        self.title_dicY = {}
        title_tuple = [(x[1],x[0]) for x in title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.len_X = len(title_dicX.keys())
        self.len_Y = 0
        self.title_box = []
        
        return self.title_tuple

    def open(self, file_name = '.txt', main_key = ''):
        
        B = open(file_name,'r').read()
        title = B.split('\n')[0]
        data = B.split('\n')[1:-1]
        data = multilistshorten([data])[0]
        title_list = title.split('\t')
        
        k = 1
        title_dic = {}
        for title in title_list:
            title_dic[title] = k
            k = k + 1
        
        j = 1
        key_M = {}
        for x in data:
            W = "".join(x).split('\t')[title_dic[main_key]-1]
            key_M[W] = j
            j = j + 1
        
        title_box = []
        position = 1
        for title in title_list:
            
            E = []            
            for key in data:
                
                integer = "".join("".join(key).split('\t')[position-1].split('.')[0])                
                decimal = "".join("".join(key).split('\t')[position-1].split('.')[-1])
                if str.isdigit(integer) == True and str.isdigit(decimal) == True:
                    E.append(float("".join(key).split('\t')[position-1]))
                else:
                    E.append("".join(key).split('\t')[position-1])
                
            title_box.append(E) 
            position = position + 1
        
        #title_box.insert(0, TCGA_M.items())
        item_number = len(title_dic.keys())
        
        title_tuple = [(x[1],x[0]) for x in title_dic.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.title_dicX = title_dic
        self.title_dicY = key_M
        self.title_box = title_box
        self.len_X = item_number
        self.len_Y = len(data)
        self.main_key = main_key
        self.main_key_list = self.title_box[self.title_dicX[self.main_key]-1]
        
        return self.title_tuple
    
    def read(self, X = '', Y = ''):
        
        P = (self.title_dicX, self.title_dicY, self.title_box)
        
        result = P[2][P[0][X]-1][P[1][Y]-1]
        
        return result
    
    def readrow(self, XY = 'X or Y', inp = ""):
        
        P = (self.title_dicX, self.title_dicY, self.title_box)
        
        if XY == 'X':
            
            result = P[2][P[0][inp]-1]
            
        elif XY == 'Y':
            
            E = listXYreverse(self.title_box)
            result = E[P[1][inp]-1]
            
        else:
            print(inp)
            print('Error at tabfilewithtitle.readrow()\nPlease define the XY value ! ')
            quit()
        
        return result
    
    
    def printtable(self, column_list = []):
        
        title_output_list = []
        value_output_list = []
        
        if column_list != []:
            
            dic = {}
            for column in column_list:
                dic[column] = self.title_dicX[column]
            
            title_output_list = dicKeysortVal(dic,'<-')
        else:
            title_output_list = dicKeysortVal(self.title_dicX,'<-')
        
        data = []
        for x in title_output_list:
            data.append([])
        
        i = 1        
        for title in title_output_list:            
            for key in self.main_key_list:
                data[i-1].append(str(self.read(title, key)))
                
            i = i + 1
        
        oup = listXYreverse(data)
        oup.insert(0,title_output_list)
        copydata = listXYreverse(oup)
        
        lenght_list = []
        for x in copydata:
            E = []
            for y in x:
                E.append(len(y))
            E.sort()
            lenght_list.append(E[-1])
        
        for x in oup:
            
            E = []
            i = 1
            for y in x:
                
                W = ' '*(lenght_list[i-1]-len(y))+y
                E.append('| '+W+' |')
                
                i = i + 1
            value_output_list.append("".join(E))
        print('\n'.join(value_output_list))
        
        return
    
    def report(self):
        
        display_list = []
        
        title_list = [x[1] for x in self.title_tuple]
        one_line = '\t'.join(title_list)
        display_list.append(one_line)
        
        data = []
        for key in self.main_key_list:
            E = []
            for title in title_list:
                result = self.read(title, key)
                E.append(str(result))
            data.append('\t'.join(E))
        two_line = '\n'.join(data)
        display_list.append(two_line)
        
        display = '\n'.join(display_list)+'\n'+str(len(data))
        
        return display
    
    def append(self, XY = 'X or Y', inp = "", value = []):
        
        if XY == 'X':
            
            if len(value) == self.len_Y:
                self.title_dicX[inp] = self.len_X+1
                self.title_box.append(value)
            else:
                print(inp)
                print(value)
                print('Error at tabfilewithtitle.append()\nPlease check the length of value !\nThat is different from the table. ')
                quit()
        
        elif XY == 'Y':
            
            if len(value) == self.len_X:
                
                if self.title_box != []:
                    self.title_dicY[inp] = self.len_Y+1
                    self.title_box = listXYreverse(self.title_box)
                    self.title_box.append(value)
                    self.title_box = listXYreverse(self.title_box)
                elif self.title_box == []:
                    self.title_dicY[inp] = self.len_Y+1
                    self.title_box.append(value)
                    self.title_box = listXYreverse(self.title_box)
            else:
                print(inp)
                print(value)
                print('Error at tabfilewithtitle.append()\nPlease check the length of value !\nThat is different from the table. ')
                quit()
        
        title_tuple = [(x[1],x[0]) for x in self.title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.len_X = len(self.title_dicX)
        self.len_Y = len(self.title_dicY)
        self.main_key_list = self.title_box[self.title_dicX[self.main_key]-1]
        
        return
    
    def delete(self, inp = ""):
        
        if inp in self.title_dicX:
            
            value_large_Then_inp = [x for x in self.title_dicX.keys() if self.title_dicX[x] > self.title_dicX[inp]]
            for key in value_large_Then_inp:
                self.title_dicX[key] = self.title_dicX[key]-1
            
            copy = self.title_box[self.title_dicX[inp]-1]
            del self.title_box[self.title_dicX[inp]-1]
            
            if self.title_box != []:
                pass
            elif self.title_box == []:
                pass
                #print('title_box is empty !!')
            
            delet = self.title_dicX.pop(inp)
            
        elif inp in self.title_dicY:
            
            value_large_Then_key = [x for x in self.title_dicY.keys() if self.title_dicY[x] > self.title_dicY[inp]]
            for key in value_large_Then_key:
                self.title_dicY[key] = self.title_dicY[key]-1
            
            self.title_box = listXYreverse(self.title_box)
            copy = self.title_box[self.title_dicY[inp]-1]
            del self.title_box[self.title_dicY[inp]-1]
            
            if self.title_box != []:
                self.title_box = listXYreverse(self.title_box)
            elif self.title_box == []:
                pass
                #print('title_box is empty !!')
            
            delet = self.title_dicY.pop(inp)
            
        else:
            print(inp)
            print('Error at tabfilewithtitle.delete()\nPlease check the input !\nThat is not in the table. ')
            quit()
        
        title_tuple = [(x[1],x[0]) for x in self.title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.len_X = len(self.title_dicX)
        self.len_Y = len(self.title_dicY)
        if len(self.title_dicY) != 0:
            self.main_key_list = self.title_box[self.title_dicX[self.main_key]-1]
        else:
            self.main_key_list = []
        
        return copy
    
    def insert(self, XY = 'X or Y', inp = "", value = [], location = 0):
        
        if XY == 'X':
            
            if len(value) == self.len_Y:
                
                value_large_Then_inp = [x for x in self.title_dicX.keys() if self.title_dicX[x] > location]
                for key in value_large_Then_inp:
                    self.title_dicX[key] = self.title_dicX[key]+1
                
                self.title_box.insert(location, value)
                
                self.title_dicX[inp] = location+1
            else:
                print(inp)
                print(value)
                print('Error at tabfilewithtitle.insert()\nPlease check the length of value !\nThat is different from the table. ')
                quit()
        
        elif XY == 'Y':
            
            if len(value) == self.len_X:
                
                value_large_Then_key = [x for x in self.title_dicY.keys() if self.title_dicY[x] > location]
                for key in value_large_Then_key:
                    self.title_dicY[key] = self.title_dicY[key]+1
                
                self.title_box = listXYreverse(self.title_box)
                self.title_box.insert(location, value)
                self.title_box = listXYreverse(self.title_box)
                
                self.title_dicY[inp] = location+1
            else:
                print(inp)
                print(value)
                print('Error at tabfilewithtitle.insert()\nPlease check the length of value !\nThat is different from the table. ')
                quit()
        
        title_tuple = [(x[1],x[0]) for x in self.title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.len_X = len(self.title_dicX)
        self.len_Y = len(self.title_dicY)
        self.main_key_list = self.title_box[self.title_dicX[self.main_key]-1]
        
        return
    
    def find(self, value = '', column = ''):
        
        result_Array = []
        if column == '':
            sort_title_list = dicKeysortVal(self.title_dicX, '<-')
            for title in sort_title_list:
                E = [x for x in self.title_dicY.keys() if self.read(title, x) == value]
                if E != []:
                    for key in E:
                        result_Array.append((title, key))
            
            if result_Array == []:
                result_Array.append(('non', 'non'))
        
        else:
            title = column
            E = [x for x in self.title_dicY.keys() if self.read(title, x) == value]
            if E != []:
                for key in E:
                    result_Array.append((title, key))
            
            if result_Array == []:
                result_Array.append(('non', 'non'))
        
        return result_Array
    
#-----------------------------------------------------------------------------------------------

def newOrderfile(title_list = ['order', 'A', 'B' ], group_list = ['a', 'b' ], outputFile_Name = '.txt'):
    
    table = tabfilewithtitle()
    
    i = 1
    title = {}
    for tit in title_list:
        title[tit] = i
        
        i = i + 1
    table.build(title, 'order')
    
    i =  1
    group_list.sort()
    for group in group_list:
        order = str(i)+'_'
        table.append('Y', order, [order]+group)
        
        i = i + 1
    
    if outputFile_Name != '.txt':
        f = open(outputFile_Name, 'w')
        f.write(table.report())
        f.close()
    else:
        pass
    
    return table

#-----------------------------------------------------------------------------------------------

def tabfileaddorder(file_name = '.txt', title_list = ['order', 'A', 'B'], outputFile_Name = '.txt'):
    
    B = open(file_name,'r').read()
    Line_list = B.split('\n')[0:-1]
    
    E = []
    i = 1
    for Line in Line_list:
        E.append('\t'.join([str(i)+'_', Line]))
        i = i + 1
    
    result = '\n'.join(E)
    title = '\t'.join(title_list)
    display = '\n'.join([title, result, str(len(E))])
    
    if outputFile_Name != '.txt':
        f = open(outputFile_Name, 'w')
        f.write(display)
        f.close()
    else:
        pass
    
    return display

#-----------------------------------------------------------------------------------------------

def dicKeysortVal(M = 'dictionary', pattern = '-> or <-'):
    result = []
    BtupleA = []
    if pattern == '->':
        
        for x in M.keys():
            BtupleA.append((M[x],x))
        
        BtupleA.sort()
        BtupleA.reverse()
        for x in BtupleA:
            result.append(x[1])
        
    elif pattern == '<-':
        
        for x in M.keys():
            BtupleA.append((M[x],x))
        
        BtupleA.sort()        
        for x in BtupleA:
            result.append(x[1])
        
    return result

#-----------------------------------------------------------------------------------------------

def ntComRev(string = '', complete = 'C or nC', reverse = 'r or nr', case = 'U or L', seq_type = 'RNA or DNA'):
    
    if string.islower() == True:
        string_type = 1
    elif string.isupper() == True:
        string_type = 0
    else:
        #print('string is hybrid! ')
        string_type = 2
    
    if complete == 'C':
        string = string.upper()
        E = []
        for x in string:
            if x == 'G':
                x = 'C'
                E.append(x)
            elif x == 'C':
                x = 'G'
                E.append(x)
            elif x == 'T':
                x = 'A'
                E.append(x)
            elif x == 'A':
                x = 'T'
                E.append(x)
            elif x == 'U':
                x = 'A'
                E.append(x)
    elif complete == 'nC':
        string = string.upper()
        E = [x for x in string]
    else:
        string = string.upper()
        E = [x for x in string]
    
    if reverse == 'nr':
        result = "".join(E)
    elif reverse == 'r':
        E.reverse()
        result = "".join(E)
    else:
        result = "".join(E)
    
    if case == 'U':
        result = result.upper()
    elif case == 'L':
        result = result.lower()
    else:
        if string_type == 1:
            result = result.lower()
        elif string_type == 0:
            result = result.upper()
        elif string_type == 2:
            result = result.upper()
    
    if seq_type == 'RNA':
        result = result.replace('T', 'U')
    elif seq_type == 'DNA':
        result = result.replace('U', 'T')
    else:
        result = result
    
    return result

#-----------------------------------------------------------------------------------------------

