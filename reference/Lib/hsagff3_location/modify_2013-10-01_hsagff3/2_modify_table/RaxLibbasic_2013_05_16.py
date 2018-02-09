'''
history = 2013.03.28, 2013.03.29, 2013.03.31, 2013.04.02, 2013.04.16, 2013.04.17, 2013.04.18, 2013.04.24, 2013.04.26,
          2013.04.28, 2013.04.30, 2013.05.09, 2013.05.13, 2013.05.15
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
    
    def read(self, ID = ''):        
        
        try:
            index = self.FastA_box[0].index(ID)
            result = self.FastA_box[-1][self.FastA_box[0].index(ID)]
        except:
            result = self.FastA_dictionary[ID]
        
        return result
    
    def rebuild(self, ID_list = []):
        
        ID_list.sort()
        
        E = []        
        for x in ID_list:
            
            i = 1            
            read_list = []
            while i < self.column_number+1-1:
                read_list.append(self.FastA_box[i-1][self.FastA_box[0].index(x)])                
                i = i + 1        
                
            E.append('>' + ' '.join(read_list) + '\n' + self.FastA_box[-1][self.FastA_box[0].index(x)])
        
        return E

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
        print 'input = '
        print E
        print 'input is too small to ListIntersection'
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
        print 'input = '
        print E
        print 'input is too small to ListIntersection'
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
        
        return

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
        
        return self.title_tuple
    
    def read(self, X = '', Y = ''):
        
        P = [self.title_dicX, self.title_dicY, self.title_box]
        
        result = P[2][P[0][X]-1][P[1][Y]-1]
        
        return result
    
    def printtable(self):
        
        title_output_list = []
        value_output_list = []        
        title_output_list = dicKeysortVal(self.title_dicX,'<-')
        
        data = []
        for x in title_output_list:
            data.append([])
        
        i = 1
        main_key_list = self.title_box[self.title_dicX[self.main_key]-1]        
        for title in title_output_list:            
            for key in main_key_list:
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
        print '\n'.join(value_output_list)
        
        return
    
    def report(self):
        
        display_list = []
        
        title_list = [x[1] for x in self.title_tuple]
        one_line = '\t'.join(title_list)
        display_list.append(one_line)
        
        main_key_list = self.title_box[self.title_dicX[self.main_key]-1]
        data = []
        for key in main_key_list:
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
                print inp
                print value
                print 'Error at tabfilewithtitle.append()\nPlease check the length of value !\nThat is different from the table. '
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
                print inp
                print value
                print 'Error at tabfilewithtitle.append()\nPlease check the length of value !\nThat is different from the table. '
                quit()
        
        title_tuple = [(x[1],x[0]) for x in self.title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.len_X = len(self.title_dicX)
        self.len_Y = len(self.title_dicY)
        
        return
    
    def delete(self, inp = ""):
        
        if self.title_dicX.has_key(inp) == True:
            
            value_large_Then_inp = [x for x in self.title_dicX.keys() if self.title_dicX[x] > self.title_dicX[inp]]
            for key in value_large_Then_inp:
                self.title_dicX[key] = self.title_dicX[key]-1
            
            copy = self.title_box[self.title_dicX[inp]-1]
            del self.title_box[self.title_dicX[inp]-1]
            
            if self.title_box != []:
                pass
            elif self.title_box == []:
                pass
                #print 'title_box is empty !!'
            
            delet = self.title_dicX.pop(inp)
            
        elif self.title_dicY.has_key(inp) == True:
            
            
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
                #print 'title_box is empty !!'
            
            delet = self.title_dicY.pop(inp)
            
        else:
            print inp
            print 'Error at tabfilewithtitle.delete()\nPlease check the input !\nThat is not in the table. '
            quit()
        
        title_tuple = [(x[1],x[0]) for x in self.title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.len_X = len(self.title_dicX)
        self.len_Y = len(self.title_dicY)
        
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
                print inp
                print value
                print 'Error at tabfilewithtitle.insert()\nPlease check the length of value !\nThat is different from the table. '
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
                print inp
                print value
                print 'Error at tabfilewithtitle.insert()\nPlease check the length of value !\nThat is different from the table. '
                quit()
        
        title_tuple = [(x[1],x[0]) for x in self.title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.len_X = len(self.title_dicX)
        self.len_Y = len(self.title_dicY)
        
        return
    
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

def ntComRev(string = '', case = 'Ucase or Lcase', reverse = 'n'):
    
    string = string.upper()
    E = []
    for x in string:
        if x == 'G':
            x = x.replace('G', 'C')
            E.append(x)
        elif x == 'C':
            x = x.replace('C', 'G')
            E.append(x)
        elif x == 'T':
            x = x.replace('T', 'A')
            E.append(x)
        elif x == 'A':
            x = x.replace('A', 'T')
            E.append(x)
        elif x == 'U':
            x = x.replace('U', 'A')
            E.append(x)
    
    if reverse == 'n':
        result = "".join(E)
    elif reverse == 'r':
        E.reverse()
        result = "".join(E)
    
    if case == 'Ucase':
        result = result.upper()
    elif case == 'Lcase':
        result = result.lower()
    
    return result

#-----------------------------------------------------------------------------------------------