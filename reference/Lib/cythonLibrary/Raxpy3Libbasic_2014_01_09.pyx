import hashlib
import pickle
import os
import time
import copy
ori_path = os.getcwd()
'''
history = 2013.03.28, 2013.03.29, 2013.03.31, 2013.04.02, 2013.04.16, 2013.04.17, 2013.04.18, 2013.04.24, 2013.04.26,
          2013.04.28, 2013.04.30, 2013.05.09, 2013.05.13, 2013.05.15, 2013.05.16, 2013.05.24, 2013.05.30, 2013.06.02,
          2013.06.05, 2013.06.10, 2013.06.14, 2013.07.01, 2013.07.09, 2013.07.10, 2013.07.30, 2013.08.02, 2013.09.06,
          2013.09.11, 2013.09.14, 2013.10.02, 2013.11.10, 2013.11.12, 2013.11.27, 2013.11.30, 2013.12.02, 2013.12.10,
          2013.12.31
'''

#-----------------------------------------------------------------------------------------------

class openFastAfile:    
    
    def __init__(self):
        self.data_number = 0
        self.FastA_dictionary = {} 
    
    def open(self, file_name = '.fa'):        
        
        FastA_dictionary = {}
        total = open(file_name,'r').read()
        total = '\n'.join(total.split('\n')[0:-1])
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

def listIntersection(E = [[]]):
    
    if len(E) < 2:
        print('input = ')
        print(E)
        print('input is too small to ListIntersection')
        print('print location: ')
        input()
        raise KeyboardInterrupt
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
        print('print location: ')
        input()
        raise KeyboardInterrupt
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
        self.title_list = []
        self.len_X = 0
        self.len_Y = 0
    
    def copy(self):
        return copy.deepcopy(self)
    
    def build(self, XY = 'X or Y', title_dic = {}, main_key = ''):
        
        if XY == 'X':
            self.main_key = main_key
            self.title_dicX = title_dic
            self.title_dicY = {}
            title_tuple = [(x[1],x[0]) for x in title_dic.items()]
            title_tuple.sort()
            self.title_tuple = title_tuple
            self.len_X = len(title_dic.keys())
            self.len_Y = 0
            self.title_box = []
            self.main_key_list = []
            self.title_list = dicKeysortVal({x:y for x, y in title_dic.items() if x != main_key}, '<-')
        elif XY == 'Y':
            self.main_key = main_key
            self.title_dicX = {main_key:1}
            self.title_dicY = title_dic
            self.title_tuple = [(1, main_key)]
            main_key_list = dicKeysortVal(title_dic, '<-')
            self.main_key_list = main_key_list
            self.title_list = []
            self.title_box = [main_key_list]
            self.len_Y = len(self.title_dicY)
            self.len_X = len(self.title_dicX)
        else:
            print('title_dic length = ', len(title_dic))
            print('main_key = ', main_key)
            print('Error at tabfilewithtitle.build()\nPlease define the XY value ! ')
            print('print location: ')
            input()
            raise KeyboardInterrupt
        
        return self.title_tuple

    def open(self, file_name = '.txt', main_key = ''):
        
        target_path = os.getcwd()
        
        h = hashlib.new('ripemd160')
        f = open(file_name, 'rb')
        h.update(f.read())
        f.close()
        
        if os.path.exists(target_path+r'/__pycache__') == False:
            os.makedirs(target_path+r'/__pycache__')
        os.chdir(target_path+r'/__pycache__')
        
        try:
            for i in range(3):
                try:
                    Pickle_box = open('tabfilewithtitle_cache', 'rb')
                    file_box = pickle.load(Pickle_box)
                    Pickle_box.close()
                    
                    if file_name in file_box[0].keys():
                        
                        if h.hexdigest() == file_box[0][file_name]:
                            
                            if file_box[1][file_name] == main_key:
                                
                                f = open(file_box[2][file_name], 'rb')
                                pickle_EE = pickle.load(f)
                                f.close()
                                self.title_tuple = pickle_EE[0]
                                self.title_dicX = pickle_EE[1]
                                self.title_dicY = pickle_EE[2]
                                self.title_box = pickle_EE[3]
                                self.main_key = pickle_EE[4]
                                self.main_key_list = pickle_EE[5]
                                self.title_list = pickle_EE[6]
                                self.len_X = pickle_EE[7]
                                self.len_Y = pickle_EE[8]
                                
                                Pickle_box.close()
                                os.chdir(target_path)
                                return self.title_tuple
                                
                            else:
                                file_box[0].pop(file_name)
                                file_box[1].pop(file_name)
                                file_box[2].pop(file_name)
                                pass
                            
                        else:
                            file_box[0].pop(file_name)
                            file_box[1].pop(file_name)
                            file_box[2].pop(file_name)
                            pass
                        
                    else:
                        pass
                    
                except:
                    time.sleep(1)
            
        except:
            file_box = [{}, {}, {}]
        
        os.chdir(target_path)
        
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
        self.title_list = dicKeysortVal({x:y for x, y in title_dic.items() if x != main_key}, '<-')
        
        #os.chdir(ori_path+r'/__pycache__')
        os.chdir(target_path+r'/__pycache__')
        f = open(file_name+'.tabfilewithtitle.picklebox', 'wb')
        pickle.dump([self.title_tuple, self.title_dicX, self.title_dicY, self.title_box,
                     self.main_key, self.main_key_list, self.title_list, self.len_X, self.len_Y], f, )
        f.close()
        
        file_box[0][file_name] = h.hexdigest()
        file_box[1][file_name] = self.main_key
        file_box[2][file_name] = file_name+'.tabfilewithtitle.picklebox'
        
        Pickle_box = open('tabfilewithtitle_cache', 'wb')
        pickle.dump(file_box, Pickle_box, )
        Pickle_box.close()
        
        os.chdir(target_path)
        return self.title_tuple
    
    def read(self, X = '', Y = ''):
        
        P = (self.title_dicX, self.title_dicY, self.title_box)
        try:
            result = P[2][P[0][X]-1][P[1][Y]-1]
        except:
            print('X =', X, ', Y =', Y)
            print("Can't catch the value !")
            input()
            raise KeyboardInterrupt
        
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
            print('print location: ')
            input()
            raise KeyboardInterrupt
        
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
    
    def report(self, file_name = '.txt'):
        
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
        
        if file_name != '.txt':
            
            file_box = [{}, {}, {}]
            
            h = hashlib.new('ripemd160')
            h.update(display.encode('cp950'))
            
            if os.path.exists(ori_path+r'/__pycache__') == False:
                os.makedirs(ori_path+r'/__pycache__')
            os.chdir(ori_path+r'/__pycache__')
            
            f = open(file_name+'.tabfilewithtitle.picklebox', 'wb')
            pickle.dump([self.title_tuple, self.title_dicX, self.title_dicY, self.title_box,
                         self.main_key, self.main_key_list, self.title_list, self.len_X, self.len_Y], f, )
            f.close()
            
            file_box[0][file_name] = h.hexdigest()
            file_box[1][file_name] = self.main_key
            file_box[2][file_name] = file_name+'.tabfilewithtitle.picklebox'
            
            Pickle_box = open('tabfilewithtitle_cache', 'wb')
            pickle.dump(file_box, Pickle_box, )
            Pickle_box.close()
        else:
            pass
        
        os.chdir(ori_path)
        
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
                print('print location: ')
                input()
                raise KeyboardInterrupt
        
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
                print('print location: ')
                input()
                raise KeyboardInterrupt
        
        else:
            print(inp)
            print(value)
            print('Error at tabfilewithtitle.append()\nPlease define the XY value ! ')
            print('print location: ')
            input()
            raise KeyboardInterrupt
        
        title_tuple = [(x[1],x[0]) for x in self.title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.len_X = len(self.title_dicX)
        self.len_Y = len(self.title_dicY)
        self.main_key_list = self.title_box[self.title_dicX[self.main_key]-1]
        self.title_list = dicKeysortVal({x:y for x, y in self.title_dicX.items() if x != self.main_key}, '<-')
        
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
            print('print location: ')
            input()
            raise KeyboardInterrupt
        
        title_tuple = [(x[1],x[0]) for x in self.title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.len_X = len(self.title_dicX)
        self.len_Y = len(self.title_dicY)
        
        if len(self.title_dicY) != 0:
            self.main_key_list = self.title_box[self.title_dicX[self.main_key]-1]
        else:
            self.main_key_list = []
        
        if len(self.title_dicX) != 0:
            self.title_list = dicKeysortVal({x:y for x, y in self.title_dicX.items() if x != self.main_key}, '<-')
        else:
            self.title_list = []
        
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
                print('print location: ')
                input()
                raise KeyboardInterrupt
        
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
                print('print location: ')
                input()
                raise KeyboardInterrupt
        
        else:
            print(inp)
            print(value)
            print('Error at tabfilewithtitle.insert()\nPlease define the XY value ! ')
            print('print location: ')
            input()
            raise KeyboardInterrupt
        
        title_tuple = [(x[1],x[0]) for x in self.title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.len_X = len(self.title_dicX)
        self.len_Y = len(self.title_dicY)
        self.main_key_list = self.title_box[self.title_dicX[self.main_key]-1]
        self.title_list = dicKeysortVal({x:y for x, y in self.title_dicX.items() if x != self.main_key}, '<-')
        
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
    
    def XYreverse(self):
        
        if self.title_dicX[self.main_key] == 1:
            pass
        elif self.title_dicX[self.main_key] != 1:
            
            inp = self.main_key
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
            
            inp = self.main_key
            value = copy
            location = 0
            value_large_Then_inp = [x for x in self.title_dicX.keys() if self.title_dicX[x] > location]
            for key in value_large_Then_inp:
                self.title_dicX[key] = self.title_dicX[key]+1
            
            self.title_box.insert(location, value)
            
            self.title_dicX[inp] = location+1
            self.main_key_list = self.title_box[self.title_dicX[self.main_key]-1]
        
        old_main_key_list = self.title_box.pop(self.title_dicX[self.main_key]-1)
        self.main_key_list = [x for x in dicKeysortVal(self.title_dicX, '<-') if x != self.main_key]
        title_box = listXYreverse(self.title_box)
        title_box.insert(self.title_dicX[self.main_key]-1, self.main_key_list)
        self.title_box = title_box
        
        dicX = {x:y for x, y in [(x, old_main_key_list.index(x)+2) for x in old_main_key_list]}
        dicX[self.main_key] = 1
        self.title_dicX = dicX
        
        dicY = {x:y for x, y in [(x, self.main_key_list.index(x)+1) for x in self.main_key_list]}
        self.title_dicY = dicY
        
        title_tuple = [(x[1],x[0]) for x in self.title_dicX.items()]
        title_tuple.sort()
        self.title_tuple = title_tuple
        self.title_list = dicKeysortVal({x:y for x, y in self.title_dic.items() if x != self.main_key}, '<-')
        self.len_X = len(self.title_dicX)
        self.len_Y = len(self.title_dicY)
        
        return
    
    def mainkeyChange(self, new_main_key = ""):
        
        if new_main_key in self.title_dicX.keys():
            self.main_key = new_main_key
            old_main_key_list = self.main_key_list
            
            M_ = {}
            P = (self.title_dicX, self.title_dicY, self.title_box)
            E = listXYreverse(self.title_box)
            for key in old_main_key_list:
                result = E[P[1][key]-1]
                M_[result[self.title_dicX[new_main_key]-1]] = self.title_dicY[key]
            self.title_dicY = M_
            
            self.title_tuple = self.title_tuple
            self.title_dicX = self.title_dicX
            self.title_box = self.title_box
            self.main_key_list = self.title_box[self.title_dicX[new_main_key]-1]
            self.title_list = dicKeysortVal({x:y for x, y in self.title_dic.items() if x != self.main_key}, '<-')
            self.len_X = len(self.title_dicX)
            self.len_Y = len(self.title_dicY)
            
        else:
            print('new_main_key = ('+new_main_key+') is not in title_dicX, please double check~!')
            print('print location: ')
            input()
            raise KeyboardInterrupt
        
        return
    
#-----------------------------------------------------------------------------------------------

def mixYOrderfile(tabfilewithtitle_list = [], outputFile_Name = '.txt'):
    
    check_title_list = []
    for TabT in tabfilewithtitle_list:
        check_title_list.append(TabT.title_tuple.__str__())
    
    if len(list(set(check_title_list))) != 1:
        print('The title of tables are different, please double check~!')
        print('print location: ')
        input()
        raise KeyboardInterrupt
    else:
        pass
    
    table = tabfilewithtitle()
    table.build(TabT.title_dicX, TabT.main_key)
    
    new_title_box = []
    i = 1
    for title in table.title_dicX.keys():
        
        E_ = []
        for TabT in tabfilewithtitle_list:
            E_ = E_ + TabT.title_box[i-1]
            
        new_title_box.append(E_)
        i = i + 1
    table.title_box = new_title_box
    
    new_title_dicY_list = []
    new_main_key_list = []
    max_len_Y = 0
    for TabT in tabfilewithtitle_list:
        
        for item in TabT.title_dicY.items():
            
            key = item[0]
            value = item[1]+max_len_Y
            
            new_title_dicY_list.append((key, value))
            new_main_key_list.append(key)
        
        max_len_Y = max_len_Y + len(TabT.title_dicY)
    table.title_dicY = {x:y for x, y in new_title_dicY_list}
    
    table.main_key_list = new_main_key_list
    table.len_Y = len(table.title_dicY)
    
    if outputFile_Name != '.txt':
        
        f = open(outputFile_Name, 'w')
        f.write(table.report())
        f.close()
        
        h = hashlib.new('ripemd160')
        f = open(outputFile_Name, 'rb')
        h.update(f.read())
        f.close()
        
        target_path = os.getcwd()
        
        if os.path.exists(target_path+r'/__pycache__') == False:
            os.makedirs(target_path+r'/__pycache__')
        os.chdir(target_path+r'/__pycache__')
        
        try:
            Pickle_box = open('tabfilewithtitle_cache', 'rb')
            file_box = pickle.load(Pickle_box)
            Pickle_box.close()
        except:
            file_box = [{}, {}, {}]
        
        f = open(outputFile_Name+'.tabfilewithtitle.picklebox', 'wb')
        pickle.dump([table.title_tuple, table.title_dicX, table.title_dicY, table.title_box,
                     table.main_key, table.main_key_list, table.len_X, table.len_Y], f, )
        f.close()
        
        file_box[0][outputFile_Name] = h.hexdigest()
        file_box[1][outputFile_Name] = table.main_key
        file_box[2][outputFile_Name] = outputFile_Name+'.tabfilewithtitle.picklebox'
        
        Pickle_box = open('tabfilewithtitle_cache', 'wb')
        pickle.dump(file_box, Pickle_box, )
        Pickle_box.close()
        
        os.chdir(target_path)
        
    else:
        pass
    
    return table

#-----------------------------------------------------------------------------------------------

def newOrderfile(title_list = ['order', 'A', 'B' ], group_list = [['a', 'b'], ], outputFile_Name = '.txt'):
    
    table = tabfilewithtitle()
    
    i = 1
    title = {}
    for tit in title_list:
        title[tit] = i
        
        i = i + 1
    table.build('X', title, 'order')
    
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
        
        h = hashlib.new('ripemd160')
        f = open(outputFile_Name, 'rb')
        h.update(f.read())
        f.close()
        
        target_path = os.getcwd()
        
        if os.path.exists(target_path+r'/__pycache__') == False:
            os.makedirs(target_path+r'/__pycache__')
        os.chdir(target_path+r'/__pycache__')
        
        try:
            Pickle_box = open('tabfilewithtitle_cache', 'rb')
            file_box = pickle.load(Pickle_box)
            Pickle_box.close()
        except:
            file_box = [{}, {}, {}]
        
        f = open(outputFile_Name+'.tabfilewithtitle.picklebox', 'wb')
        pickle.dump([table.title_tuple, table.title_dicX, table.title_dicY, table.title_box,
                     table.main_key, table.main_key_list, table.len_X, table.len_Y], f, )
        f.close()
        
        file_box[0][outputFile_Name] = h.hexdigest()
        file_box[1][outputFile_Name] = table.main_key
        file_box[2][outputFile_Name] = outputFile_Name+'.tabfilewithtitle.picklebox'
        
        Pickle_box = open('tabfilewithtitle_cache', 'wb')
        pickle.dump(file_box, Pickle_box, )
        Pickle_box.close()
        
        os.chdir(target_path)
        
    else:
        pass
    
    return table

#-----------------------------------------------------------------------------------------------

def tabfileaddorder(file_name = '.txt', title_list = "['order', 'A', 'B'] or 'X'", outputFile_Name = '.txt'):
    
    if title_list == 'X':
        B = open(file_name,'r').read()
        Line_list = B.split('\n')[0:-1]
        
        E = []
        i = 1
        E.append('\t'.join(['order', Line_list[0]]))
        for Line in Line_list[1:]:
            E.append('\t'.join([str(i)+'_', Line]))
            i = i + 1
        
        result = '\n'.join(E)
        display = '\n'.join([result, str(len(E))])

    else:
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

def dicKeysortVal(dictionary = {}, pattern = '-> or <-'):
    result = []
    BtupleA = []
    if pattern == '->':
        
        for x in dictionary.keys():
            BtupleA.append((dictionary[x],x))
        
        BtupleA.sort()
        BtupleA.reverse()
        for x in BtupleA:
            result.append(x[1])
        
    elif pattern == '<-':
        
        for x in dictionary.keys():
            BtupleA.append((dictionary[x],x))
        
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
        result = result.replace('t', 'u')
    elif seq_type == 'DNA':
        result = result.replace('U', 'T')
        result = result.replace('u', 't')
    else:
        result = result
    
    return result

#-----------------------------------------------------------------------------------------------

def ntmatch(A_seq = '', B_seq = ''):
    
    if len(A_seq) >= len(B_seq):
        
        template = A_seq
        seq = B_seq
        marker1 = 'A'
        marker2 = 'B'
        
    elif len(B_seq) >= len(A_seq):
        
        template = B_seq
        seq = A_seq
        marker1 = 'B'
        marker2 = 'A'
    
    def seqreturn(marker, loac):
        if marker == 'A':
            seq = B_seq
            template = A_seq
            
            new_seq = '-'*loac + seq
            new_template = template + '-'*(len(new_seq)-len(template))
            
            if len(new_template) != len(new_seq):
                new_seq = new_seq + '-'*(len(template)-len(new_seq))
            
            new_A = new_template
            new_B = new_seq
            
        elif marker == 'B':
            seq = A_seq
            template = B_seq
            
            new_seq = '-'*loac + seq
            new_template = template + '-'*(len(new_seq)-len(template))
            
            if len(new_template) != len(new_seq):
                new_seq = new_seq + '-'*(len(template)-len(new_seq))
            
            new_A = new_seq
            new_B =  new_template
        return new_A, new_B
    
    def locationscore(template, seq, marker = 'A or B'):
        score_list = []
        location = 0
        while location < (len(template)/1):
            
            seq_location = 0
            score = 0
            
            new_seq = '-'*location + seq
            new_template = template + '-'*(len(new_seq)-len(template))
            
            i = 1
            for x in seq:
                
                if x == new_template[location+i-1]:
                    score = score + 1
                    
                i = i + 1
            
            if score > int(len(template)/1):
                R = seqreturn(marker, location)
                result = (R[0], R[1])
                return result
            
            score_list.append((score, location, marker))
            
            location = location + 1
        return score_list
    
    score_list_1 = locationscore(template, seq, marker1)
    if isinstance(score_list_1,tuple) == True:
        return score_list_1
    score_list_2 = locationscore(seq, template, marker2)
    if isinstance(score_list_2,tuple) == True:
        return score_list_2
    
    score_list = score_list_1 + score_list_2
    
    score_list.sort()
    marker = score_list[-1][2]
    loac = score_list[-1][1]
    
    R = seqreturn(marker, loac)
    
    result = (R[0], R[1])
    
    return result

#-----------------------------------------------------------------------------------------------

class LoopingTime:
    
    def __init__(self, ):
        self.looptitle = ""
        self.total_loop = 0
        self.mark = '|'
        self.passingtime_list = []
        self.needTime = ""
        self.switch = 'off'
    
    def open(self, looptitle = "", total_list = []):
        self.looptitle = looptitle
        self.total_loop = len(total_list)
    
    def showEndTime(self, passingtime = 'seconds from loopend cut loopin', switch = 'off'):
        
        self.passingtime_list.append(passingtime)
        if len(self.passingtime_list) < self.total_loop:
            ending = '\r'
        elif len(self.passingtime_list) == self.total_loop:
            ending = '\n'
        else:
            print('Looptime Error')
            raise KeyboardInterrupt
        
        mean = self.meanSD(self.passingtime_list)[0]
        self.needTime = mean*self.total_loop
        if switch == 'on':
            display = '{0}\t{1:25}\t{2}/100\tTotal time(m): {3:.8}\tStill time(m): {4:.8}'
            display = display.format(self.looptitle,
                                     self.mark*int(20/self.total_loop*len(self.passingtime_list)),
                                     str(int(100/self.total_loop*len(self.passingtime_list))),
                                     str(self.needTime/60),
                                     str(mean*(self.total_loop-len(self.passingtime_list))/60))
        else:
            display = '{0}\t{1:25}\t{2}/100'
            display = display.format(self.looptitle,
                                     self.mark*int(20/self.total_loop*len(self.passingtime_list)),
                                     str(int(100/self.total_loop*len(self.passingtime_list))))
        
        #print(display, end = ending)

    def clean(self,):
        self.looptitle = ""
        self.total_loop = 0
        self.mark = '|'
        self.passingtime_list = []
        self.needTime = ""

    def meanSD(self, E = [0.0]):
        
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
            mean = float(i/len(E))
            
            for x in E:
                j = j + (float(x)-mean)**2
            SD = float(j/len(E))**(0.5)
        
        return mean, SD
