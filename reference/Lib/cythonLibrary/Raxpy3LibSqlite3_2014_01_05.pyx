import Raxpy3Libbasic_2013_12_31 as Lib
import sqlite3 as DB
import sys as sys

ori_path = Lib.os.getcwd()
'''
history = 2014.01.05
'''

#-----------------------------------------------------------------------------------------------
#tabfilewithtitle and sqlite3 only use 'str' and 'float' types
#-----------------------------------------------------------------------------------------------

def textTOsqlite(file_name = '.txt', sql_name = '.sqldb', table_name = ' '):
    
    TabT = Lib.tabfilewithtitle()
    TabT.open(file_name, 'order')
    
    title_tuple = [x[1] for x in TabT.title_tuple if x[1] != 'order']
    group_list = [TabT.readrow('Y', x)[1:] for x in TabT.main_key_list]
    
    SQLDB = DB.connect(sql_name)
    cr_obj = SQLDB.cursor()
    cr_obj.execute('''create table '''+table_name+''' '''+str(tuple(title_tuple))+''';''')
    
    for group in group_list:
        cr_obj.execute('''insert into '''+table_name+''' values ('''+','.join('?'*len(title_tuple))+''');''', group)
        
    SQLDB.commit()
    SQLDB.close()
    
    return 
    
def sqliteTOtext(sql_name = '.sqldb', table_name = ' ', file_name = '.txt', limit = 0):
    
    SQLDB = DB.connect(sql_name)
    cr_obj = SQLDB.cursor()
    
    title = cr_obj.execute('''pragma table_info('''+table_name+''')''').fetchall()
    title_list = [x[1] for x in title]
    
    E_ = []
    for row in cr_obj.execute('''select * from '''+table_name+''';'''):
        E_.append(list(row))
    
    Lib.newOrderfile(['order']+title_list, E_, file_name)
    
    return

def textAddOrder(inp_file_name = '.txt', title_list = "['order', 'A', 'B'] or 'X'", oup_file_name = '.txt'):
    Lib.tabfileaddorder(inp_file_name, title_list, oup_file_name)
    

def dirTOsqlite(source_dir = '/../..'):
    
    sql_name = source_dir.split('/')[-1]+'.sqldb'
    file_list = [x for x in Lib.os.listdir(source_dir) if x.split('.')[-1] == 'txt']
    
    print(file_list)
    
    Lib.os.chdir(source_dir)
    E_ = []
    for file_name in file_list:
        first_line = open(file_name, 'r').readline()
        if first_line.split('\t')[0] != 'order':
            E_.append(file_name)
    
    for file_name in file_list:
        if file_name in E_:
            print('There is no title in '+file_name+'. And processing excluded !')
        else:
            textTOsqlite(file_name, sql_name, file_name.split('.')[0])
            print(file_name+' finished !')
    
    return

def sqliteTOdir(source_dir_sql_name = '/../..sqldb'):
    
    #if '.sqldb' in source_dir_sql_name and '/' in source_dir_sql_name:
    #    
    #    Lib.os.chdir('/'.join(source_dir_sql_name.split('/')[:-1]))
    #    
    #    sql_name = source_dir_sql_name.split('/')[-1]
    #    SQLDB = DB.connect(sql_name)
    #    cr_obj = SQLDB.cursor()
    #    
    #    TT = cr_obj.execute('''.table''').fetchall()
    #    
    #    
    #    
    #    print(TT)
    #    
    #    
    #else:
    #    print('Please enter the directory with .sqldb name !')
    #    pass
    
    return
    
if __name__ == '__main__':
    
    if sys.argv[1] == 'textAddOrder':
        textAddOrder(sys.argv[2] , 'X', sys.argv[2].split('.')[0]+'2.txt')
    else:
        print('not run')
    
    quit()