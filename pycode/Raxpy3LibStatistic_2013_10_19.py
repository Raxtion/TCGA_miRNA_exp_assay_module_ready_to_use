'''
history = 2013.08.02, 2013.08.04, 2013.10.19
'''

#-----------------------------------------------------------------------------------------------
import math, RaxLib
#-----------------------------------------------------------------------------------------------
RaxLib = RaxLib

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
    
    aa = (float(A_5p) + float(A_3p))
    bb = (float(A_5p) + float(A_3p))
    cc = (float(B_5p) + float(B_3p))
    dd = (float(B_5p) + float(B_3p))
    if aa == 0.0 or bb == 0.0 or cc == 0.0 or dd == 0.0:
        critical_value = 'Denominator = 0'
        return critical_value
    
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
        mean = float(i/len(E))
        
        for x in E:
            j = j + (float(x)-mean)**2
        SD = float(j/len(E))**(0.5)
    
    return mean, SD

#-----------------------------------------------------------------------------------------------
#error function from http://stackoverflow.com/questions/457408/is-there-an-easily-available-implementation-of-erf-for-python
#Error at x less than 1.5*10^-7
def erf(x): 
    # save the sign of x
    sign = 1 if x >= 0 else -1
    x = abs(x)

    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
    return sign*y # erf(-x) = -erf(x)

#-----------------------------------------------------------------------------------------------

class normaldist:
    
    def __init__(self):
        self.mean = 0.0
        self.N = 0.0
        self.SD = 0.0
        self.top = 0.0
        self.low = 0.0
    
    def open(self, E = [0.0]):
        self.N = len(E)
        self.mean, self.SD = meanSD(E)
    
    def pdf(self, x = 0.0):
        y = (1/(self.SD*2*math.pi**(0.5)))*math.e**((-1*0.5)*((x-self.mean)/self.SD)**2)
        return y
    
    '''
    def antipdf(self, x = 0.0):
        y = -1*self.SD*(2*math.log(x*self.SD*(2*math.e)**(0.5)))**(0.5)+self.mean
        return y
    '''
    
    def cdf(self, x = 0.0):
        y = (0.5)*(1+(erf((x-self.mean)/(self.SD*(2**(0.5))))))
        return y
    
    def probAtoB(self, A = 0.0, B = 0.0):
        y0 = self.cdf(A)
        y1 = self.cdf(B)
        probability = y1 - y0
        return probability
    
    def itemScoreTrans(self, inp = [], top = 0.0, low = 0.0):
        oup = [(self.cdf(x)*(top-low)/100)+low for x in inp]
        return oup
    
#-----------------------------------------------------------------------------------------------

