import os, Raxpy3Libbasic_2013_07_30, RaxLib, time

oldmodule = Raxpy3Libbasic_2013_07_30
newmodule = RaxLib

t = time.time()
try:
    TabT1 = oldmodule.tabfilewithtitle()
    TabT1.open('TEST1_COAD_tumor.txt', 'order')
finally:
    t1 = time.time()-t
print(t1)

t = time.time()
try:
    TabT2 = newmodule.tabfilewithtitle()
    TabT2.open('TEST1_COAD_tumor.txt', 'order')
finally:
    t2 = time.time()-t
print(t2)

t = time.time()
try:
    f = open('TEST1_COAD_tumor.txt', 'r').read()
finally:
    t3 = time.time()-t
print(t3)


