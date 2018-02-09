from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

#ext_modules = [Extension('RaxLib', ['Raxpy3Libbasic_2013_07_30.pyx'])]
#ext_modules = [Extension('RaxTCGALib', ['Raxpy3LibdataTCGA_2013_07_26.pyx'])]
#ext_modules = [Extension('RaxTCGALib', ['Raxpy3LibdataTCGA_2013_08_18.pyx'])]
#ext_modules = [Extension('RaxLib', ['Raxpy3Libbasic_2013_09_14.pyx'])]
#ext_modules = [Extension('RaxLib', ['Raxpy3Libbasic_2013_10_02.pyx'])]
#ext_modules = [Extension('RaxLib', ['Raxpy3Libbasic_2013_11_12.pyx'])]
#ext_modules = [Extension('RaxLib', ['Raxpy3Libbasic_2013_11_27.pyx'])]
#ext_modules = [Extension('RaxLib', ['Raxpy3Libbasic_2013_12_02.pyx'])]
#ext_modules = [Extension('RaxLib', ['Raxpy3Libbasic_2013_12_31.pyx'])]
#ext_modules = [Extension('RaxTCGALib', ['Raxpy3LibdataTCGA_2013_08_18.pyx'])]
#ext_modules = [Extension('Raxsqlit3Lib', ['Raxpy3LibSqlite3_2014_01_05.pyx'])]
#ext_modules = [Extension('RaxLib', ['Raxpy3Libbasic_2014_01_09.pyx'])]
ext_modules = [Extension('RaxTCGALib', ['Raxpy3LibdataTCGA_2013_12_12.pyx'])]



setup(name = 'my cython libarary test',
      cmdclass = {'build_ext':build_ext},
      ext_modules = ext_modules)

#it can't run on Windows = =" shit!
#'python3 setup.py build_ext --inplace' on Linux is work