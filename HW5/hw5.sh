#This script just complies the f2py code and runs hw5.py in ipython
#Jonathan Varkovitzky
#March 2, 2012

rm jameson.so

f2py -c -m jameson jameson.f90
ipython hw5.py