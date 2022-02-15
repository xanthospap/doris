#! /usr/bin/python

## 
##  Enclose a C header file, within an 'extern C' macro (within ifdef guards)
##  Apply to C-header files to make them compatible within C++
##  All output is directed to stdout
## 

import sys, os

if not os.path.isfile(sys.argv[1]):
    print('ERROR. Failed to locate file {:}'.format(sys.argv[1]))

head = """#ifdef __cplusplus
extern "C" {
#endif
"""

tail = """    #ifdef __cplusplus
}
#endif
"""

print(head)
with open(sys.argv[1], 'r') as fin:
    for line in fin.readlines():
        print(line.strip())
print(tail)
