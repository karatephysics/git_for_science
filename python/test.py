#!/usr/bin/env python

import os
import os.path
import sys
import datetime
import subprocess
import shutil
import time

from logging import getLogger,StreamHandler,Formatter,DEBUG
logger = getLogger(__name__)
handler = StreamHandler(sys.stderr)
handler.setLevel(DEBUG)
formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(process)d - %(message)s')
handler.setFormatter(formatter)
logger.setLevel(DEBUG)
logger.addHandler(handler)

FLUTE_EXE = ' ../orion_flute '
FLUTE_ADDRESS = ' -S: -a:239.1.1.2 -p:40000 -i:192.168.10.10 '
FLUTE_SEND_PARAM = ' -n:10 -l:info  -E:1000 -N:1000 -r:300000 -B:data'

fdt = 1
toi = 1

f = open('/Users/Tetsu/Desktop/share/python/fdt.txt','r')  

for row in f:
	print row
f.close

f = open('/Users/Tetsu/Desktop/share/python/fdt.txt','w')  
f.write(str(fdt+3))
f.close()
toi = fdt
logger.debug('fdt =%d' % (fdt))

"""
print ' -f:%d -F:"%d;;;;;;;data.json;;;"' % (fdt, toi)
"""

cmdline = FLUTE_EXE + FLUTE_ADDRESS + FLUTE_SEND_PARAM + ' -f:%d -F:"%d;;;;;;;data.json;;;"' % (fdt, toi)
logger.debug(cmdline)
retcode = subprocess.call(cmdline, shell=True)
toi += 1
fdt += 1

cmdline = FLUTE_EXE + FLUTE_ADDRESS + FLUTE_SEND_PARAM + ' -f:%d -F:"%d;;;;;;;data.json;;;"' % (fdt, toi)
logger.debug(cmdline)
retcode = subprocess.call(cmdline, shell=True)
toi += 1
fdt += 1
