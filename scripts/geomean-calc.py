from __future__ import division
import sys, os
import numpy as np 
# import matplotlib.pyplot as plt
from collections import defaultdict
# from pathlib import Path



d = defaultdict(list)
directory = "../results/"
for root, dirs, files in sorted(os.walk(directory)):
	# print root[11:]
	avg = 0.0
	gavg = 1.0
	for file in sorted(files):
		filename = root +"/"+ file
		for line in open(filename):
			if "CPU 0 cumulative IPC:" in line:
				ipc = []
				flag = 0
				for x in line[22:]:
					if x == ' ':
						flag = 1
					if flag == 0:
						ipc.append(x)
				num = float(''.join(map(str,ipc)))
				avg = avg + num
				gavg = gavg*num
				# d[file[:3]].append(root[11:] +":"+ str(num))
				ind = file.find('B.')
				# d[file[:ind+1]].append(num)
		continue
	print root[11:],avg/20.0,gavg**(1/20.0)		
