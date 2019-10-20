import sys, os
import numpy as np 
import matplotlib.pyplot as plt
from collections import defaultdict
# from pathlib import Path



d = defaultdict(list)
directory = "../results/"
for root, dirs, files in sorted(os.walk(directory)):
	# print root[11:]
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
				# d[file[:3]].append(root[11:] +":"+ str(num))
				d[file[:3]].append(num)
		continue		

for key, value in sorted(d.iteritems()) :
	print key
	prefetchers = ['berti','bingo','bouquet','default','enhancing','multi-lop','next-line','pangloss','sangam','sangam++','t-skid']
	x = [1,2,3,4,5,6,7,8,9,10,11]
	fig = plt.figure()
	ax = plt.axes()
	plt.xticks(x, prefetchers,rotation=30, fontsize=8)
	plt.plot(x, value)
	plt.xlabel('Prefetchers from DPC3',fontsize=14)
	plt.ylabel('Cumulative IPC',fontsize=14)
	title = "Benchmark number: "+key
	fig.suptitle(title,fontsize=15)

	plt.savefig(key+'.pdf',bbox_inches = "tight")
