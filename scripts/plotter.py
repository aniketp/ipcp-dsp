from __future__ import division
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
				ind = file.find('B.')
				d[file[:ind+1]].append(num)
		continue		
cnt = 0
for key, value in sorted(d.iteritems()) :
	cnt = cnt + 1
	print key
	prefetchers = ['default','next-line','berti','bingo','bouquet','enhancing','multi-lop','pangloss','sangam','sangam++','t-skid']
	x = [1,2,3,4,5,6,7,8,9,10]
	min_value = min(value[1:])
	max_value = max(value[1:])
	# avg_value = sum(value)/len(value)
	# min_value = min_value - 0.001
	# min_value = round(min_value, 2) - 0.001
	# max_value = max_value + 0.01
	print min_value, max_value
	print value
	fig = plt.figure()
	ax = plt.axes()
	# axes = plt.gca()
	# axes.set_xlim([xmin,xmax])
	ax.set_ylim([min_value,max_value])
	plt.xticks(x, prefetchers[1:],rotation=30, fontsize=8)
	# plt.yticks(np.arange(min_value, max_value, 0.001))
	plt.bar(x, value[1:], width=0.3)
	# plt.xlabel('Prefetchers from DPC3',fontsize=14)
	if cnt%4 == 1:
		plt.ylabel('Cumulative IPC',fontsize=14)
	# title = "Benchmark: "+key
	# fig.suptitle(title,fontsize=15)
	# plt.show()
	plt.savefig(key+'.pdf',bbox_inches = "tight")
	# sys.exit(0)