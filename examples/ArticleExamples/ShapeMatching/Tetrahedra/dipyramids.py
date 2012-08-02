#!/usr/bin/python

import glob
import subprocess
import string

for i in range(47,65): 
	sum=0
	total=0
	arg = './../../../example_data/tetrahedra/P' + str(i) + '_*.txt'
	files = glob.glob(arg)
	sum1 = 0.0
	sum2 = 0.0
	sum3 = 0.0
	
	for j in files:
		cmd = "./dipyramids"
		pproc = subprocess.Popen([cmd, "512", j, "0.7", "0.75"], stdout=subprocess.PIPE)
		ans = pproc.stdout.read()
		ans = string.split(ans)
		sum1 = sum1 + float(ans[0])
		sum2 = sum2 + float(ans[1])
		sum3 = sum3 + float(ans[2])
		total = total+1
	
	print str(i) + ' ' + str(sum1/total) + ' ' + str(sum2/total) + ' ' +  str(sum3/total)

