#!/usr/bin/env python

import numpy as np

# Generate gaussion noise
# (Mean, STD, numbers
noise = np.random.normal(0.042,0.178,49568)

print np.mean(noise)
print np.std(noise)


# Add noise to data
inpf = open('syndataloc_1km10_1dinitcheck.txt','r')
outf = open('gns.txt','w')
k=-1
for i in range(1,36):
	hd1 = inpf.readline()
	hd2 = inpf.readline()
	outf.write(hd1)
	outf.write(hd2)

	hd1tmp = hd1.split()
	print hd1tmp[0]
	
	for j in range(1,int(hd1tmp[0])+1):
		data = inpf.readline()
		dtmp = data.split()

		tt = float(dtmp[0])
		slon = float(dtmp[1])
		slan = float(dtmp[2])
		sdep = float(dtmp[3])
		evnid = int(dtmp[4])
		elon = float(dtmp[5])
		elan = float(dtmp[6])
		edep = float(dtmp[7])
		ot = float(dtmp[8])
		
		k = k + 1
		tt = tt + noise[k]
		
		format1 = '%15.11f %9.5f %8.5f %6.3f %4d %9.5f %8.5f %7.4f %9.6f\n'
		outf.write(format1 % (tt,slon,slan,sdep,evnid,elon,elan,edep,ot))

inpf.close()
outf.close()
