#!/usr/bin/env python

import numpy as np
from addgns4loc_config import *

# Generate gaussion noise
# (Mean, STD, numbers)
noisex = np.random.normal(nsmean,nslcstd/100.0,nevn)
noisey = np.random.normal(nsmean,nslcstd/110.0,nevn)
noisez = np.random.normal(nsmean,nslcstd,nevn)
noiset = np.random.normal(nsmean,nsotstd,nevn)

print np.mean(noisex),np.mean(noisey),np.mean(noisez),np.mean(noiset)
print np.std(noisex),np.std(noisey),np.std(noisez),np.std(noiset)


# Add noise to p data
inpf = open(inpfp,'r')
outf = open(outfp,'w')
for i in range(1,nstap+1):
	hd1 = inpf.readline()
	hd2 = inpf.readline()
	outf.write(hd1)
	outf.write(hd2)

	hd1tmp = hd1.split()
	# print hd1tmp[0]
	
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
		ot = 0.0
		
		elon = elon + noisex[evnid-1]
		if(elon < minlon) or (elon > maxlon):
			elon = elon - 2*noisex[evnid-1]

		elan = elan + noisey[evnid-1]
		if(elan < minlan) or (elan > maxlan):
			elan = elan - 2*noisey[evnid-1]

		edep = edep + noisez[evnid-1]
		if(edep < minz) or (edep > maxz):
			edep = edep - 2*noisez[evnid-1]
			
		ot = ot + noiset[evnid-1]
		
		format1 = '%15.11f %9.5f %8.5f %6.3f %4d %9.5f %8.5f %7.4f %9.6f\n'
		outf.write(format1 % (tt,slon,slan,sdep,evnid,elon,elan,edep,ot))

inpf.close()
outf.close()


# Add noise to s data
inpf = open(inpfs,'r')
outf = open(outfs,'w')
for i in range(1,nstas+1):
	hd1 = inpf.readline()
	hd2 = inpf.readline()
	outf.write(hd1)
	outf.write(hd2)

	hd1tmp = hd1.split()
	# print hd1tmp[0]
	
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
		ot = 0.0
		
		elon = elon + noisex[evnid-1]
		if(elon < minlon) or (elon > maxlon):
			elon = elon - 2*noisex[evnid-1]

		elan = elan + noisey[evnid-1]
		if(elan < minlan) or (elan > maxlan):
			elan = elan - 2*noisey[evnid-1]

		edep = edep + noisez[evnid-1]
		if(edep < minz) or (edep > maxz):
			edep = edep - 2*noisez[evnid-1]

		ot = ot + noiset[evnid-1]
		
		format2 = '%15.11f %9.5f %8.5f %6.3f %4d %9.5f %8.5f %7.4f %9.6f\n'
		outf.write(format1 % (tt,slon,slan,sdep,evnid,elon,elan,edep,ot))

inpf.close()
outf.close()
