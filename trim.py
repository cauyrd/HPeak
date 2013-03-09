#|/usr/bin/env python


import numpy as np
import sys
import gzip as gz

if '.gz' in sys.argv[1]:
	ifp = gz.open(sys.argv[1])
else:
	ifp = open(sys.argv[1])
fraglen=[]
chr = ['chr'+str(i) for i in range(1,23)] + ['chrX']+ ['chrY']
while True:
    line1 = ifp.readline().rstrip()
    if not line1:
        break
    line2 = ifp.readline().rstrip()
    chrname1=line1.split()[0]
    chrname2=line2.split()[0]
    if chrname1 in chr and chrname1 == chrname2:
        pos = [int(line1.split()[1]),int(line1.split()[2]),int(line2.split()[1]),int(line2.split()[2])]
        dis = max(pos)-min(pos)
        fraglen.append(dis)
ifp.close()
#ofp1 = open(sys.argv[1]+'.summary.txt','w')

fraglen=np.array(sorted(fraglen))
chop=float(sys.argv[2])     #### choose a chop percentage, will cut both tails by the chosen percentage
i1 = len(fraglen) * chop/100 + 0.5 -1
i2 = len(fraglen) * (100-chop)/100 + 0.5 -1
k1 = np.floor(i1)
f1 = i1- np.floor(i1)
bottom = (1-f1) * fraglen[k1] + f1 * fraglen[k1+1]
k2 = np.floor(i2)
f2 = i2 - np.floor(i2)
top = (1-f2) * fraglen[k2] + f2 * fraglen[k2+1]   ### 95th percentile of frag length
#mu1 = np.mean(fraglen)                           #### mean before chopping
#sigma1 = np.mean(fraglen)		      ### std before chopping
fraglen=fraglen[(fraglen>bottom)*(fraglen<top)]  ### trim off top 5 and bottom 5 

#mu = np.mean(fraglen)
#sigma = np.std(fraglen)

#print >> ofp1, 'Chop by',chop, 'percent top and bottom'
#print >> ofp1, 'Before choppoing', 'mean is', mu1, 'and std is', sigma1
#print >> ofp1, 'After choppoing' 'mean is', mu, 'and std is', sigma
#ofp1.close()


ofp = open(sys.argv[1]+'.trimmed.bed','w')
ifp = open(sys.argv[1], 'r')

while True:
    line1 = ifp.readline().rstrip()
    if not line1:
        break
    line2 = ifp.readline().rstrip()
    chrname1=line1.split()[0]
    chrname2=line2.split()[0]
    if chrname1 in chr and chrname1 == chrname2:
        pos = [int(line1.split()[1]),int(line1.split()[2]),int(line2.split()[1]),int(line2.split()[2])]
        dis = max(pos)-min(pos)
        if dis >= bottom and dis <= top:
            print >> ofp,line1.rstrip()
            print >> ofp,line2.rstrip()


