import sys
ifp = open(sys.argv[1])
ofp = open(sys.argv[1]+'.pe.bed','w')
line1 = ifp.readline()
name1 = line1.split()[3].split('/')[0]
for line2 in ifp:
    name2 = line2.split()[3].split('/')[0]
    if name2 == name1:
        print >> ofp,line1.rstrip()
        print >> ofp,line2.rstrip()
    else:
        line1 = line2
        name1 = name2
ifp.close()


