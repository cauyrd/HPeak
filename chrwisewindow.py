#!/usr/bin/python
'''
File: chrwisewindow.py
Author: Rendong Yang
Description: the file is a modified version of chrwisewindw.pl for handling paired-end data
Email: cauyrd@gmail.com
Date: Mon Dec  5 14:31:37 EST 2011
'''
import sys

def get_position(pair1, pair2):
    position = [pair1[1],pair2[1],pair1[2],pair2[2]]
    position = map(int,position)
    position.sort()
    start = position[0]
    length = position[-1]-position[0]
    return (start,length)

def get_update_bin(bins, fragments, resolution):
    '''this function is used to calculate the bins coverage for each chromosome'''
    for i in range(len(fragments)):
        # get the bin index of start of each fragment
        (start_div,start_mod) = divmod(fragments[i][0],resolution)
        (fraglen_div,fraglen_mod) = divmod(fragments[i][1],resolution)
        if start_mod == 0:
            start_index = start_div - 1
            start_len = 1
        else:
            start_index = start_div 
            start_len = (resolution - start_mod)+1
            if fraglen_div == 0:
                if start_len > fraglen_mod:
                    bins[start_index] += fraglen_mod/float(resolution)
                    continue
                else:
                    bins[start_index] += start_len/float(resolution)
                    bins[start_index+1] += (fraglen_mod-start_len)/float(resolution)
                    continue
        bins[start_index] += start_len/float(resolution)
        full_cover_bin = (fragments[i][1] - start_len)/resolution
        for j in range(start_index+1,start_index+1+full_cover_bin):
            bins[j] += 1
        try:
            bins[j+1] += (fragments[i][1]-start_len)%resolution/float(resolution)
        except NameError:
            bins[start_index+1] += (fragments[i][1]-start_len)%resolution/float(resolution)
        except:
            print 'Error in fill_cover_bin computing!'
            sys.exit(0)
import re
import os

sp = sys.argv[1]
if sp.lower() == 'human' or 'hg19':
    nchr = 24
    fold1 = 'hg19'
elif sp.lower() == 'hg18':
	nchr = 24
	fold1 = 'hg18'
elif sp.lower() == 'mouse' or 'mm9':
    nchr = 21
    fold1 = 'mm9'
else:
    print 'species error!'
    exit(0)
infilename = sys.argv[2]
INFILE_pre = open(infilename)
infilename2 = INFILE_pre.readline().rstrip()
INFILE_pre.close()

INFILE = open(infilename2)
format = sys.argv[3]
readlenth = int(sys.argv[4])
resolution = int(sys.argv[5])
frontname = sys.argv[6]+'.windowhitscount.'
selectname = sys.argv[6]+'.select.'

format = sys.argv[3]
# get the all chromosome size of selected reference genome
tmp = len(sys.argv[0].split('/'))
path = '/'.join(sys.argv[0].split('/')[0:(tmp-1)])
CHROMOFILE = open(path+'data/'+fold1+'/chromoends.txt')
chromosize = [0]*nchr
for i,line in enumerate(CHROMOFILE):
    if sp.lower() == 'bed':
        chromosize[i] = int(line.rstrip()) - 1
    else:
        chromosize[i] = int(line.rstrip())

# get the fragments coverage from paired-end data
fragments = {}
chr = ['chr'+str(i) for i in range(1,23)] + ['chrX']+ ['chrY']
while True:
    pair1 = INFILE.readline().split()
    if not pair1:
        break
    pair2 = INFILE.readline().split()
    if pair1[0] != pair2[0] or pair1[0] not in chr:
        continue
    p = re.compile('chr(.*)')
    m = p.match(pair1[0])
    if m.group(1) == 'X':
        chr_index = nchr-2
    elif m.group(1) == 'Y':
        chr_index = nchr-1
    else:
        chr_index = int(m.group(1))-1
    if chr_index not in fragments:
        fragments[chr_index] = [get_position(pair1,pair2)]
    else:
        fragments[chr_index].append(get_position(pair1,pair2))

# output the window coverage for each chromosome in windowhitcount and select file
threshold =1 
for i in range(nchr):
    upperlimit = chromosize[i]/resolution
    if chromosize[i] % resolution == 0:
        chr_bins = [0]*upperlimit
    else:
        chr_bins = [0]*(upperlimit+1)
    #chr_bins = [0]*upperlimit if chromosize[i]%resolution==0 else [0]*(upperlimit+1)
    OUTFILE = open(frontname+'chr'+str(i+1)+'.txt','w')
    OUTSELECT = open(selectname+'chr'+str(i+1)+'.txt','w')
    if i in fragments.keys():
        get_update_bin(chr_bins,fragments[i],resolution)
    for j in range(upperlimit):
        start = resolution*j + 1
        end = start + resolution -1
        if (i == (nchr-2)):
            print >> OUTFILE, 'chrX\t'+str(start)+'\t'+str(end)+'\t'+str(chr_bins[j])
        elif (i == (nchr - 1)):
            print >> OUTFILE, 'chrY\t'+str(start)+'\t'+str(end)+'\t'+str(chr_bins[j])
        else:
            print >> OUTFILE, 'chr'+str(i+1)+'\t'+str(start)+'\t'+str(end)+'\t'+str(chr_bins[j])
        if chr_bins[j] > threshold:
            print >> OUTSELECT, str(j)+'\t'+str(chr_bins[j])

