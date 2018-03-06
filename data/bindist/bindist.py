# 786

import pysam, os, re, copy, subprocess, tempfile, glob, sys
import collections 
import numpy as np
import pandas as pd
import cPickle as pickle

path = sys.argv[1]
out = sys.argv[2]
print >>sys.stderr, 'Reading', path, 'to', out

# needs sorted by name
clouds = {} #collections.defaultdict(lambda: collections.defaultdict(int))
# barcode -> [pos -> count]
# prev_name = ''
pckl = out + '.pickle'
if os.path.exists(pckl):
    clouds = pickle.load(open(pckl))
else:
    with pysam.AlignmentFile(path) as sam:
        prev_chr = ''
        for read in sam.fetch():
            if read.flag & 0x100: # "alternate": more than 2 mappings...
                continue
            if not read.cigartuples: # not mapped
                continue
            if not read.has_tag('BX'):
                continue
            if not read.has_tag('XC'):
                continue
            if read.get_tag('XF') != 0:
                continue
            if read.reference_id == -1:
                continue
            # if read.qname == prev_name:
            #     continue    
            if read.next_reference_id != read.reference_id:
                continue
            if read.next_reference_start < read.reference_start:
                continue
            if read.reference_name != prev_chr:
                print >>sys.stderr, out, '@', read.reference_name, '...'
                prev_chr = read.reference_name
            # if read.reference_name != 'chr22':
                # continue
            
            barcode = read.get_tag('BX')
            cloud = read.get_tag('XC')
            cloud = '{}/{}'.format(barcode, cloud)
            if cloud not in clouds:
                clouds[cloud] = collections.defaultdict(int)
            dx = clouds[cloud]
            
            dx[read.reference_start] += 1
            if abs(read.next_reference_start - read.reference_start) < 1000: # Filter out crap
                dx[read.next_reference_start] += 1
            # prev_name = read.qname  
    pickle.dump(clouds, open(pckl, "wb" ))

print >>sys.stderr, out, 'Reading done'
print >>sys.stderr, out, 'Cloud count', len(clouds)

# clouds = pickle.load(open("_x_read"))


BIN_SIZE = 1000
bins = [[] for i in range(1000)]
# print >>sys.stderr, 'Using bin size', BIN_SIZE
for cn, c in clouds.items():
    pos = sorted(c.keys())
    reads = {p - pos[0]: c[p] for p in pos}
    bn = [0 for _ in range((pos[-1] - pos[0]) // BIN_SIZE + 1)]
    for p, n in reads.items():
        bn[p // BIN_SIZE] += n
    for i, x in enumerate(bn):
        bins[i].append(x) # = np.array(x)

dists = {}
for i in range(len(bins)):
    c = collections.Counter(bins[i])
    total = sum(w for _, w in c.items())
    if total == 0: continue
    dists[i] = np.array([100.0*c[x]/float(total) for x in range(20)])
    for x in c.keys():
        if x >= 20: dists[i][19] += 100.0*c[x]/float(total)
df = pd.DataFrame.from_dict(dists).reset_index()
print df.head()


df.to_csv(out)

print >>sys.stderr, out, 'All done'

