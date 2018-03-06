# 786

from __future__ import print_function
import sys
import pysam
from collections import defaultdict
import numpy as np
from bisect import bisect_left

def search(y, donor):
	m = bisect_left(donor[0], y)
	if m < len(donor) and donor[0][m] == y:
		return donor[1][m]
	else: 
		return donor[1][m - 1] + (y - donor[0][m - 1])

def parse(path):
	print("parsing {}".format(path), file=sys.stderr)
	reference = []
	deleted_ranges = set() # here no thing should map!

	with open(path) as f:
		for l in f:
			if l[0] == '#': continue
			l = l.strip().split()
			event = l[2]
			ref = l[3]
			alt = l[4]
			idx = int(l[1])
			extra = l[7].split(';')
			if 'small_indel' in event or 'sv_indel' in event:
				if len(ref) > len(alt): # deletion
					deleted_ranges |= set(range(idx + 1, idx + len(ref) - len(alt)))
					reference.append((idx + 1, 'D', len(ref) - len(alt) )) 
				else:
					reference.append((idx + 1, 'I', len(alt) - len(ref) )) # insert N bases
			elif 'sim_dup' in event:
				l = int(extra[1].split('=')[1])
				target = int(extra[3].split('=')[1].split(':')[1])
				reference.append((target, 'I', l))
			elif 'sim_inv' in event:
				y = int(extra[1].split('=')[1])
				reference.append((idx, 'R', y))
			elif 'sim_trans' in event:
				te = [idx]
				for x in range(5):
					l = next(f)
					l = l.strip().split()
					te += [int(l[1])]
				reference.append((te[0], 'T', te[1], te[2], te[3], te[4], te[5]))
	donor = {0: 0}
	reference = sorted(reference)
	donor_idx = 0
	for Q in reference:		
		if Q[1] == 'I':
			donor_idx += Q[2]
			donor[Q[0] + donor_idx] = Q[0]
		elif Q[1] == 'D':
			donor_idx -= Q[2]
			donor[Q[0] + donor_idx] = Q[0]
		elif Q[1] == 'R':
			#print(Q[0], Q[2], Q[2]-Q[0])
			for e in range(Q[0], Q[2]):
				donor[e + donor_idx] = Q[0] + Q[2] - e - 1
			donor[Q[2]+donor_idx] = Q[2]
		elif Q[1] == 'T':
			a = list(range(Q[2], Q[3])) + list(range(Q[4], Q[5]))
			b = list(range(Q[4], Q[5])) + list(range(Q[2], Q[3]))
			donor[Q[0] + donor_idx] = Q[0]
			for x, y in zip(b, a):
				donor[y + donor_idx] = x
			donor[Q[6] + donor_idx] = Q[6]

	donor = sorted(donor.items())
	donor = ([x[0] for x in donor], [x[1] for x in donor], deleted_ranges)
	return donor

donor1 = parse(sys.argv[1])
donor2 = parse(sys.argv[2])

print("parsing sam {}".format(sys.argv[3]), file=sys.stderr)

correct = defaultdict(lambda: [0,0])
# wrong_rc = 0
with pysam.AlignmentFile(sys.argv[3]) as sam:
	for read in sam.fetch():
		name = read.qname.split('_')
		pos = read.reference_start

		if read.flag & 0x100 or not read.cigar: # more than 2 mappings...
			continue
		if read.reference_id == -1:
			continue


 # read reverse strand (0x10)
 # first in pair (0x40)
		pair_pos = 0 if (read.flag & 0x40) != 0 else 1 # first in strand
		rev_comp = 1 if (read.flag & 0x10) != 0 else 0 # paired

		d = name[0]
		opos = list(map(int, name[1:3]))
		opos = opos[pair_pos]
		orev = list(map(int, name[3:5]))
		orev = orev[pair_pos]
		oopos = opos
		if d[-1] == '1': opos = search(opos, donor1)
		if d[-1] == '2': opos = search(opos, donor2)

		# if orev != rev_comp:
		# 	wrong_rc += 1
			# print(orev, rev_comp)
			# print(read)
			# print("strand")
			# exit(1)
		if abs(pos - opos) < 10:
			correct[read.mapq][0] += 1
		else:
			print(opos, orev, read.tostring(sam))
		correct[read.mapq][1] += 1
		#if sum(y[1] for y in correct.values()) > 100000: break
for mq, (c, t) in sorted(correct.items()):
	print("{} {} {} {}".format(mq, c, t, 100.0*c/float(t)), file=sys.stderr)
c = sum(y[0] for y in correct.values())
t = sum(y[1] for y in correct.values())
print("{} {} {} {}".format('total', c, t, 100.0*c/float(t)), file=sys.stderr)