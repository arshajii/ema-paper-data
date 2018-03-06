import os
import sys
import pysam
import subprocess
import shutil
import getopt
import shlex
from itertools import izip
from collections import defaultdict, Counter
import gzip

def fastqIter(fq):
  linesPerRead = 4
  def grouped(iterator):
    while True:
      vals = tuple(next(iterator, None) for _ in xrange(linesPerRead))
      if None not in vals:
        yield vals
      else:
        raise StopIteration

  assert os.path.isfile(fq)
  with open(fq) as f:
    for (qname_ln, bases, c, bqs) in grouped(f):
      # strip aux fields
      qname_w = qname_ln.split()[0]
      # strip /1,/2
      qname = qname_w.split('/')[0]
      # strip @
      qname = qname[1:]
      txt = ''.join([qname_ln, bases, c, bqs])
      yield (qname, txt)

  raise StopIteration

def demultiplex(inputFq_path, outputDir_path):

  def getFname(qname):
    fname = qname.split(':')[0]
    if not (fname.endswith('fq') or fname.endswith('fastq')):
      return None
    else:
      return fname
  def getHandle(fname):
    print 'creating new split group:', fname
    f = open(os.path.join(outputDir_path, fname), 'w')
    return f

  cfname = None
  f = None
  seen_set = set()
  for (qname, txt) in fastqIter(inputFq_path):
    # obtain target fname from qname
    nfname = getFname(qname)
    if nfname == None:
      print 'warning: read {0} not assigned target fastq group'.format(nfname)
      continue
    # open new handle if moved onto a new fastq group
    if cfname != nfname:
      if f != None:
        f.close()
      f = getHandle(nfname)
      # all reads for the same fastq group should be together
      assert nfname not in seen_set
      seen_set.add(nfname)
    cfname = nfname

    f.write(txt)

  # close tail file
  if f != None:
    f.close()

def main(argv):

  # argument parsing
  if len(argv) < 3:
    print 'fastqsplit.py <input fastq path> <output directory>'
    sys.exit(1)
  inputFq_path = argv[1]
  outputDir_path = argv[2]

  if not os.path.isfile(inputFq_path):
    print 'error: specified input must be valid file'
    sys.exit(2)
  if not os.path.isdir(outputDir_path):
    print 'error: specified output must be path to existing directory'
    sys.exit(3)

  demultiplex(inputFq_path, outputDir_path)


if __name__ == '__main__':
  main(sys.argv)

