#!/usr/bin/env python3

import argparse
import sys
import re
import math

parser = argparse.ArgumentParser(
	description='Motif finder.')
parser.add_argument('--motifs', required=True, type=str,
	metavar='<path>', help='file of motifs in MEME format')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='file of sequences in FASTA format')
parser.add_argument('--fp_rate', required=False, type=float, default=1.0,
	metavar='<float>', help='threshold [%(default).3f]')
arg = parser.parse_args()

def read_motifs(file):
	motifs = {}
	ntps = '(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)'
	with open(file) as fp:
		motif_id = None
		for line in fp.readlines():
			match = re.search('MOTIF (\S+)', line)
			if match != None:
				motif_id = match[1]
			match = re.search(ntps, line)
			if match != None:
				if motif_id not in motifs: motifs[motif_id] = []
				p = {'A':float(match[1]), 'C':float(match[2]), 'G':float(match[3]), 'T':float(match[4])}
				motifs[motif_id].append(p)
	return motifs

def read_fasta(filename):
	name = None
	seqs = []
	
	fp = None
	if filename == '-':
		fp = sys.stdin
	elif filename.endswith('.gz'):
		fp = gzip.open(filename, 'rt')
	else:
		fp = open(filename)

	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

def search_seq(dna, pwm):
	t = 1e-4
	found = []
	lp = len(pwm)
	for i in range(len(dna) - lp + 1):
		seq = dna[i:i+lp]
		p = 1
		for j in range(len(seq)):
			p *= pwm[j][seq[j]]
		if p > t: found.append(i)
	print(found)
	sys.exit()
	return found

def site_freq(pwm, rate):
	sum_H = 0
	for i in range(len(pwm)):
		h = 0
		for nt in pwm[i]:
			if pwm[i][nt] == 1.0:
				h += 2
			elif pwm[i][nt] != 0:
				h -= pwm[i][nt] * math.log2(pwm[i][nt])
		sum_H += h
	return 1 / 2 ** sum_H

mdb = read_motifs(arg.motifs)
for name, seq in read_fasta(arg.fasta):
	match = re.search('^(\S+)', name)
	id = match[1]
	for mid in mdb:
		f = site_freq(mdb[mid], arg.fp_rate)
		t = f * 1000
		print(mid, f, t)
	#	motifs = search_seq(seq, mdb[mid])
	#sys.exit()



