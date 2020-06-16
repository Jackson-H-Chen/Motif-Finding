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
parser.add_argument('--tolerance', required=False, type=float, default=0.9,
	metavar='<float>', help='threshold [%(default).3f]')
arg = parser.parse_args()

SIM = {
	'A': {'A':0.97, 'C':0.01, 'G':0.01, 'T':0.01},
	'C': {'A':0.01, 'C':0.97, 'G':0.01, 'T':0.01},
	'G': {'A':0.01, 'C':0.01, 'G':0.97, 'T':0.01},
	'T': {'A':0.01, 'C':0.01, 'G':0.01, 'T':0.97},
	'R': {'A':0.49, 'C':0.01, 'G':0.49, 'T':0.01},
	'Y': {'A':0.01, 'C':0.49, 'G':0.01, 'T':0.49},
	'M': {'A':0.49, 'C':0.49, 'G':0.01, 'T':0.01},
	'K': {'A':0.01, 'C':0.01, 'G':0.49, 'T':0.49},
	'W': {'A':0.49, 'C':0.01, 'G':0.01, 'T':0.49},
	'S': {'A':0.01, 'C':0.49, 'G':0.49, 'T':0.01},
	'B': {'A':0.01, 'C':0.33, 'G':0.33, 'T':0.33},
	'D': {'A':0.33, 'C':0.01, 'G':0.33, 'T':0.33},
	'H': {'A':0.33, 'C':0.33, 'G':0.01, 'T':0.33},
	'V': {'A':0.33, 'C':0.33, 'G':0.33, 'T':0.01},
	'N': {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25},
}
	
REG = {
	'A': 'A',
	'C': 'C',
	'G': 'G',
	'T': 'T',
	'R': '[AG]',
	'Y': '[CT]',
	'M': '[AC]',
	'K': '[GT]',
	'W': '[AT]',
	'S': '[GC]',
	'B': '[CGT]',
	'D': '[AGT]',
	'H': '[ACT]',
	'V': '[ACG]',
	'N': '.',
}

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

def search_pwm(dna, pwm, t):
	found = []
	lp = len(pwm)
	for i in range(len(dna) - lp + 1):
		seq = dna[i:i+lp]
		p = 1
		for j in range(len(seq)):
			p *= pwm[j][seq[j]]
		if p > t: found.append(i)
	return found

def search_re(dna, rst):
	found = []
	for match in re.finditer(rst, dna):
		found.append(match.start())
	return found

def motif_entropy(pwm):
	sum_H = 0
	for i in range(len(pwm)):
		h = 0
		for nt in pwm[i]:
			if pwm[i][nt] == 1.0:
				h += 2
			elif pwm[i][nt] != 0:
				h -= pwm[i][nt] * math.log2(pwm[i][nt])
		sum_H += h
	return sum_H

def dkl(p, q):
	d = 0
	for i in p:
		if p[i] != 0:
			d += p[i] * math.log2(p[i]/q[i])
	return d

def make_regex(pwm):
	ntstr = ''
	restr = ''
	for i in range(len(pwm)):
		min_d = 1e6
		min_nt = None
		for nt in SIM:
			d = dkl(pwm[i], SIM[nt])
			if d < min_d:
				min_d = d
				min_nt = nt
		ntstr += min_nt
		restr += REG[min_nt]
	
	return restr, ntstr


if __name__ == '__main__':
	motifs = read_motifs(arg.motifs)
	for name, seq in read_fasta(arg.fasta):
		match = re.search('^(\S+)', name)
		id = match[1]
		for mid in motifs:
			pwm = motifs[mid]
			h = motif_entropy(pwm)
			rst, nst = make_regex(pwm)
			freq = 1 / 2 ** h
			tol = arg.tolerance ** len(pwm)
			t = freq * tol
			pwm_sites = search_pwm(seq, pwm, t)
			re_sites = search_re(seq, rst)
			print(mid, nst, 'regex', re_sites, 'pwm', pwm_sites)



