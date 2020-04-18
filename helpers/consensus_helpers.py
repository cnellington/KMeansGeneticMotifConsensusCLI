"""
Caleb Ellington
2020/02/27

simple consensus MEME algorithm and helpers
"""


import numpy as np
import math


dna = {'A': 0,
       'T': 1,
       'C': 2,
       'G': 3}

q = [0.35, 0.35, 0.15, 0.15]


def parse_fasta(filepath):
    with open(filepath, 'r') as file:
        rows = file.readlines()
        data = []
        i = 1
        seq = ""
        while i < len(rows):
            if '>' in rows[i]:
                data.append(seq)
                seq = ""
            else:
                seq += rows[i].strip('\n')
            i += 1
        data.append(seq)
        return data

class SimpleConsensus:

    def __init__(self, samples, k, T):
        self.samples = samples
        self.k = k
        self.T = T

    def run(self):
        best_motifs = []

        for t in range(1, len(self.samples)):
            motifs = []
            if t == 1:
                # First iteration, create T initial motifs from first 2 sequences
                seq1 = self.samples[t-1]
                seq2 = self.samples[t]
                for i in range(len(seq1)-self.k):
                    for j in range(len(seq2)-self.k):
                        x1 = seq1[i:i+self.k]
                        x2 = seq2[j:j+self.k]
                        motif = [x1, x2]
                        motifs.append(motif)
            else:
                # Append one k-mer per iteration
                seq = self.samples[t]
                for i in range(len(seq)-self.k):
                    x = seq[i:i+self.k]
                    for motif in best_motifs:
                        motifs.append(motif + [x])
            # Only retain T best motifs
            best_motifs = self._get_best_motifs(motifs, self.T)

        return best_motifs

    def _get_best_motifs(self, motifs, T):
        scores = np.array([self._get_information_content(motif) for motif in motifs])
        best = np.argsort(-scores)[:T]
        best_motifs = [motifs[i] for i in best]
        return best_motifs

    def _get_information_content(self, motif):
        pwm = self._get_pwm(motif)
        ic = 0
        for i in range(len(pwm)):
            for j in range(4):
                wbk = pwm[i][j]
                qb = q[j]
                if wbk == 0:
                    continue
                ic += wbk*math.log10(wbk/qb)
        return ic

    def _get_pwm(self, motif):
        n = len(motif)
        k = len(motif[0])
        pwm = np.zeros((k, 4))
        for i in range(n):
            for j in range(k):
                base = dna[motif[i][j]]
                pwm[j][base] += 1/n
        return pwm

