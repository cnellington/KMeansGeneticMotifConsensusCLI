"""
Caleb Ellington
2020/02/25

simple consensus motif finder
"""

import argparse
from os.path import isfile
from helpers.consensus_helpers import *


"""
Usage: kmeans.py <expression_data.txt> <initial_clusters.txt> [--verbose]
"""
def main():
    parser = argparse.ArgumentParser(description="simple consensus motif finding")

    parser.add_argument(
        "--fasta",
        "-fa",
        action="store",
        default="data/simple.fa",
        help="FASTA format sequences to scan for motifs",
    )

    parser.add_argument(
        "--motif-length",
        "-k",
        action="store",
        type=int,
        default=4,
        help="k-mer length of motif candidates",
    )

    parser.add_argument(
        "--trials",
        "-T",
        action="store",
        type=int,
        default=1,
        help="Number of motifs to find",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=False,
        help="Verbose run",
    )

    args = parser.parse_args()
    if not isfile(args.fasta):
        print("Invalid filepath")
        exit(1)

    samples = parse_fasta(args.fasta)
    meme = SimpleConsensus(samples, args.motif_length, args.trials)
    motifs = meme.run()
    for motif in motifs:
        print(meme._get_information_content(motif))
    print(f"Best Motif: {motifs[0]}")
    print(f"Best Score: {meme._get_information_content(motifs[0])}")

    # Q 2a
    # motif = ['ACG', 'ACT', 'ATC']
    # pwm = meme._get_pwm(motif)
    # ic = meme._get_information_content(motif)
    # print("PWM: ")
    # for base in dna.keys():
    #     line = base+"\t"
    #     for i in range(len(pwm)):
    #         line += str(round(pwm[i][dna[base]], 2))
    #         line += "\t"
    #     print(line)
    # print(f"IC: {ic}")

if __name__ == "__main__":
    main()
