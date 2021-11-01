import os
import sys

from Bio import SeqIO


# Author = Bill winnett
# email = wwinnett@iastate.edu


def get_prob(raw, first):
    prob = {}

    for seq in raw:
        for p in range(1, 5):
            if seq[p] == first:
                for i in range(p, p + 20 if p + 20 >= len(seq) - 1 else len(seq) - 1):
                    letter = seq[i]

                    if i not in prob:
                        prob[i] = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}

                    prob[i][letter] += 1
                break
            break
        continue

    for dic in prob:

        dic_total = sum(prob[dic].values())

        for elem in prob[dic]:
            prob[dic][elem] = prob[dic][elem] / int(dic_total)

    return prob


def derive_lcs(prob):
    r = 2

    tech_seq = ""

    for dic in prob:
        if prob[dic][max(prob[dic], key=lambda x: prob[dic][x])] > .6:
            tech_seq += max(prob[dic], key=lambda x: prob[dic][x])
        else:
            break
    print("Technical Seq: ", tech_seq)
    return tech_seq


def get_n_ratios(raw, target):
    rawnum = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}

    for seq in raw:
        if target in seq:
            str(seq).replace(target, '')

        for i in range(0, len(seq)):
            result = seq[i]

            rawnum[result] += 1

    total = sum(rawnum.values())
    return (rawnum['A'] / total, rawnum['C'] / total, rawnum['G'] / total, rawnum['T'] / total)


def get_num_target(raw, target):
    count = 0
    for read in raw:
        if target in read:
            count += 1

    return count


test_file = sys.argv[1]

fastq_reads = []
for record in SeqIO.parse(test_file, "fastq"):
    fastq_reads += [record.seq]

prob = get_prob(fastq_reads, "G")
tech_seq = derive_lcs(prob)

out_loc = input("Where do you want the file to go (path):")
if not os.path.isfile(out_loc):
    current_file = open(out_loc, "x")
current_file = open(out_loc, "w")

for record in SeqIO.parse(test_file, "fastq"):
    # if tech_seq in record.seq:
    indexes = [0, 0]

    for i in range(0, 5):

        if record.seq[i] == tech_seq[0] and record.seq[i + 1] == tech_seq[1] and record.seq[i + 2] == tech_seq[2]:
            for x in range(0, len(tech_seq)):
                if record.seq[i + x] != tech_seq[x] or x == len(tech_seq) - 1:
                    indexes = [i, i + x + 1]

    new_seq = record.seq[indexes[0]:indexes[1]]

    if len(new_seq) == 0:
        continue

    current_file.write(str("@" + record.id + " length=" + str(indexes[1] - indexes[0]) + "\n" + new_seq + "\n"))

    qual_score_new = ""
    for j in range(indexes[0], indexes[1]):
        qual_score_new += chr(record.letter_annotations["phred_quality"][j] + 33)

    # print(qual_score_new)
    current_file.write(str("+" + record.id + " length=" + str(indexes[1] - indexes[0]) + "\n" + qual_score_new + "\n"))

numkmers = 0
for seq in fastq_reads:
    numkmers += len(seq) - len(tech_seq) + 1

prob_mult = 1
ratios = get_n_ratios(fastq_reads, tech_seq)
print(ratios)
for i in range(len(tech_seq)):

    if tech_seq[i] == "A":
        prob_mult *= ratios[0]

    elif tech_seq[i] == "C":
        prob_mult *= ratios[1]

    elif tech_seq[i] == "G":
        prob_mult *= ratios[2]

    elif tech_seq[i] == "T":
        prob_mult *= ratios[3]

target_seq_observed = get_num_target(fastq_reads, tech_seq)

print("Target sequence is found in ", target_seq_observed, "out of ", len(fastq_reads),
      "sequences. Or ", target_seq_observed / len(fastq_reads), "%")
print("Expected number of target sequence is: " + str(prob_mult * numkmers))
print("Actual number of target sequence is: ", target_seq_observed)

print("Observed over/under abundance: ", target_seq_observed / (prob_mult * numkmers))
