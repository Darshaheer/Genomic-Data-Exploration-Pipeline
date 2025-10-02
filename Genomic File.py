# Installing Libraries
from Bio import Entrez, SeqIO
from collections import defaultdict
from Bio.PDB import PDBParser
import requests
import matplotlib.pyplot as plt
import seaborn as sns 

# Downloading FASTA FILE
Entrez.email = "shaheerdar085@gmail.com"
handle = Entrez.einfo()
record = Entrez.read(handle)
handle.close()
print(record['DbList'])

handles = Entrez.esearch(db="nucleotide", term="PFCRT[Gene] AND Plasmodium falciparum[Organism]", retmax="40")
rec_list = Entrez.read(handles)
handles.close()
print(rec_list['Count'])
print(len(rec_list['IdList']))
print(rec_list['IdList'])

id_list = rec_list['IdList']
handlee = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb")
records = list(SeqIO.parse(handlee, "gb"))
handlee.close()

for rec in records:
    if rec.id == "KM288867":
        break
    print(rec.name)
    print(rec.description)
    print(str(rec.seq))

# FASTA FILE EXPLORATION
recu = list(SeqIO.parse("Long.txt", "fasta"))
print(recu[0].id, recu[0].description)
print(recu[0].seq)

print("Complement:", recu[0].seq.complement())
print("Reverse Complement:", recu[0].seq.reverse_complement())

mrna = recu[0].seq.transcribe()
print("Transcribed:", mrna)

protein = recu[0].seq.translate()
print("Translated:", protein)

def gc_contents(sequence):
    return (sequence.count('G') + sequence.count('C'))/ len(sequence) * 100

countt = 0
for r in recu:
    print("ID:", r.id, "GC%:", gc_contents(r.seq))
    countt += 1
    if countt == 5:
        break

def nucleotide_count(sequences):
    return {nucleotide: sequences.count('nucleotide') for nucleotide in 'ATGC'}
    for i, s in enumerate(recu[:5]):
        print("The Nucleotide Count of Sequence", i, "is:", nucleotide_count(s.seq))

def longest_sequence(file_path, file_format):
    sequence = list(SeqIO.parse(file_path, file_format))
    longest_seq = max(sequence, key=lambda recu:len(recu.seq))
    return longest_seq.id, longest_seq.seq
print("The longest sequence in the file is::", longest_sequence("Long.txt", "fasta"))

def sequence_distribution_plot(file_path, file_format):
    lenghts = [len(rec.seq) for rec in SeqIO.parse(file_path, file_format)]
    plt.hist(lenghts, bins=20, edgecolor='black')
    plt.title("Sequence Distribution Plot")
    plt.xlabel("Sequence Lenght")
    plt.ylabel("Frequency")
    plt.show()

sequence_distribution_plot("Long.txt", "fasta")

output = SeqIO.SeqRecord(recu[0].seq, id='12345', description='My Extraction')
SeqIO.write(output, "output.txt", "fasta")
print("The file has been written")

# FASTQ FILE EXPLORATION
docs = SeqIO.parse("SRR390728_1.fastq", "fastq")
doc = next(docs)

print(doc)
print(doc.id, doc.description, doc.seq, doc.letter_annotations)

count = defaultdict(int)
n_cnt = defaultdict(int)

for dox in docs:
    for letter in dox.seq:
        count[letter] += 1
    for i, letters in enumerate(dox.seq):
        pos = i + 1
        if letters == "N":
            n_cnt[pos] += 1

tot = sum(count.values())
seq_len = max(n_cnt.keys())
positions = range(1, seq_len + 1)

for letter, cnt in count.items():
    print("%s: %.2f %d" % (letter, 100. * cnt/tot, cnt))

fig, ax = plt.subplots(figsize=(16,9))
ax.plot(positions, [n_cnt[x] for x in positions])
ax.set_title("Frequency Distribution of N Base")
ax.set_xlabel("Positions on the Reads")
ax.set_ylabel("Counts of N Base")
ax.set_xlim(1, seq_len)
plt.show()

report = SeqIO.parse("SRR390728_1.fastq", "fastq")
cnt_qual = defaultdict(int)
qual_pos = defaultdict(list)


for rep in reports:
    for i, qual in enumerate(rep.letter_annotations['phred_quality']):
        if i < 25:
            cnt_qual[qual] += 1
    for i, quals in enumerate(rep.letter_annotations['phred_quality']):
        pos = i + 1
        if i < 25 or quals == 40:
            continue
        qual_pos[pos].append(quals)

tots = sum(cnt_qual.values())
vps = []
poses = list(qual_pos.keys())
poses.sort()
for poss in poses:
    vps.append(qual_pos[pos])

for qual, count in cnt_qual.items():
    print("%d: %.2f %d" % (qual, 100. * count/tots, count))

fig, ax = plt.subplots(figsize=(16,9))
sns.bloxplot(data=vps, ax=ax)
plt.set_title("Frequency of Quality Scores Across Positions")
plt.set_xlabel("Position on Reads")
plt.set_ylabel("Phred Scores")
plt.xticklabels([str(x) for x in range(26, max(qual_pos.keys()) + 1)])
plt.show()


        