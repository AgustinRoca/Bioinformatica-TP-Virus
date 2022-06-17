from Bio.Blast import NCBIWWW
from Bio import SeqIO

fasta = open('sequences/results/protein.fasta','r')
result_handle = NCBIWWW.qblast('blastp', 'swissprot', fasta)
with open(f'sequences/results/blast_remote.out', 'w') as save_file:
    blast_results = result_handle.read() 
    save_file.write(blast_results)