from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def probable_orf(proteins: dict):
    '''Returns the most probable ORF. This supposes the longest protein is the most probable one'''
    longest_protein = 0
    orf = 0
    for k, v in proteins.items():
        if len(v) > longest_protein:
            longest_protein = len(v)
            orf = k
    return orf

def get_longest_protein_by_orf(record: SeqRecord) -> dict:
    '''Returns a dict with keys +1, +2, +3, -1, -2, -3. And values the longest protein to that ORF'''
    proteins = {}
    for strand, nucleotids in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(1, 4):
            longest_protein = 0

            # frame = 1 -> no neuclotids skipped
            # frame = 2 -> first 2 neuclotids skipped
            # frame = 3 -> first neuclotid skipped
            record_to_translate = nucleotids[(-frame + 4) % 3:]
            
            # Make the sequence length a multiple of 3 (biopython shows warning if not)
            record_to_translate = record_to_translate[:3*(len(record_to_translate) // 3)] 

            for protein in record_to_translate.translate().split('*'): # Translate the sequence and split by Stop codons (*)
                protein = protein[protein.find('M'):] # Look for start of protein if exists (M = Start)
                if protein.startswith('M') and len(protein) > longest_protein: # Register only the longest protein in the ORF
                    orf = strand * frame
                    orf = f'+{orf}' if orf > 0 else f'{orf}'
                    proteins[orf] = protein
                    longest_protein = len(protein)
    return proteins

def get_all_proteins_by_orf(record: SeqRecord) -> dict:
    '''Returns a dict with keys +1, +2, +3, -1, -2, -3. And values the all the proteins to that ORF'''
    proteins = {}
    for orf in ['+1', '+2', '+3', '-1', '-2', '-3']:
        proteins[orf] = []
    for strand, nucleotids in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(1, 4):
            longest_protein = 0

            # frame = 1 -> no neuclotids skipped
            # frame = 2 -> first 2 neuclotids skipped
            # frame = 3 -> first neuclotid skipped
            record_to_translate = nucleotids[(-frame + 4) % 3:]
            
            # Make the sequence length a multiple of 3 (biopython shows warning if not)
            record_to_translate = record_to_translate[:3*(len(record_to_translate) // 3)] 

            for protein in record_to_translate.translate().split('*'): # Translate the sequence and split by Stop codons (*)
                protein = protein[protein.find('M'):] # Look for start of protein if exists (M = Start)
                if protein.startswith('M'): # Register only the longest protein in the ORF
                    orf = strand * frame
                    orf = f'+{orf}' if orf > 0 else f'{orf}'
                    proteins[orf].append(protein)
    return proteins

def get_len_protein(prot_tuple):
    return len(prot_tuple[0])

# TODO: Be able to change gb_file path and output by command line
if __name__ == '__main__':
    gb_file = 'sequences/unknown_virus.gb'
    output_path = 'sequences/results/protein.fasta'

    for gb_record in SeqIO.parse(open(gb_file,'r'), 'genbank') :
        proteins = get_all_proteins_by_orf(gb_record)
        protein_qty = 0
        all_proteins = []
        for orf,proteins_in_orf in proteins.items():
            for protein in proteins_in_orf:
                if len(protein)*3/len(gb_record.seq) > 0.01:
                    all_proteins.append((protein, orf))

    chosen_proteins = sorted(all_proteins, key=get_len_protein, reverse=True)[:10]
    records = []
    for i, (protein, orf) in enumerate(chosen_proteins):
        print(f"{orf}: {len(protein)*3/len(gb_record.seq) * 100:.2f}% of gen used for protein")
        selected_protein = protein
        record = SeqRecord(protein, description=f'Protein #{i+1} translated from {gb_record.id} from ORF {orf}.')
        records.append(record)
    SeqIO.write(records, output_path, 'fasta')


