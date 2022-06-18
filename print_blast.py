import xml.etree.ElementTree as ET
from numpy import identity

import pandas as pd


xml_tree = ET.parse('sequences/results/blast.out')
xml_root = xml_tree.getroot()
proteins = list(xml_root.find('BlastOutput_iterations').iter('Iteration'))
data = []
for protein in proteins:
    hits = list(protein.find('Iteration_hits').iter('Hit'))[:20]
    for hit in hits:
        protein_number = int(protein.find('Iteration_query-def').text.split('#')[1].split(' ')[0])
        orf = protein.find('Iteration_query-def').text.split('ORF ')[1].split(' ')[0]
        matched_protein_id = hit.find('Hit_def').text.split(' ')[0]
        matched_protein_name = hit.find('Hit_def').text.split('=')[1].split(';')[0].split('[')[0]
        matched_protein_species = hit.find('Hit_def').text.split('[')[1].split(']')[0]
        identity_percentage = int(hit.find('Hit_hsps').find('Hsp').find('Hsp_identity').text)/int(protein.find('Iteration_query-len').text)
        bitscore = float(hit.find('Hit_hsps').find('Hsp').find('Hsp_bit-score').text)
        evalue = float(hit.find('Hit_hsps').find('Hsp').find('Hsp_evalue').text)
        if bitscore > 40:
            data.append((
                protein_number, 
                orf,
                matched_protein_id, 
                matched_protein_name, 
                matched_protein_species, 
                identity_percentage,
                bitscore,
                evalue
            ))

df = pd.DataFrame(data, 
    columns=[
        'Protein',
        'ORF', 
        'Matched Protein ID', 
        'Matched Protein name', 
        'Matched Protein Species', 
        'Identity %',
        'Bitscore',
        'E-value'
    ])

df.to_csv('sequences/results/results.csv')

