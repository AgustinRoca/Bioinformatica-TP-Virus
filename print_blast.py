import xml.etree.ElementTree as ET


xml_tree = ET.parse('sequences/results/blast.out')
xml_root = xml_tree.getroot()
proteins = list(xml_root.find('BlastOutput_iterations').iter('Iteration'))
for protein in proteins[:1]:
    print('=====================================================================')
    print(protein.find('Iteration_query-def').text)
    print()
    hits = list(protein.find('Iteration_hits').iter('Hit'))[:10]
    for hit in hits:
        print(f"\t{hit.find('Hit_def').text}")
        print(f"\t{int(hit.find('Hit_hsps').find('Hsp').find('Hsp_identity').text)/int(protein.find('Iteration_query-len').text) * 100:.2f}% of exact match")
        print()
