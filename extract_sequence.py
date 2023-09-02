from Bio import SeqIO
from input import species

# Define functions to extract sequences
def extract_sequences(id_list, fasta_file, out_file):
    with open(id_list, 'r') as id_handle, open(out_file, 'w') as out_handle:
        ids = [line.strip() for line in id_handle.readlines()]
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if name in ids:
                SeqIO.write(fasta, out_handle, 'fasta')

extract_sequences('id_list.txt', "./protein_sequence/{}.fasta".format(species), 'protein_for_target_id.fasta')