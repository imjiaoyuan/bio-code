from modules.preparation import data_download
from modules.preparation import decompress_data
from modules.hmmer import run_command
from modules.extract_id import extract_id
from modules.extract_seq import extract_sequences
import config

hmmer_result = './result.out'
id_list = './id_list.txt'
target_species_proteins = "./protein_seq/{}.fasta".format(config.species)
target_gene_proteins = 'protein_for_target_id.fasta'
hmmbuild_command = 'hmmbuild model.hmm PFseed.txt'
hmmsearch_command = 'hmmsearch model.hmm {} > {}'.format(target_species_proteins, hmmer_result)

data_download(config.PF_number)
decompress_data()
run_command(hmmbuild_command)
run_command(hmmsearch_command)
extract_id(hmmer_result, id_list, config.evaluation_threshold)
extract_sequences(id_list, target_species_proteins, target_gene_proteins)