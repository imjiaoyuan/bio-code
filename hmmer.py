import subprocess
import input

def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True) 
    output, error = process.communicate() 

    if process.returncode != 0:
        print(f'Error occurred: {error}')
    else:
        print(f'Success!\n{output}')

# 创建你的命令
hmmbuild_command = 'hmmbuild build.hmm PFseed.txt'
hmmsearch_command = 'hmmsearch build.hmm ./protein_sequence/{}.fasta > result.out'.format(input.species)

# 运行你的命令
run_command(hmmbuild_command)
run_command(hmmsearch_command)