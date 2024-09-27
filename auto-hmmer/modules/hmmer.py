def run_command(command):
    import subprocess

    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True) 
    output, error = process.communicate() 

    if process.returncode != 0:
        print(f'Error occurred: {error}')
    else:
        print(f'Success!\n{output}')
