import requests
import input

hmm = "http://pfam-legacy.xfam.org/family/{}/hmm".format(input.PF_number)  # 文件的URL地址

hmm_content = requests.get(hmm)  # Send a GET request to obtain file content

hmm_content = hmm_content.text  # Obtain the text content of the file

hmm_file = open("build.hmm", "w")

hmm_file.write(hmm_content)  # Write the obtained content to a file

hmm_file.close()

PFseed = "http://pfam-legacy.xfam.org/family/{}/alignment/seed/format?format=stockholm&alnType=seed&order=t&case=l&gaps=default&download=1".format(input.PF_number)

PFseed_content = requests.get(PFseed)  # Send a GET request to obtain file content

PFseed_content = PFseed_content.text  # Obtain the text content of the file

PFseed_file = open("PFseed.txt", "w")

PFseed_file.write(PFseed_content)  # Write the obtained content to a file

PFseed_file.close()