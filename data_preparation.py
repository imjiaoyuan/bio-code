import requests
import input

hmm = "http://pfam-legacy.xfam.org/family/{}/hmm".format(input.PF_number)  # 文件的URL地址

hmm_content = requests.get(hmm)  # 发送GET请求获取文件内容

hmm_content = hmm_content.text  # 获取文件的文本内容

hmm_file = open("build.hmm", "w")

hmm_file.write(hmm_content)

hmm_file.close()

PFseed = "http://pfam-legacy.xfam.org/family/{}/alignment/seed/format?format=stockholm&alnType=seed&order=t&case=l&gaps=default&download=1".format(input.PF_number)

PFseed_content = requests.get(PFseed)  # 发送GET请求获取文件内容

PFseed_content = PFseed_content.text  # 获取文件的文本内容

PFseed_file = open("PFseed.txt", "w")

PFseed_file.write(PFseed_content)

PFseed_file.close()