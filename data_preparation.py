import requests
import input

hmm = "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{}?annotation=hmm".format(input.PF_number)

hmm_content = requests.get(hmm)

hmm_content = hmm_content.text

hmm_file = open("build.hmm", "w")

hmm_file.write(hmm_content)

hmm_file.close()

PFseed = "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{}/?annotation=alignment:seed&download".format(input.PF_number)

PFseed_content = requests.get(PFseed)

PFseed_content = PFseed_content.text

PFseed_file = open("PFseed.txt", "w")

PFseed_file.write(PFseed_content)

PFseed_file.close()