def data_download(PF_number):
    import wget
    hmm_url = "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{}?annotation=hmm".format(PF_number)
    wget.download(hmm_url, 'hmm.gz')

    PFseed_url = "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{}/?annotation=alignment:seed&download".format(PF_number)
    wget.download(PFseed_url, 'PFseed.gz')

def decompress_data():
    import gzip

    with gzip.open('hmm.gz', 'rb') as f:
        content = f.read().decode('utf-8')
    with open('./model.hmm', 'w') as f:
        f.write(content)

    with gzip.open('PFseed.gz', 'rb') as f:
        content = f.read().decode('utf-8')
    with open('./PFseed.txt', 'w') as f:
        f.write(content)