import numpy as np
import json
def softmax(s, t=1.0):
    x = np.array(s)
    res = np.exp(x / t) / np.sum(np.exp(x / t), axis=0)
    return res.tolist()

def generate_GENE(population, gene_length=50):
    ge = np.random.choice(4, gene_length).tolist()
    GE = []
    for i in xrange(population):
        GE.append(np.random.choice(4, gene_length).tolist())
    return GE

def load_genes(dir_GE1, dir_GE2):
    def convert_str2int_list(s):
        return [int(i) for i in s]

    with open(dir_GE1, 'r') as f:
        GE1 = json.load(f)
    with open(dir_GE2, 'r') as f:
        GE2 = json.load(f)
    GE1, GE2 = GE1.keys(), GE2.keys()
    GE1 = [convert_str2int_list(i) for i in GE1]
    GE2 = [convert_str2int_list(i) for i in GE2]
    return GE1, GE2