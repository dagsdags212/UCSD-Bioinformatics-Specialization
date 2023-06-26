import sys
from central_dogma import reverse_complement, translate, find_coding_region

def parse_genome(path):
    genome = ''
    with open(path, 'r') as fh:
        for line in fh:
            genome += line.rstrip()
    fh.close()
    return genome

if __name__ == '__main__':
    peptide = 'VKLFPWFNQY'
    genome = parse_genome(sys.argv[1])
    result = find_coding_region(genome, peptide)
    print(result)
