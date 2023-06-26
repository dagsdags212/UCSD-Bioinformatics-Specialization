from itertools import product
from random import shuffle
from deBruijn import *

def kmer_composition(text: str, k: int) -> list[str]:
    """Generate the k-mer composition of a string text."""
    composition: list = []
    for i in range(len(text)-k+1):
        composition.append(text[i:i+k])
    composition.sort()
    return composition

def reconstruct_genome_from_path(genome_path: list[str]) -> str:
    """Reconstruct a genome from a collection of substrings/k-mers
    such that the last k-1 symbols ofthe current k-mer matches the
    first k-1 symbols of the next k-mer."""
    genome: str = genome_path[0]
    for i in range(1, len(genome_path)):
        if genome_path[i-1][1:] == genome_path[i][:-1]:
            genome += genome_path[i][-1]
        else:
            raise ValueError('Neighboring kmers must have a matching prefix and suffix')
    return genome

def find_k_universal_string(k: int) -> str:
    ls = ([''.join(l) for l in product('01', repeat=k)])
    dbg = DeBruijinGraph(ls, k)
    path = list(dbg.eulerian_path())

    kstring = path[0]
    for i in range(1, len(path)-(k-1)):
        kstring += path[i][-1]
    return kstring

