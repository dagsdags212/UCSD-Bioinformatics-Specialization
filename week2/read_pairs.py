from deBruijn import Node, DeBruijinGraph
import random
import sys

"""
Exercise Break: Generate the (3,2)-mer composition of TAATGCCATGGGATGTT in lexicographic order.
Include repeats, and return your answer as a list on a single line.
As a hint to help you with formatting, your answer should begin "(AAT|CAT) (ATG|ATG)..."
"""
def generate_read_pair_composition(text: str, k: int, d: int) -> list[str]:
    """Return a list of (k,d)-mers from a given text."""
    composition = []
    for i in range(len(text)-2*k-d+1):
        k1 = text[i:i+k]
        k2 = text[i+k+d:i+2*k+d]
        composition.append(k1+k2)
    composition.sort()
    return [f"({kdmer[:k]}|{kdmer[k:]})" for kdmer in composition]

def generate_string_from_gapped_genome(kdmers: list[str], sep: str = '-') -> str:
    k = len(kdmers[0].split(sep)[0])
    prefix, suffix = '', ''
    for i, kdmer in enumerate(kdmers):
        a, b = kdmer.split(sep)
        if i+1 == len(kdmers):
            prefix += a
            suffix += b
        else:
            prefix += a[:-1]
            suffix += b[:-1]
    return prefix + suffix[-(k+1):]

def string_spelled_by_gapped_patterns(kdmers: list[str], k: int, d: int, sep: str = '-'):
    kmers1, kmers2 = [], []
    for kdmer in kdmers:
        k1, k2 = kdmer.split(sep)
        kmers1.append(k1)
        kmers2.append(k2)
    prefix = ''.join([k1[0] for k1 in kmers1[:-1]]) + kmers1[-1]
    suffix = ''.join([k2[0] for k2 in kmers2[:-1]]) + kmers2[-1]

    for i in range(k+d+1, len(prefix)):
        if prefix[i] != suffix[i-k-d]:
            print("There is no string spelled by the gapped patterns")
            return False

    return prefix + suffix[-(k+d):]

"""
Solve the String Reconstruction from Read-Pairs Problem.

    1. Split the read-pairs into its k-mers and instantiate into nodes
    2. Concatenate nodes and create deBruijn graph
    3. Generate Eulerian path from graph
"""

def split_reads(read_pair: str, sep: str = '|') -> list[str]:
    k1, k2 = read_pair.split(sep)
    return k1, k2

def generate_dbg_from_paired_reads(paired_reads: list[str]) -> dict[tuple[str, str], list[tuple]]:
    reads = list(map(split_reads, paired_reads))
    graph = {(read[0], read[1]): [] for read in reads}

    for key in graph.keys():
        s1, s2 = key[0][1:], key[1][1:]
        for read in reads:
            if key == read:
                continue
            p1, p2 = read[0][:-1], read[1][:-1]
            if s1 == p1 and s2 == p2:
                graph[key].append(read)

    return graph

def eulerian_path_from_paired_reads(paired_reads: list[str]):
    g = generate_dbg_from_paired_reads(paired_reads)
    # find a random starting node and a fixed ending node
    for k, v in g.items():
        if v == []:
            end = k
    start = random.choice(list(g.keys()))
    # initialize stack and path
    stack = [start]
    path = [end]
    graph = g.copy()
    # walk path
    while stack:
        current = stack[-1]
        if len(graph[current]) > 0:
            child = graph[current].pop()
            if child != end:
                stack.append(child)
        else:
            n = stack.pop()
            path = [n] + path
    return path

def parse_file(path):
    with open(path, 'r') as fh:
        k, d = map(int, fh.readline().strip().split())
        paired_reads = fh.readline().strip().split()
    fh.close()
    return k, d, paired_reads

def main() -> None:
    k, d, paired_reads = parse_file(sys.argv[1])
    path = eulerian_path_from_paired_reads(paired_reads)

    while len(path) != len(paired_reads):
        path = eulerian_path_from_paired_reads(paired_reads)

    reads = [f"{read[0]}-{read[1]}" for read in path]
    seq = string_spelled_by_gapped_patterns(reads, k, d)
    print(seq)

if __name__ == '__main__':
    main()
