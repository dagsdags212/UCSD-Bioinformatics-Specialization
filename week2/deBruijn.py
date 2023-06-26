class Node:
    def __init__(self, kmer: str) -> None:
        self.kmer = kmer
        self.incoming: int = 0
        self.outgoing: int = 0

    def __hash__(self) -> int:
        return hash(self.kmer)

    def __str__(self) -> str:
        return self.kmer

    def __repr__(self) -> str:
        return f"[{self.incoming},{self.outgoing}]{self.kmer}"

    def is_balanced(self) -> bool:
        """Check if count for incoming and outgoing edges are equal."""
        return self.incoming == self.outgoing

    def is_semi_balanced(self) -> bool:
        """Check if counts for incoming and outgoing edges is with 1."""
        return abs(self.incoming - self.outgoing) == 1

class DeBruijinGraph:
    @staticmethod
    def chop(read, k):
        """Return composition of a string, with each element having length k."""
        return [read[i:i+k] for i in range(len(read)-(k-1))]

    def __init__(self, reads: list[str], k: int) -> None:
        self.G = {}
        self.nodes = {}
        self.k = k

        # build deBruijn graph
        for read in reads:
            for kmer in self.chop(read, k):
                prefix, suffix = kmer[:-1], kmer[1:]
                nodeL, nodeR = None, None
                if prefix in self.nodes:
                    nodeL = self.nodes[prefix]
                else:
                    nodeL = self.nodes[prefix] = Node(prefix)
                if suffix in self.nodes:
                    nodeR = self.nodes[suffix]
                else:
                    nodeR = self.nodes[suffix] = Node(suffix)
                nodeL.outgoing += 1
                nodeR.incoming += 1
                self.G.setdefault(nodeL, [])
                self.G[nodeL].append(nodeR)

        # check for balanced nodes
        self.num_semi = 0
        self.num_balanced = 0
        self.num_unbalanced = 0
        for node in self.nodes.values():
            if node.is_balanced():
                self.num_balanced += 1
            elif node.is_semi_balanced():
                if node.incoming == node.outgoing + 1:
                    self.tail = node
                if node.incoming == node.outgoing - 1:
                    self.head = node
                self.num_semi += 1
            else:
                self.num_unbalanced += 1

    def node_count(self) -> int:
        """Return number of nodes."""
        return len(self.nodes)

    def edge_count(self) -> int:
        """Return number of edges."""
        return len(self.G)

    def has_eulerian_path(self) -> bool:
        """Return True if all nodes are balanced excpeced for starting and ending nodes."""
        return self.num_unbalanced == 0 and self.num_semi == 2

    def has_eulerian_cycle(self) -> bool:
        """Return True if all nodes are balanced."""
        return self.num_unbalanced == 0 and self.num_semi == 0

    def is_eulerian(self) -> bool:
        """Return True if graph has Eulerian path or cycle."""
        return self.has_eulerian_path() or self.has_eulerian_cycle()

    def eulerian_path(self) -> list[str]:
        """Find and return Eulerian path or cycle."""
        assert self.is_eulerian()
        g = self.G.copy()
        path = []

        def __visit(node: Node) -> None:
            while len(g[node]) > 0:
                dest = g[node].pop()
                __visit(dest)
            path.append(node)

        if self.has_eulerian_path():
            assert self.head is not None
            assert self.tail is not None
            g.setdefault(self.tail, [])

        __visit(self.head) if self.has_eulerian_path() else __visit(next(iter(g.keys())))
        path = path[::-1]

        if self.has_eulerian_path():
            start = path.index(self.head)
            path = path[start:] + path[:start]

        return list(map(str, path))

def parse_file(path) -> tuple[int, list[str]]:
    with open(path, 'r') as fh:
        k = int(fh.readline().strip())
        kmers = fh.readline().strip().split()
    fh.close()
    return k, kmers

def path_to_genome(path) -> str:
    genome = path[0]
    for i in range(1, len(path)):
        genome += path[i][-1]
    return genome

def main() -> None:
    file_num = 8
    path = f'/home/dagsdags/home/courses/ucsd-bioinformatics-specialization/course-2/datasets/string_reconstruction/input_{file_num}.txt'
    k, kmers = parse_file(path)

    dbg = DeBruijinGraph(reads=kmers, k=k)
    path = dbg.eulerian_path()
    genome = path_to_genome(path)
    print(genome)

if __name__ == '__main__':
    main()
