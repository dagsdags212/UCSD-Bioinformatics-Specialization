from pydantic.dataclasses import dataclass
from pydantic import BaseModel


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

# create Node class
class Node:
    def __init__(self, seq: str):
        self.seq = seq
        self.k = len(seq)
        self.edges = []
        self.incoming = 0
        self.outgoing = 0

    def __hash__(self) -> int:
        return hash(self.seq)

    def __str__(self) -> str:
        return self.seq

    def __repr__(self) -> str:
        return self.seq

    def __eq__(self, other_node) -> bool:
        return self.seq == other_node.seq

    def __lt__(self, other_node) -> bool:
        return self.seq < other_node.seq

    def __gt__(self, other_node) -> bool:
        return self.seq > other_node.seq

    def get_suffix(self) -> str:
        """Return node suffix of length l."""
        return self.seq[1:]

    def get_prefix(self) -> str:
        """Return node prefix of length l."""
        return self.seq[:-1]

    def is_parent(self, other_node) -> bool:
        return self.get_suffix() == other_node.get_prefix()

    def is_child(self, other_node) -> bool:
        return self.get_prefix() == other_node.get_suffix()

    def is_neighbor(self, other_node) -> bool:
        if self.get_suffix() == other_node.get_prefix():
            print(f'{self.seq}->{other_node.seq}')
            return True
        if other_node.get_suffix() == self.get_prefix():
            print(f'{other_node.seq}->{self.seq}')
            return True
        return False

    def is_balanced(self) -> bool:
        return self.incoming == self.outgoing

# create Edge class
class DirectedEdge:
    def __init__(self, src: Node, dest: Node):
        if src.get_suffix() != dest.get_prefix():
            raise ValueError('Source node suffix must match the destination node prefix')
        self.src = src
        self.dest = dest
        self.src.outgoing += 1
        self.dest.incoming += 1
        self.src.edges.append(self)
        self.dest.edges.append(self)

    def __str__(self) -> str:
        return f'{self.src}->{self.dest}'

    def __repr__(self) -> str:
        return f'{self.src}->{self.dest}'

    def get_src(self) -> Node:
        return self.src

    def get_dest(self) -> Node:
        return self.dest

class OverlapGraph:
    def __init__(self, kmer_collection: list[str]):
        self.nodes = [Node(kmer) for kmer in kmer_collection]
        self.graph = self._construct_graph()

    def _construct_graph(self) -> dict:
        graph = {node.seq: [] for node in self.nodes}
        for n1 in self.nodes:
            for n2 in self.nodes:
                if n1 == n2:
                    continue
                if n2.is_child(n1):
                    graph[n1.seq].append(n2)
        # delete first and last nodes
        for n in list(graph.keys()):
            if graph[n] == []:
                del graph[n]
        return graph

    def get_nodes(self):
        return self.nodes

    def get_edges(self):
        edge_list = []
        for node, edges in self.graph.items():
            for edge in edges:
                if edge not in edge_list:
                    edge_list.append(edge)
        return edge_list

class DBG:
    def __init__(self, seq: str, k: int) -> None:
        self.seq = seq

        if isinstance(seq, str):
            self.k = k
            self.nodes: list[str] = self._construct_node_list()
            self.edges: list[str] = self._construct_edge_list()
            self.graph: dict = self._construct_graph()
        if isinstance(seq, list):
            self.k = len(seq[0])
            nodes = set([edge[:-1] for edge in seq] + [edge[1:] for edge in seq])
            self.nodes = sorted(list(map(Node, nodes)))
            self.edges = seq

        self.graph: dict = self._construct_graph()


    def __str__(self) -> str:
        """Format graph."""
        f = ''
        for k, v in self.graph.items():
            f += f"{k} :  {' '.join([node.seq for node in v])}\n"
        return f

    def _construct_node_list(self) -> list[str]:
        """Generate a unique (k-1)-mers of sequence to represent the glued nodes in graph."""
        nodes = []
        for i in range(len(self.seq)-self.k+2):
            node = Node(self.seq[i:i+self.k-1]) ##
            if node not in nodes:
                nodes.append(node)
        return sorted(nodes)

    def _construct_edge_list(self) -> list[str]:
        """Generate edges."""
        edges = []
        for i in range(len(self.seq)-self.k+1):
            edges.append(self.seq[i:i+self.k])
        return edges

    def _construct_graph(self) -> list[str]:
        """Construct a DeBruijin graph."""
        graph = {node: [] for node in self.nodes}
        for edge in self.edges:
            for node in self.nodes:
                if Node(edge[1:]) == node:
                    graph[Node(edge[:-1])].append(node)
        for k in graph.copy():
            if graph[k] == []:
                del graph[k]
        return graph

    def is_eulerian(self) -> bool:
        for node in self.nodes:
            if not node.is_balanced():
                return False
        return True

    def composition(self, k):
        pass

    def find_path(self):
        pass


