import itertools
import random
from pathlib import Path
from typing import TypedDict, List, Optional

# ===== Eulerian Path Problem ===== #
class Node:
    def __init__(self, value) -> None:
        if isinstance(value, (int, str)):
            self._value = value
            self.incoming = 0
            self.outgoing = 0
            self.edges = []
        else:
            raise ValueError('Value must be of type int or str. Given value is of type ' + str(type(value)))

    def __hash__(self) -> int:
        count = 0
        for i, char in enumerate(str(self._value)):
            count += ord(char) * (i+5) ** 2
        return hash(count)

    def __str__(self) -> str:
        return str(self._value)

    def __repr__(self) -> str:
        return str(self._value)

    def __len__(self) -> None:
        pass

    def __add__(self, other: 'Node') -> 'Edge':
        return Edge(self, other)

    def __radd__(self, other: 'Node') -> 'Edge':
        return Edge(other, self)

    def __eq__(self, other: 'Node') -> bool:
        return hash(self) == hash(other)

    def __gt__(self, other: 'Node') -> None:
        pass

    def __lt__(self, other: 'Node') -> None:
        pass

    def __getitem__(self, idx) -> str:
        return str(self._value)[idx]

    def __lshift__(self, other: 'Node') -> bool:
        if isinstance(other.value(), str) and isinstance(self.value(), str):
            return self.prefix() == other.suffix()
        raise TypeError('Node object must be initialized with a string value')

    def __rshift__(self, other: 'Node') -> bool:
        if isinstance(other.value(), str) and isinstance(self.value(), str):
            return self.suffix() == other.prefix()
        raise TypeError('Node object must be initialized with a string value')

    def __contains__(self, other: list) -> bool:
        return self.value in other

    @property
    def value(self) -> str | int:
        return self._value

    @property
    def prefix(self) -> str:
        if isinstance(self._value, str):
            return self._value[:-1]
        raise TypeError('Node object must be initialized with a string value')

    @property
    def suffix(self) -> str:
        if isinstance(self._value, str):
            return self._value[1:]
        raise TypeError('Node object must be initialized with a string value')

    def is_balanced(self) -> bool:
        return self.incoming == self.outgoing

class Edge:
    def __init__(self, src: Node, dest: Node) -> None:
        src.outgoing += 1
        src.edges += [self]
        self._src = src

        dest.incoming += 1
        dest.edges += [self]
        self._dest = dest

    def __hash__(self) -> int:
        count = 0
        for i in self.src:
            count += ord(i)
        for j in self.dest:
            count += ord(j)
        return hash(count)

    def __repr__(self):
        return f"{self.src}->{self.dest}"

    def __str__(self):
        return f"{self.src}->{self.dest}"

    def __del__(self):
        self.src.outgoing -= 1
        self.dest.incoming -= 1

    def __eq__(self, other) -> bool:
        return hash(self) == hash(other)

    @property
    def src(self) -> Node:
        return self._src

    @property
    def dest(self) -> Node:
        return self._dest


class AdjacencyMatrix(TypedDict):
    Node: list[Node] | list

class DirectedGraph:
    def __init__(self, adjList: AdjacencyMatrix):
        self.adj_list = adjList
        self.nodes = self._generate_nodes()
        self.edges = self._generate_edges()
        self.g = self._generate_graph()
        self._balance_nodes()


    def _generate_nodes(self) -> list[Node]:
        parents = set(self.adj_list.keys())
        children = set([item for sl in self.adj_list.values() for item in sl])
        nodes = map(Node, (parents | children))
        return list(nodes)

    def _generate_edges(self)-> List[Edge]:
        edges = []
        for src, dests in self.adj_list.items():
            for dest in dests:
                edges += [Edge(Node(src), Node(dest))]
        return edges

    def _generate_graph(self) -> dict:
        graph = {}
        for src, dests in self.adj_list.items():
            parent = Node(src)
            parent.outgoing += len(dests)
            graph[Node(src)] = list(map(Node, dests))
        return graph

    def _balance_nodes(self, display=False) -> None:
        keys = list(self.g.keys())
        values = [d for dest in self.g.values() for d in dest]

        for src in self.g.keys():
            src.outgoing += len(self.g[src])
            src.incoming += values.count(src)

        for dests in self.g.values():
            for dest in dests:
                if dest in keys:
                    dest.incoming += len(self.g[dest])
                    dest.outgoing += values.count(dest)
                else:
                    dest.incoming += 1

        if display: print('Balancing done!')

    def find_path(self) -> List[Node]:
        def get_start():
            for src in self.g.keys():
                if src.outgoing > src.incoming:
                    return src

        def get_end():
            for dests in self.g.values():
                for dest in dests:
                    if dest.outgoing == 0:
                        return dest

        start, end = get_start(), get_end()

        path = [end]
        stack = [start]
        graph = self.g.copy()

        while stack:
            current = stack[-1]
            if len(graph[current]) > 0:
                child = graph[current].pop(0)
                if child != end:
                    stack.append(child)
            else:
                n = stack.pop()
                path = [n] + path
        print(' -> '.join([str(p.value) for p in path]))
        return path

keys = [0, 1, 2, 3, 6, 7, 8, 9]
values = [[2], [3], [1], [0,4], [3,7], [8], [9], [6]]
adjList = {k: v for k, v in zip(keys, values)}


class DeBruijin:
    def __init__(self, kmers, k):
        self.kmers = kmers
        self.k = k
        self.nodes = [Node(kmer) for kmer in kmers]
        self.edges = self._generate_edges()
        self.g = self._generate_graph()
        self._balance_nodes()
        self.path = self.find_path()
        self.genome = self.construct_genome()

    def __repr__(self) -> str:
        f = ""
        for key, value in self.g.items():
            line = ""
            line += f"{key}: "
            for v in value:
                line += f"{v} "
            f += line + "\n"
        return f


    def _generate_edges(self) -> list[Edge]:
        edges = []
        for n1 in self.nodes:
            for n2 in self.nodes:
                if n1 == n2: continue
                if n1.suffix == n2.prefix:
                    e = Edge(n1, n2)
                    if e not in edges:
                        edges.append(e)
                    else:
                        del e
        return edges

    def _generate_graph(self) -> dict:
        graph = {}
        for edge in self.edges:
            if edge.src not in graph:
                graph[edge.src] = [edge.dest]
            else:
                graph[edge.src].append(edge.dest)
        return graph

    def _balance_nodes(self) -> None:
        values = [i for sl in self.g.values() for i in sl]

        for key in self.g.keys():
            key.outgoing = len(self.g[key])
            key.incoming = values.count(key)

        for values in self.g.values():
            for value in values:
                value.outgoing = len(self.g[value]) if value in self.g else 0
                value.incoming = values.count(value)

    def find_starting_node(self) -> Node | None:
        for node in self.g.keys():
            if node.incoming < node.outgoing:
                return node

    def find_ending_node(self) -> Node | None:
        for nodes in self.g.values():
            for node in nodes:
                if node.incoming > node.outgoing:
                    return node

    def find_path(self, display: bool = False) -> list[Node]:
        # initialize stack with the starting node, and path with the ending node
        start = self.find_starting_node()
        end = self.find_ending_node()
        stack = [start]
        path = [end]
        graph = self.g.copy()

        while stack:
            current = stack[-1]
            if len(graph[current]) > 0:
                child = graph[current].pop()
                if child != end:
                    stack.append(child)
            else:
                n = stack.pop()
                path = [n] + path
        # display eulerian path
        if display:
            print('->'.join([str(p.value) for p in path]))
        self.path = path
        return path

    def construct_genome(self) -> str | None:
        if self.path:
            genome = self.path[0].value
            for node in self.path[1:]:
                genome += node[-1]
            return genome

    def has_eulerian_path(self) -> bool:
        unbalanced_count = 0
        for node in self.nodes:
            if not node.is_balanced():
                unbalanced_count += 1
                if unbalanced_count > 2:
                    return False
        return True if unbalanced_count <= 2 else False

    def has_eulerian_cycle(self) -> bool:
        for node in self.nodes:
            if not node.is_balanced():
                return False
        return True


num = 8
filename = f"input_{num}.txt"
directory = 'datasets/string_reconstruction'
PATH = Path.cwd() / directory / filename

def parse_file(path):
    with open(path, 'r') as fh:
        k = fh.readline().strip()
        kmers = fh.readline().strip().split()
    fh.close()
    return k, kmers

k, kmers = parse_file(PATH)
graph = DeBruijin(kmers, 4)

print(graph.edges)
