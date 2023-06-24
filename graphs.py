import itertools
import random
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
            count += ord(char) * (i+1)
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

    def __repr__(self):
        return f"{self.src}->{self.dest}"

    def __str__(self):
        return f"{self.src}->{self.dest}"

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
        self.g = self._generate_graph()
        self.edges = self._generate_edges()
        self.genome = None

    def _generate_graph(self) -> dict:
        graph = {node: [] for node in self.nodes}
        for parent in graph.keys():
            for child in self.nodes:
                if parent == child:
                    continue
                if parent.suffix == child.prefix:
                    graph[parent] += [child]
        return graph

    def _generate_edges(self) -> list[Edge]:
        edges = []
        for parent, children in self.g.items():
            for child in children:
                edges.append(Edge(parent, child))
        return edges

    def find_starting_node(self) -> Node | None:
        for node in self.nodes:
            if node.incoming < node.outgoing:
                return node

    def find_ending_node(self) -> Node | None:
        for node in self.nodes:
            if node.incoming > node.outgoing:
                return node

    def find_path(self, display: bool = False) -> list[Node]:
        # check if the graph has a Eulerian path
        if not self.has_eulerian_path():
            print('A Eulerian path does not exist!')
            return False

        # initialize stack with the starting node, and path with the ending node
        start = self.find_starting_node()
        end = self.find_ending_node()
        path = [end]
        stack = [start]
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
        return path

    def construct_genome(self):
        if not self.has_eulerian_path():
            return None
        return self.find_path()


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


kmers = ['CTTA', 'ACCA', 'TACC', 'GGCT', 'GCTT', 'TTAC']
k = 4
graph = DeBruijin(kmers, 4)

print(graph.find_path())
print(graph.construct_genome())

