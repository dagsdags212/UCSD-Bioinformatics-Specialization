from deBruijn import *

def maximal_nonbranching_paths(graph):
    paths = []
    for node in graph.nodes:
        if not node.is_balanced():
            if node.outgoing > 0:
                for child in graph[node]:
                    nonbranching_path = [node]
                    while child.is_balanced():
                        nonbranching_path.append(node)
                        child = graph[child]
                    paths.append(nonbranching_path)
    return paths
