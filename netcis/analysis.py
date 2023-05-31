


def graph_properties(G):
    print(f"number of nodes: {G.number_of_nodes()}")
    print(f"number of edges: {G.number_of_edges()}")
    num_inserts = 0
    for node in G.nodes:
        num_inserts += G.nodes[node]['counts']
    print(f"number of insertions: {num_inserts}")
    
    
