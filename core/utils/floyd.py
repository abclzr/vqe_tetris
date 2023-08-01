import pdb

def floyd_warshall(graph):
    """
    Finds the shortest paths between all pairs of vertices in a weighted graph using the Floyd-Warshall algorithm.
    
    Parameters:
        graph (list of lists): The graph representation as an adjacency matrix. The weight of edge (u, v) is stored
                               in graph[u][v]. If there's no edge between u and v, the weight should be set to infinity
                               (or a sufficiently large value).
                               
    Returns:
        list of lists: A 2D matrix containing the shortest distances between all pairs of vertices.
    """
    num_vertices = len(graph)
    
    # Initialize the distance matrix with the same values as the adjacency matrix.
    dist = [[1 if x == 1 else 0x777777 for x in row] for row in graph]
    for i in range(num_vertices):
        dist[i][i] = 0
    
    # Main algorithm
    for k in range(num_vertices):
        for i in range(num_vertices):
            for j in range(num_vertices):
                # If there's a shorter path from vertex i to vertex j through vertex k, update the distance.
                dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j])
    
    return dist

