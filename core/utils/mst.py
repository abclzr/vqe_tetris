class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        root_x, root_y = self.find(x), self.find(y)
        if root_x == root_y:
            return False

        if self.rank[root_x] < self.rank[root_y]:
            self.parent[root_x] = root_y
        elif self.rank[root_x] > self.rank[root_y]:
            self.parent[root_y] = root_x
        else:
            self.parent[root_y] = root_x
            self.rank[root_x] += 1
        return True


def kruskal_mst(graph):
    def sort_edges(edges):
        return sorted(edges, key=lambda x: x[2])

    num_vertices = len(graph)
    union_find = UnionFind(num_vertices)
    mst_edges = []
    edges = []

    # Convert the graph to a list of edges (source, destination, weight)
    for i in range(num_vertices):
        for j in range(i + 1, num_vertices):
            if graph[i][j] != 0:
                edges.append((i, j, graph[i][j]))

    edges = sort_edges(edges)

    for edge in edges:
        src, dest, weight = edge
        if union_find.union(src, dest):
            mst_edges.append(edge)

    return mst_edges


# Example usage:
if __name__ == "__main__":
    # Example graph represented as an adjacency matrix
    graph = [
        [0, 2, 0, 6, 0],
        [2, 0, 3, 8, 5],
        [0, 3, 0, 0, 7],
        [6, 8, 0, 0, 9],
        [0, 5, 7, 9, 0]
    ]

    mst_edges = kruskal_mst(graph)
    total_weight = sum(weight for _, _, weight in mst_edges)
    print("Minimum Spanning Tree Edges:", mst_edges)
    print("Total Weight of MST:", total_weight)
