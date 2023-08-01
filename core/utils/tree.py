import queue

class TreeNode:
    def __init__(self, idx, parent):
        self.idx = idx
        self.parent = parent

class Tree:
    def __init__(self, edges, root):
        self.node_list = []
        self.root = root
        self.flood_fill(root, edges)
        # the node_list have all nodes in the order of distance from root from low to high
        # reverse the node_list so it's convenient to run the Pauli String
        self.node_list.reverse()
    
    def flood_fill(self, root, edges):
        visited = [root]
        self.node_list.append(TreeNode(root, -1))
        Q = queue.Queue()
        Q.put(root)
        while not Q.empty():
            node = Q.get()
            for edge in edges:
                if node in edge:
                    adj = edge[0] + edge[1] - node
                    if not adj in visited:
                        visited.append(adj)
                        self.node_list.append(TreeNode(adj, node))
                        Q.put(adj)
                        