import collections

def topological_sort(graph):
    """
    Topological sort of a DAG.

    Args:
        graph: A dictionary representing the graph where keys are nodes and values are lists of their neighbors.

    Returns:
        A list of nodes in topological order, or None if the graph is not a DAG.
    """
    in_degree = {node: 0 for node in graph}
    for node in graph:
        for neighbor in graph[node]:
            in_degree[neighbor] += 1

    queue = collections.deque([node for node in in_degree if in_degree[node] == 0])
    result = []

    while queue:
        node = queue.popleft()
        result.append(node)

        for neighbor in graph[node]:
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                queue.append(neighbor)

    if len(result) != len(graph):
        return None  # Graph has a cycle
    else:
        return result

def dag_shortest_paths(graph, weights, source):
    """
    Finds the shortest paths from a source node to all other nodes in a DAG.

    Args:
        graph: A dictionary representing the graph.
        weights: A dictionary of edge weights, where keys are tuples (u, v) and values are the weights.
        source: The source node.

    Returns:
        A dictionary of shortest path distances from the source to each node.
    """
    distances = {node: float('inf') for node in graph}
    distances[source] = 0
    topological_order = topological_sort(graph)

    for u in topological_order:
        for v in graph[u]:
            if distances[v] > distances[u] + weights[(u, v)]:
                distances[v] = distances[u] + weights[(u, v)]

    return distances

def dag_longest_paths(graph, weights, source):
    """
    Finds the longest paths from a source node to all other nodes in a DAG.

    Args:
        graph: A dictionary representing the graph.
        weights: A dictionary of edge weights.
        source: The source node.

    Returns:
        A dictionary of longest path distances from the source to each node.
    """
    distances = {node: float('-inf') for node in graph}
    distances[source] = 0
    topological_order = topological_sort(graph)

    for u in topological_order:
        for v in graph[u]:
            if distances[v] < distances[u] + weights[(u, v)]:
                distances[v] = distances[u] + weights[(u, v)]

    return distances

def all_paths(graph, source, target):
    """
    Finds all paths between two nodes in a DAG.

    Args:
        graph: A dictionary representing the graph.
        source: The source node.
        target: The target node.

    Returns:
        A list of all paths from the source to the target.
    """
    paths = []
    current_path = [source]

    def dfs(u, target, current_path, paths):
        if u == target:
            paths.append(current_path[:])  # Append a copy of the current path
            return

        for v in graph[u]:
            current_path.append(v)
            dfs(v, target, current_path, paths)
            current_path.pop()  # Backtrack

    dfs(source, target, current_path, paths)
    return paths

def resource_constrained_shortest_path(graph, weights, resources, source, target, max_resource):
    """
    Finds the resource-constrained shortest path in a DAG.

    Args:
        graph: A dictionary representing the graph.
        weights: A dictionary of edge weights.
        resources: A dictionary of edge resource consumption.
        source: The source node.
        target: The target node.
        max_resource: The maximum resource limit.

    Returns:
        The shortest path distance subject to the resource constraint, or None if no such path exists.
    """
    distances = {}
    predecessors = {}

    for node in graph:
        for resource_level in range(max_resource + 1):
            distances[(node, resource_level)] = float('inf')
            predecessors[(node, resource_level)] = None

    distances[(source, 0)] = 0

    topological_order = topological_sort(graph)

    for u in topological_order:
        for v in graph[u]:
            resource_cost = resources[(u, v)]
            for i in range(max_resource - resource_cost + 1):
                if distances[(u, i)] != float('inf'):
                    new_resource_level = i + resource_cost
                    if distances[(v, new_resource_level)] > distances[(u, i)] + weights[(u, v)]:
                        distances[(v, new_resource_level)] = distances[(u, i)] + weights[(u, v)]
                        predecessors[(v, new_resource_level)] = (u, i)

    min_dist = float('inf')
    best_resource_level = -1

    for i in range(max_resource + 1):
        if distances[(target, i)] < min_dist:
            min_dist = distances[(target, i)]
            best_resource_level = i

    if best_resource_level == -1:
        return None  # No feasible path
    else:
        # Reconstruct the path
        path = []
        curr = (target, best_resource_level)
        while curr:
            path.append(curr[0])
            curr = predecessors[curr]

        return min_dist, path[::-1]  # Return min distance and reversed path

def incremental_table_dag_shortest_path(graph, weights, source, target):
    """
    Finds the shortest path from source to target using incremental table building.

    Args:
        graph: A dictionary representing the graph.
        weights: A dictionary of edge weights.
        source: The source node.
        target: The target node.

    Returns:
        The shortest path distance from source to target.
    """
    distances = {node: float('inf') for node in graph}
    predecessors = {node: None for node in graph}
    visited = {node: False for node in graph}

    distances[source] = 0
    queue = collections.deque([source])

    while queue:
        u = queue.popleft()
        if u == target:
            break
        visited[u] = True

        for v in graph[u]:
            if distances[v] > distances[u] + weights[(u, v)]:
                distances[v] = distances[u] + weights[(u, v)]
                predecessors[v] = u
                if not visited[v]:
                    queue.append(v)

    return distances[target]
