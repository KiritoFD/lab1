# Path Planning Algorithms for Directed Acyclic Graphs (DAGs)

A Directed Acyclic Graph (DAG) is a graph with directed edges and no cycles. Due to their acyclic nature, DAGs are particularly suitable for path planning problems and can be solved with efficient algorithms. This document presents several algorithms for path planning in DAGs with standardized pseudocode following the conventions used in "Introduction to Algorithms" (CLRS).

## Problem Definition

The path planning problem in a DAG can be defined as:
- Given a DAG G = (V, E), where V is the set of vertices and E is the set of edges
- Each edge (u, v) may have a weight w(u, v)
- Given a source vertex s and a destination vertex t, find an "optimal" path from s to t

Depending on the optimization goal, we may seek the shortest path, longest path, or a path satisfying specific constraints.

## Topological Sorting: The Foundation

The key to efficient path planning in a DAG is utilizing the graph's topological order. A topological sort arranges all vertices in a linear sequence such that for every directed edge (u, v), vertex u appears before vertex v in the sequence.

### Kahn的拓扑排序算法

```
TOPOLOGICAL-SORT(G)
    // 输入：有向无环图 G = (V, E)
    // 输出：顶点的拓扑排序序列
    
    L ← 空列表，用于存储已排序的元素
    S ← 所有无入边的节点集合
    
    while S ≠ ∅ do
        从S中移除一个节点n
        L ← L ∪ {n}
        for 每个邻接点 m ∈ Adj[n] do
            从图中移除边 (n,m)
            if m的入度 = 0 then
                S ← S ∪ {m}
    
    if G中仍存在边 then
        return error // 图中存在环
    else
        return L
```

**时间复杂度**：O(|V| + |E|)，其中|V|是顶点数，|E|是边数。

## 1. DAG最短路径算法

### 算法描述

```
DAG-SHORTEST-PATH(G, w, s)
    // 输入：DAG图G = (V, E)，权重函数w，源顶点s
    // 输出：从s到所有顶点的最短路径距离
    
    对G中的顶点进行拓扑排序
    INITIALIZE-SINGLE-SOURCE(G, s)
    
    for 按拓扑排序顺序遍历每个顶点u do
        for 每个邻接点v ∈ Adj[u] do
            RELAX(u, v, w)
```

其中辅助过程为：

```
INITIALIZE-SINGLE-SOURCE(G, s)
    for G.V中的每个顶点v do
        d[v] ← ∞
        π[v] ← NIL
    d[s] ← 0
```

```
RELAX(u, v, w)
    if d[v] > d[u] + w(u, v) then
        d[v] ← d[u] + w(u, v)
        π[v] ← u
```

**时间复杂度**：O(|V| + |E|)，因为我们每个顶点和边都只处理一次。

## 2. DAG最长路径算法

在一般图中，寻找最长路径是NP难问题。但在DAG中，我们可以通过修改最短路径算法来高效地找到最长路径。

### 算法描述

```
DAG-LONGEST-PATH(G, w, s)
    // 输入：DAG图G = (V, E)，权重函数w，源顶点s
    // 输出：从s到所有顶点的最长路径距离
    
    对G中的顶点进行拓扑排序
    
    for G.V中的每个顶点v do
        d[v] ← -∞
        π[v] ← NIL
    d[s] ← 0
    
    for 按拓扑排序顺序遍历每个顶点u do
        for 每个邻接点v ∈ Adj[u] do
            if d[v] < d[u] + w(u, v) then
                d[v] ← d[u] + w(u, v)
                π[v] ← u
```

**时间复杂度**：O(|V| + |E|)，与最短路径算法相同。

## 3. DAG中的路径枚举

用于枚举两个顶点之间的所有路径：

```
ALL-PATHS(G, s, t)
    // 输入：DAG图G = (V, E)，源顶点s，目标顶点t
    // 输出：从s到t的所有路径
    
    paths ← 空列表
    current_path ← [s]
    DFS-ALL-PATHS(G, s, t, current_path, paths)
    return paths

DFS-ALL-PATHS(G, u, t, current_path, paths)
    if u = t then
        将current_path的副本添加到paths中
        return
    
    for 每个邻接点v ∈ Adj[u] do
        将v添加到current_path末尾
        DFS-ALL-PATHS(G, v, t, current_path, paths)
        从current_path中移除v
```

**时间复杂度**：在最坏情况下，如果从s到t的路径数量呈指数级增长，复杂度为O(|V|!)。但在实际中，这通常要低得多。

## 4. 资源受限的路径规划

用于寻找满足额外约束的路径，例如不超过资源限制：

```
RESOURCE-CONSTRAINED-SHORTEST-PATH(G, w, r, s, t, R)
    // 输入：DAG图G = (V, E)，权重函数w，资源函数r，
    //      源顶点s，目标顶点t，资源限制R
    // 输出：总资源消耗≤R的最短s-t路径
    
    对G中的顶点进行拓扑排序
    
    for G.V中的每个顶点v和每个资源等级i（从0到R）do
        d[v, i] ← ∞
        π[v, i] ← NIL
    d[s, 0] ← 0
    
    for 按拓扑排序顺序遍历每个顶点u do
        for 每个邻接点v ∈ Adj[u] do
            for 每个资源等级i（从0到R - r(u, v)）do
                if d[v, i + r(u, v)] > d[u, i] + w(u, v) then
                    d[v, i + r(u, v)] ← d[u, i] + w(u, v)
                    π[v, i + r(u, v)] ← (u, i)
    
    // 在所有可行资源等级中找到最小距离
    min_dist ← ∞
    best_res ← -1
    for i ← 0到R do
        if d[t, i] < min_dist then
            min_dist ← d[t, i]
            best_res ← i
    
    if best_res = -1 then
        return "不存在可行路径"
    else
        return CONSTRUCT-PATH(π, s, t, best_res)
```

**时间复杂度**：O(|V| + |E| × R)，其中R是资源限制。

## 5. 优化的表格构建

对于大型DAG，增量式表格构建可能更高效：

```
INCREMENTAL-TABLE-DAG-SHORTEST-PATH(G, w, s, t)
    // 输入：DAG图G = (V, E)，权重函数w，源顶点s，目标顶点t
    // 输出：从s到t的最短路径
    
    for G.V中的每个顶点v do
        d[v] ← ∞
        π[v] ← NIL
        visited[v] ← FALSE
    d[s] ← 0
    Q ← {s}
    
    while Q非空 do
        u ← EXTRACT-MIN(Q)
        if u = t then
            break
        visited[u] ← TRUE
        
        for 每个邻接点v ∈ Adj[u] do
            if d[v] > d[u] + w(u, v) then
                d[v] ← d[u] + w(u, v)
                π[v] ← u
                if 未访问v then
                    INSERT(Q, v)
    
    return d[t], CONSTRUCT-PATH(π, s, t)
```

**时间复杂度**：O(|V'| + |E'|)，其中|V'|和|E'|是实际处理的顶点和边的数量。

## Conclusion

The acyclic nature of DAGs enables many path planning problems to be solved efficiently in linear or near-linear time. The key algorithms we've presented are:

1. **DAG shortest path**: O(|V| + |E|)
2. **DAG longest path**: O(|V| + |E|)
3. **Path enumeration**: Worst case O(|V|!), typically much lower in practice
4. **Resource-constrained planning**: O(|V| + |E| × R)
5. **Incremental table building**: O(|V'| + |E'|)

These algorithms can be further optimized for specific applications by combining techniques or using parallel processing for large graphs.