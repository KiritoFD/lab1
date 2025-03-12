#ifndef DAG_ALGORITHMS_H
#define DAG_ALGORITHMS_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>

// 图的数据结构定义
typedef struct Edge {
    int dest;           // 目标顶点
    int weight;         // 边的权重
    int resource;       // 边的资源消耗（用于资源受限路径规划）
    struct Edge* next;  // 链表中的下一条边
} Edge;

typedef struct {
    int vertices;       // 顶点数量
    Edge** adjLists;    // 邻接表
} Graph;

// 用于存储路径的数据结构
typedef struct {
    int* vertices;
    int length;
    int capacity;
} Path;

typedef struct {
    Path** paths;
    int count;
    int capacity;
} PathList;

// 堆节点定义
typedef struct {
    int vertex;
    int distance;
} HeapNode;

// 最小堆数据结构
typedef struct {
    HeapNode* array;
    int size;
    int capacity;
    int* positions;  // 存储顶点在堆中的位置
} MinHeap;

// 最短路径结果结构
typedef struct {
    int distance;
    Path* path;
} ShortestPathResult;

// 资源受限最短路径结果结构
typedef struct {
    int distance;
    Path* path;
} ResourceConstrainedResult;

// 基本图操作
Graph* createGraph(int vertices);
void addEdge(Graph* graph, int src, int dest, int weight);
void addResourceEdge(Graph* graph, int src, int dest, int weight, int resource);
void freeGraph(Graph* graph);

// 路径操作
Path* createPath(int initialCapacity);
void addToPath(Path* path, int vertex);
void copyPath(Path* dest, Path* src);
void removeLastFromPath(Path* path);
void freePath(Path* path);
void printPath(Path* path);

// 路径列表操作
PathList* createPathList(int initialCapacity);
void addToPathList(PathList* list, Path* path);
void freePathList(PathList* list);
void printAllPaths(PathList* list);

// 最小堆操作
MinHeap* createMinHeap(int capacity, int vertices);
void swapHeapNodes(HeapNode* a, HeapNode* b);
void minHeapify(MinHeap* heap, int index);
bool isEmpty(MinHeap* heap);
HeapNode extractMin(MinHeap* heap);
void decreaseKey(MinHeap* heap, int vertex, int distance);
bool isInHeap(MinHeap* heap, int vertex);
void insert(MinHeap* heap, int vertex, int distance);
void freeMinHeap(MinHeap* heap);

// 拓扑排序
int* topologicalSort(Graph* graph);

// DAG最短路径
int* dagShortestPath(Graph* graph, int source);
Path* reconstructPath(int* predecessor, int source, int target);

// DAG最长路径
int* dagLongestPath(Graph* graph, int source);

// 所有路径枚举
PathList* allPaths(Graph* graph, int source, int target);
void dfsAllPaths(Graph* graph, int u, int target, Path* currentPath, PathList* pathList);

// 资源受限最短路径
ResourceConstrainedResult* resourceConstrainedShortestPath(Graph* graph, int source, int target, int resourceLimit);

// 增量式表格构建的最短路径
ShortestPathResult* incrementalTableDagShortestPath(Graph* graph, int source, int target);

#endif // DAG_ALGORITHMS_H
