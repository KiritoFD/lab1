#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>

typedef struct Edge {
    int dest;
    int weight;
    struct Edge* next;
} Edge;

typedef struct Graph {
    int vertices;
    Edge** adjLists;
} Graph;

typedef struct HeapNode {
    int vertex;
    int distance;
} HeapNode;

typedef struct MinHeap {
    int size;
    int capacity;
    int* positions;
    HeapNode* array;
} MinHeap;

typedef struct ShortestPathResult {
    int distance;
    int* path;
} ShortestPathResult;

MinHeap* createMinHeap(int capacity, int vertices) {
    MinHeap* heap = (MinHeap*)malloc(sizeof(MinHeap));
    heap->size = 0;
    heap->capacity = capacity;
    heap->positions = (int*)malloc(vertices * sizeof(int));
    heap->array = (HeapNode*)malloc(capacity * sizeof(HeapNode));
    for (int i = 0; i < vertices; i++) {
        heap->positions[i] = -1;
    }
    return heap;
}

void swapHeapNodes(HeapNode* a, HeapNode* b) {
    HeapNode temp = *a;
    *a = *b;
    *b = temp;
}

void decreaseKey(MinHeap* heap, int vertex, int distance) {
    // 获取顶点在堆中的索引
    int i = heap->positions[vertex];
    
    // 更新距离
    heap->array[i].distance = distance;
    
    // 自下而上调整堆
    while (i > 0 && heap->array[i].distance < heap->array[(i - 1) / 2].distance) {
        // 更新位置
        heap->positions[heap->array[i].vertex] = (i - 1) / 2;
        heap->positions[heap->array[(i - 1) / 2].vertex] = i;
        
        // 交换节点
        swapHeapNodes(&heap->array[i], &heap->array[(i - 1) / 2]);
        
        // 移动到父节点
        i = (i - 1) / 2;
    }
}

bool isInHeap(MinHeap* heap, int vertex) {
    return heap->positions[vertex] != -1;
}

void insert(MinHeap* heap, int vertex, int distance) {
    // 检查堆是否已满
    if (heap->size == heap->capacity) {
        fprintf(stderr, "堆已满\n");
        return;
    }
    
    // 增加堆大小
    heap->size++;
    int i = heap->size - 1;
    
    // 插入新节点
    heap->array[i].vertex = vertex;
    heap->array[i].distance = distance;
    heap->positions[vertex] = i;
    
    // 自下而上调整堆
    while (i > 0 && heap->array[i].distance < heap->array[(i - 1) / 2].distance) {
        // 更新位置
        heap->positions[heap->array[i].vertex] = (i - 1) / 2;
        heap->positions[heap->array[(i - 1) / 2].vertex] = i;
        
        // 交换节点
        swapHeapNodes(&heap->array[i], &heap->array[(i - 1) / 2]);
        
        // 移动到父节点
        i = (i - 1) / 2;
    }
}

void freeMinHeap(MinHeap* heap) {
    if (heap) {
        free(heap->array);
        free(heap->positions);
        free(heap);
    }
}

bool isEmpty(MinHeap* heap) {
    return heap->size == 0;
}

HeapNode extractMin(MinHeap* heap) {
    if (isEmpty(heap)) {
        HeapNode node = {-1, INT_MAX};
        return node;
    }
    
    // 存储根节点
    HeapNode root = heap->array[0];
    
    // 将最后一个节点移到根节点并减小堆大小
    HeapNode lastNode = heap->array[heap->size - 1];
    heap->array[0] = lastNode;
    
    // 更新位置
    heap->positions[root.vertex] = -1;
    heap->positions[lastNode.vertex] = 0;
    
    // 减小堆大小
    heap->size--;
    
    // 调整堆
    minHeapify(heap, 0);
    
    return root;
}

ShortestPathResult* incrementalTableDagShortestPath(Graph* graph, int source, int target) {
    int vertices = graph->vertices;
    
    // 分配距离和前驱数组
    int* distance = (int*)malloc(vertices * sizeof(int));
    int* predecessor = (int*)malloc(vertices * sizeof(int));
    bool* visited = (bool*)malloc(vertices * sizeof(bool));
    
    if (!distance || !predecessor || !visited) {
        fprintf(stderr, "内存分配失败\n");
        if (distance) free(distance);
        if (predecessor) free(predecessor);
        if (visited) free(visited);
        return NULL;
    }
    
    // 初始化
    for (int i = 0; i < vertices; i++) {
        distance[i] = INT_MAX;
        predecessor[i] = -1;
        visited[i] = false;
    }
    distance[source] = 0;
    
    // 创建最小堆
    MinHeap* heap = createMinHeap(vertices, vertices);
    
    // 插入源顶点到堆中
    insert(heap, source, 0);
    
    // 主循环
    while (!isEmpty(heap)) {
        // 提取最小距离顶点
        HeapNode minNode = extractMin(heap);
        int u = minNode.vertex;
        
        // 如果到达目标顶点，可以提前结束
        if (u == target) {
            break;
        }
        
        visited[u] = true;
        
        // 对所有邻接顶点进行松弛操作
        Edge* edge = graph->adjLists[u];
        while (edge) {
            int v = edge->dest;
            
            // 只有未访问的顶点才需要处理
            if (!visited[v]) {
                int weight = edge->weight;
                
                // 松弛操作
                if (distance[u] != INT_MAX && distance[v] > distance[u] + weight) {
                    distance[v] = distance[u] + weight;
                    predecessor[v] = u;
                    
                    // 更新堆
                    if (isInHeap(heap, v)) {
                        decreaseKey(heap, v, distance[v]);
                    } else {
                        insert(heap, v, distance[v]);
                    }
                }
            }
            
            edge = edge->next;
        }
    }
    
    ShortestPathResult* result = (ShortestPathResult*)malloc(sizeof(ShortestPathResult));
    if (!result) {
        fprintf(stderr, "内存分配失败\n");
        free(distance);
        free(predecessor);
        free(visited);
        freeMinHeap(heap);
        return NULL;
    }
    
    result->distance = distance[target];
    
    // 重建路径
    if (distance[target] == INT_MAX) {
        result->path = NULL;
    } else {
        result->path = reconstructPath(predecessor, source, target);
    }
    
    // 释放内存
    free(distance);
    free(predecessor);
    free(visited);
    freeMinHeap(heap);
    
    return result;
}
