#include <stdio.h>
#include <stdlib.h>
#include "dag_algorithms.h"

// 打印测试结果的辅助函数
void printDistances(int* distances, int vertices, const char* label) {
    printf("\n%s:\n", label);
    for (int i = 0; i < vertices; i++) {
        if (distances[i] == INT_MAX) {
            printf("顶点 %d: 不可达\n", i);
        } else {
            printf("顶点 %d: %d\n", i, distances[i]);
        }
    }
}

// 测试函数
void testBasicDAGAlgorithms() {
    printf("==== 基本DAG算法测试 ====\n");
    
    // 创建一个简单的DAG
    // 0 -> 1 -> 3
    // |    |
    // v    v
    // 2 -> 4
    int vertices = 5;
    Graph* g = createGraph(vertices);
    
    // 添加边和权重
    addEdge(g, 0, 1, 2);  // 0->1, 权重2
    addEdge(g, 0, 2, 3);  // 0->2, 权重3
    addEdge(g, 1, 3, 5);  // 1->3, 权重5
    addEdge(g, 1, 4, 1);  // 1->4, 权重1
    addEdge(g, 2, 4, 6);  // 2->4, 权重6
    
    printf("图已创建，具有 %d 个顶点和以下边：\n", vertices);
    printf("0->1 (2), 0->2 (3), 1->3 (5), 1->4 (1), 2->4 (6)\n\n");
    
    // 拓扑排序
    int* topoOrder = topologicalSort(g);
    if (topoOrder) {
        printf("拓扑排序结果: ");
        for (int i = 0; i < vertices; i++) {
            printf("%d ", topoOrder[i]);
        }
        printf("\n");
        free(topoOrder);
    } else {
        printf("图中存在环，不是DAG\n");
    }
    
    // 最短路径
    int source = 0;
    int* shortestDistances = dagShortestPath(g, source);
    
    if (shortestDistances) {
        printDistances(shortestDistances, vertices, "从顶点0开始的最短路径距离");
        free(shortestDistances);
    }
    
    // 最长路径
    int* longestDistances = dagLongestPath(g, source);
    
    if (longestDistances) {
        printDistances(longestDistances, vertices, "从顶点0开始的最长路径距离");
        free(longestDistances);
    }
    
    // 清理图
    freeGraph(g);
}

void testAllPaths() {
    printf("\n==== 所有路径枚举测试 ====\n");
    
    // 创建一个简单的DAG
    int vertices = 5;
    Graph* g = createGraph(vertices);
    
    // 添加边
    addEdge(g, 0, 1, 1);
    addEdge(g, 0, 2, 1);
    addEdge(g, 1, 3, 1);
    addEdge(g, 1, 4, 1);
    addEdge(g, 2, 3, 1);
    addEdge(g, 2, 4, 1);
    
    printf("图已创建，具有多条从顶点0到顶点4的路径\n");
    
    // 枚举所有路径
    int source = 0;
    int target = 4;
    
    PathList* paths = allPaths(g, source, target);
    
    printf("从顶点 %d 到顶点 %d 的所有路径：\n", source, target);
    printAllPaths(paths);
    
    // 清理
    freePathList(paths);
    freeGraph(g);
}

void testResourceConstrainedPath() {
    printf("\n==== 资源受限路径规划测试 ====\n");
    
    // 创建一个DAG，边有权重和资源消耗
    int vertices = 5;
    Graph* g = createGraph(vertices);
    
    // 添加边：(源，目标，权重，资源消耗)
    addResourceEdge(g, 0, 1, 2, 3);  // 0->1: 权重2，资源3
    addResourceEdge(g, 0, 2, 1, 5);  // 0->2: 权重1，资源5
    addResourceEdge(g, 1, 3, 3, 2);  // 1->3: 权重3，资源2
    addResourceEdge(g, 1, 4, 5, 1);  // 1->4: 权重5，资源1
    addResourceEdge(g, 2, 3, 2, 1);  // 2->3: 权重2，资源1
    addResourceEdge(g, 2, 4, 4, 3);  // 2->4: 权重4，资源3
    addResourceEdge(g, 3, 4, 1, 2);  // 3->4: 权重1，资源2
    
    printf("图已创建，边具有权重和资源消耗\n");
    
    int source = 0;
    int target = 4;
    int resourceLimit = 7;
    
    printf("源顶点: %d, 目标顶点: %d, 资源限制: %d\n", source, target, resourceLimit);
    
    ResourceConstrainedResult* result = resourceConstrainedShortestPath(g, source, target, resourceLimit);
    
    if (result) {
        if (result->path) {
            printf("找到满足资源限制的最短路径，总距离: %d\n", result->distance);
            printPath(result->path);
            freePath(result->path);
        } else {
            printf("在资源限制 %d 下，没有找到从顶点 %d 到顶点 %d 的路径\n", 
                   resourceLimit, source, target);
        }
        free(result);
    }
    
    // 清理
    freeGraph(g);
}

void testIncrementalTableMethod() {
    printf("\n==== 增量式表格构建最短路径测试 ====\n");
    
    // 创建一个大一点的DAG
    int vertices = 8;
    Graph* g = createGraph(vertices);
    
    // 添加边
    addEdge(g, 0, 1, 3);
    addEdge(g, 0, 2, 1);
    addEdge(g, 1, 3, 2);
    addEdge(g, 1, 4, 4);
    addEdge(g, 2, 3, 5);
    addEdge(g, 2, 5, 2);
    addEdge(g, 3, 4, 1);
    addEdge(g, 3, 6, 7);
    addEdge(g, 4, 7, 3);
    addEdge(g, 5, 6, 4);
    addEdge(g, 6, 7, 1);
    
    printf("创建了一个较大的DAG，测试增量式表格构建方法\n");
    
    int source = 0;
    int target = 7;
    
    ShortestPathResult* result = incrementalTableDagShortestPath(g, source, target);
    
    if (result) {
        if (result->path) {
            printf("从顶点 %d 到顶点 %d 的最短路径距离: %d\n", source, target, result->distance);
            printPath(result->path);
            freePath(result->path);
        } else {
            printf("无法到达顶点 %d\n", target);
        }
        free(result);
    }
    
    // 清理
    freeGraph(g);
}

int main() {
    printf("DAG路径规划算法演示\n");
    printf("====================\n\n");
    
    testBasicDAGAlgorithms();
    testAllPaths();
    testResourceConstrainedPath();
    testIncrementalTableMethod();
    
    printf("\n所有测试完成!\n");
    
    return 0;
}
