#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "dag_algorithms.h"

// 辅助函数：反转字符串
void reverseString(char* str) {
    int len = strlen(str);
    for (int i = 0; i < len / 2; i++) {
        char temp = str[i];
        str[i] = str[len - i - 1];
        str[len - i - 1] = temp;
    }
}

// 辅助函数：检查重复序列
void findRepeats(const char* query, const char* reference, FILE* outputFile) {
    int queryLen = strlen(query);
    int refLen = strlen(reference);

    for (int len = 50; len <= 100; len++) { // 检查长度为50到100的重复序列
        for (int i = 0; i <= queryLen - len; i++) {
            char subQuery[101];
            strncpy(subQuery, query + i, len);
            subQuery[len] = '\0';

            // 检查正向重复
            for (int j = 0; j <= refLen - len; j++) {
                char subRef[101];
                strncpy(subRef, reference + j, len);
                subRef[len] = '\0';

                if (strcmp(subQuery, subRef) == 0) {
                    fprintf(outputFile, "%d,%d,1,否,%s,%s\n", i + 1, len, subQuery, subRef);
                }
            }

            // 检查反向重复
            char revSubQuery[101];
            strncpy(revSubQuery, subQuery, len);
            revSubQuery[len] = '\0';
            reverseString(revSubQuery);

            for (int j = 0; j <= refLen - len; j++) {
                char subRef[101];
                strncpy(subRef, reference + j, len);
                subRef[len] = '\0';

                if (strcmp(revSubQuery, subRef) == 0) {
                    fprintf(outputFile, "%d,%d,1,是,%s,%s\n", i + 1, len, subQuery, subRef);
                }
            }
        }
    }
}

// 辅助函数：从文件读取图数据
Graph* readGraphFromFile(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "无法打开文件: %s\n", filename);
        return NULL;
    }

    int vertices, edges;
    if (fscanf(file, "%d %d", &vertices, &edges) != 2) {
        fprintf(stderr, "文件格式错误: 顶点和边数\n");
        fclose(file);
        return NULL;
    }

    Graph* graph = createGraph(vertices);
    if (!graph) {
        fclose(file);
        return NULL;
    }

    for (int i = 0; i < edges; i++) {
        int src, dest, weight, resource;
        if (fscanf(file, "%d %d %d %d", &src, &dest, &weight, &resource) == 4) {
            addResourceEdge(graph, src, dest, weight, resource);
        } else if (fscanf(file, "%d %d %d", &src, &dest, &weight) == 3) {
            addEdge(graph, src, dest, weight);
        } else {
            fprintf(stderr, "文件格式错误: 边数据\n");
            freeGraph(graph);
            fclose(file);
            return NULL;
        }
    }

    fclose(file);
    return graph;
}

// 辅助函数：将路径信息写入文件
void writePathToFile(FILE* file, Path* path) {
    if (!path) {
        fprintf(file, "路径: NULL\n");
        return;
    }

    fprintf(file, "路径: 长度 = %d, 序列 = [", path->length);
    for (int i = 0; i < path->length; i++) {
        fprintf(file, "%d", path->vertices[i]);
        if (i < path->length - 1) {
            fprintf(file, ", ");
        }
    }
    fprintf(file, "]\n");
}

// 辅助函数：将所有路径信息写入文件
void writeAllPathsToFile(FILE* file, PathList* paths) {
    fprintf(file, "所有路径 (共 %d 条):\n", paths->count);
    for (int i = 0; i < paths->count; i++) {
        fprintf(file, "%d. ", i + 1);
        writePathToFile(file, paths->paths[i]);
    }
}

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
void testBasicDAGAlgorithms(Graph* graph) {
    printf("==== 基本DAG算法测试 ====\n");
    
    int vertices = graph->vertices;
    
    printf("图已创建，具有 %d 个顶点\n", vertices);
    
    // 拓扑排序
    int* topoOrder = topologicalSort(graph);
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
    int* shortestDistances = dagShortestPath(graph, source);
    
    if (shortestDistances) {
        printDistances(shortestDistances, vertices, "从顶点0开始的最短路径距离");
        free(shortestDistances);
    }
    
    // 最长路径
    int* longestDistances = dagLongestPath(graph, source);
    
    if (longestDistances) {
        printDistances(longestDistances, vertices, "从顶点0开始的最长路径距离");
        free(longestDistances);
    }
}

void testBasicDAGAlgorithmsToFile(FILE* file, Graph* graph) {
    fprintf(file, "==== 基本DAG算法测试 ====\n");
    
    int vertices = graph->vertices;
    
    fprintf(file, "图已创建，具有 %d 个顶点\n", vertices);
    
    // 拓扑排序
    int* topoOrder = topologicalSort(graph);
    if (topoOrder) {
        fprintf(file, "拓扑排序结果: ");
        for (int i = 0; i < vertices; i++) {
            fprintf(file, "%d ", topoOrder[i]);
        }
        fprintf(file, "\n");
        free(topoOrder);
    } else {
        fprintf(file, "图中存在环，不是DAG\n");
    }
    
    // 最短路径
    int source = 0;
    int* shortestDistances = dagShortestPath(graph, source);
    
    if (shortestDistances) {
        fprintf(file, "从顶点0开始的最短路径距离:\n");
        for (int i = 0; i < vertices; i++) {
            if (shortestDistances[i] == INT_MAX) {
                fprintf(file, "顶点 %d: 不可达\n", i);
            } else {
                fprintf(file, "顶点 %d: %d\n", i, shortestDistances[i]);
            }
        }
        free(shortestDistances);
    }
    
    // 最长路径
    int* longestDistances = dagLongestPath(graph, source);
    
    if (longestDistances) {
        fprintf(file, "从顶点0开始的最长路径距离:\n");
        for (int i = 0; i < vertices; i++) {
            if (longestDistances[i] == INT_MIN) {
                fprintf(file, "顶点 %d: 不可达\n", i);
            } else {
                fprintf(file, "顶点 %d: %d\n", i, longestDistances[i]);
            }
        }
        free(longestDistances);
    }
}

void testAllPaths(Graph* graph) {
    printf("\n==== 所有路径枚举测试 ====\n");
    
    int vertices = graph->vertices;
    
    printf("图已创建，具有 %d 个顶点\n", vertices);
    
    // 枚举所有路径
    int source = 0;
    int target = vertices - 1;
    
    PathList* paths = allPaths(graph, source, target);
    
    printf("从顶点 %d 到顶点 %d 的所有路径：\n", source, target);
    printAllPaths(paths);
    
    // 清理
    freePathList(paths);
}

void testAllPathsToFile(FILE* file, Graph* graph) {
    fprintf(file, "\n==== 所有路径枚举测试 ====\n");
    
    int vertices = graph->vertices;
    
    fprintf(file, "图已创建，具有 %d 个顶点\n", vertices);
    
    // 枚举所有路径
    int source = 0;
    int target = vertices - 1;
    
    PathList* paths = allPaths(graph, source, target);
    
    fprintf(file, "从顶点 %d 到顶点 %d 的所有路径：\n", source, target);
    writeAllPathsToFile(file, paths);
    
    // 清理
    freePathList(paths);
}

void testResourceConstrainedPath(Graph* graph) {
    printf("\n==== 资源受限路径规划测试 ====\n");
    
    int vertices = graph->vertices;
    
    printf("图已创建，具有 %d 个顶点\n", vertices);
    
    int source = 0;
    int target = vertices - 1;
    int resourceLimit = 10;
    
    printf("源顶点: %d, 目标顶点: %d, 资源限制: %d\n", source, target, resourceLimit);
    
    ResourceConstrainedResult* result = resourceConstrainedShortestPath(graph, source, target, resourceLimit);
    
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
}

void testResourceConstrainedPathToFile(FILE* file, Graph* graph) {
    fprintf(file, "\n==== 资源受限路径规划测试 ====\n");
    
    int vertices = graph->vertices;
    
    fprintf(file, "图已创建，具有 %d 个顶点\n", vertices);
    
    int source = 0;
    int target = vertices - 1;
    int resourceLimit = 10;
    
    fprintf(file, "源顶点: %d, 目标顶点: %d, 资源限制: %d\n", source, target, resourceLimit);
    
    ResourceConstrainedResult* result = resourceConstrainedShortestPath(graph, source, target, resourceLimit);
    
    if (result) {
        if (result->path) {
            fprintf(file, "找到满足资源限制的最短路径，总距离: %d\n", result->distance);
            writePathToFile(file, result->path);
            freePath(result->path);
        } else {
            fprintf(file, "在资源限制 %d 下，没有找到从顶点 %d 到顶点 %d 的路径\n", 
                   resourceLimit, source, target);
        }
        free(result);
    }
}

void testIncrementalTableMethod(Graph* graph) {
    printf("\n==== 增量式表格构建最短路径测试 ====\n");
    
    int vertices = graph->vertices;
    
    printf("图已创建，具有 %d 个顶点\n", vertices);
    
    int source = 0;
    int target = vertices - 1;
    
    ShortestPathResult* result = incrementalTableDagShortestPath(graph, source, target);
    
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
}

void testIncrementalTableMethodToFile(FILE* file, Graph* graph) {
    fprintf(file, "\n==== 增量式表格构建最短路径测试 ====\n");
    
    int vertices = graph->vertices;
    
    fprintf(file, "图已创建，具有 %d 个顶点\n", vertices);
    
    int source = 0;
    int target = vertices - 1;
    
    ShortestPathResult* result = incrementalTableDagShortestPath(graph, source, target);
    
    if (result) {
        if (result->path) {
            fprintf(file, "从顶点 %d 到顶点 %d 的最短路径距离: %d\n", source, target, result->distance);
            writePathToFile(file, result->path);
            freePath(result->path);
        } else {
            fprintf(file, "无法到达顶点 %d\n", target);
        }
        free(result);
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "用法: %s <query_file> <output_file>\n", argv[0]);
        return 1;
    }

    const char* queryFile = argv[1];
    const char* outputFile = argv[2];

    // 读取query文件
    FILE* queryF = fopen(queryFile, "r");
    if (!queryF) {
        fprintf(stderr, "无法打开query文件: %s\n", queryFile);
        return 1;
    }
    char* query = NULL;
    size_t queryLen = 0;
    fseek(queryF, 0, SEEK_END);
    long fsize = ftell(queryF);
    fseek(queryF, 0, SEEK_SET);

    query = (char*)malloc(fsize + 1);
    fread(query, fsize, 1, queryF);
    query[fsize] = '\0';
    fclose(queryF);

    // 读取reference文件
    FILE* referenceF = fopen("/home/xy/lab1/reference.txt", "r");
    if (!referenceF) {
        fprintf(stderr, "无法打开reference文件: /home/xy/lab1/reference.txt\n");
        free(query);
        return 1;
    }
    char* reference = NULL;
    size_t referenceLen = 0;
    fseek(referenceF, 0, SEEK_END);
    long refsize = ftell(referenceF);
    fseek(referenceF, 0, SEEK_SET);

    reference = (char*)malloc(refsize + 1);
    fread(reference, refsize, 1, referenceF);
    reference[refsize] = '\0';
    fclose(referenceF);

    // 打开输出文件
    FILE* outputF = fopen(outputFile, "w");
    if (!outputF) {
        fprintf(stderr, "无法打开输出文件: %s\n", outputFile);
        free(query);
        free(reference);
        return 1;
    }

    // 写入表头
    fprintf(outputF, "位置,长度,重复次数,是否反向重复,原始序列,重复实例\n");

    // 查找并写入重复序列
    findRepeats(query, reference, outputF);

    // 清理
    free(query);
    free(reference);
    fclose(outputF);

    printf("结果已写入文件: %s\n", outputFile);

    return 0;
}
