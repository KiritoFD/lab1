#include <iostream>
#include <vector>

// CUDA核函数，并行计算多个序列的滚动哈希值
__global__ void rollingHashKernel(const char* sequences, int num_sequences, int sequence_length,
                                  int prime, int mod, const int* base_map, int* hash_values) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num_sequences) {
        long long hash_val = 0;
        for (int i = 0; i < sequence_length; ++i) {
            hash_val = (hash_val * prime + base_map[sequences[idx * sequence_length + i]]) % mod;
        }
        hash_values[idx] = static_cast<int>(hash_val);
    }
}

// CPU函数，用于调用CUDA核函数
void calculateRollingHashes(const std::vector<std::string>& sequences, int prime, int mod,
                            const std::unordered_map<char, int>& base_map, std::vector<int>& hash_values) {
    int num_sequences = sequences.size();
    int sequence_length = sequences[0].length();

    // 将base_map转换为可以在CUDA中使用的int数组
    std::vector<int> cuda_base_map(256, 0);
    for (int i = 0; i < 256; ++i) {
        cuda_base_map[i] = base_map.count(static_cast<char>(i)) ? base_map.at(static_cast<char>(i)) : 0;
    }

    // 将序列数据复制到连续的内存区域
    std::vector<char> sequence_data;
    for (const auto& seq : sequences) {
        sequence_data.insert(sequence_data.end(), seq.begin(), seq.end());
    }

    // 在GPU上分配内存
    char* d_sequences;
    int* d_hash_values;
    int* d_base_map;
    cudaMalloc(&d_sequences, num_sequences * sequence_length * sizeof(char));
    cudaMalloc(&d_hash_values, num_sequences * sizeof(int));
    cudaMalloc(&d_base_map, 256 * sizeof(int));

    // 将数据复制到GPU
    cudaMemcpy(d_sequences, sequence_data.data(), num_sequences * sequence_length * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_hash_values, hash_values.data(), num_sequences * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_base_map, cuda_base_map.data(), 256 * sizeof(int), cudaMemcpyHostToDevice);

    // 配置CUDA核函数
    int threads_per_block = 256;
    int num_blocks = (num_sequences + threads_per_block - 1) / threads_per_block;

    // 调用CUDA核函数
    rollingHashKernel<<<num_blocks, threads_per_block>>>(d_sequences, num_sequences, sequence_length, prime, mod, d_base_map, d_hash_values);

    // 将结果复制回CPU
    cudaMemcpy(hash_values.data(), d_hash_values, num_sequences * sizeof(int), cudaMemcpyDeviceToHost);

    // 释放GPU内存
    cudaFree(d_sequences);
    cudaFree(d_hash_values);
    cudaFree(d_base_map);
}
