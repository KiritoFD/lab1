# distutils: language = c++
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from libcpp.unordered_map cimport unordered_map as libcpp_unordered_map

cdef extern from "core_cuda.cu" namespace std:
    void calculateRollingHashes(libcpp_vector[libcpp_string]& sequences, int prime, int mod,
                                libcpp_unordered_map[char, int]& base_map, libcpp_vector[int]& hash_values)

def calculate_rolling_hashes(list sequences, int prime, int mod, dict base_map):
    cdef libcpp_vector[libcpp_string] cpp_sequences
    cdef libcpp_vector[int] cpp_hash_values
    cdef libcpp_unordered_map[char, int] cpp_base_map

    for seq in sequences:
        cpp_sequences.push_back(seq.encode('utf-8'))

    for key, value in base_map.items():
        cpp_base_map[key] = value

    cpp_hash_values.resize(len(sequences))
    calculateRollingHashes(cpp_sequences, prime, mod, cpp_base_map, cpp_hash_values)

    return [cpp_hash_values[i] for i in range(len(sequences))]
