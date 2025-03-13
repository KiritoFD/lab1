# DNA Repeat Finder

A high-performance application for finding repeat patterns in DNA sequences, optimized for AMD Ryzen 9 7940HX processors.

## Project Structure

```
/home/xy/lab1/
├── src/               # Core C implementation
│   ├── main.c         # Main program entry point
│   ├── dna_common.c   # Common DNA utilities
│   ├── dna_io.c       # I/O operations for DNA sequences
│   ├── dna_traditional.c # Traditional repeat finding implementation
│   ├── dna_graph.c    # Graph-based repeat finding implementation
│   └── dna_sequence_utils.c # Sequence utilities
│
├── include/           # Header files
│   ├── dna_common.h   # Common DNA utilities declarations
│   ├── dna_io.h       # I/O operations declarations
│   ├── dna_traditional.h # Traditional approach declarations
│   ├── dna_graph.h    # Graph-based approach declarations
│   └── cpu_optimize.h # R9 7940HX specific optimizations
│
├── tools/            # Additional tools
│   └── generate_test_sequences.c # Test sequence generator (optional)
│
├── algorithms/       # DAG algorithm implementations (optional)
│   ├── dag_algorithms.c
│   ├── dag_algorithms.h
│   └── dag_path_planning.md
│
├── cpp/             # C++ implementation (alternative)
│   └── dna_repeat_finder.cpp
│
├── bin/             # Executable files (created during build)
├── obj/             # Object files (created during build)
├── Makefile         # Build configuration
└── README.md        # This file
```

## Building the Project

```bash
# Standard build
make

# Debug build
make debug

# Clean and rebuild
make clean && make
```

## Running the Program

```bash
# Using default files (reference.txt and query.txt)
./bin/dna_repeat_finder

# Specifying input files
./bin/dna_repeat_finder path/to/reference.txt path/to/query.txt
```
![alt text]({25DD791A-D1DF-4291-A9DD-2407276C24AD}.png)

## 算法描述与分析

### DNA重复片段查找算法

本实现用于查找满足以下条件的DNA重复片段：
1. 在参考序列(reference)中只出现一次，但在查询序列(query)中连续重复出现的片段
2. 在参考序列中出现的序列，在查询序列中出现其反向互补序列

#### 一、基础算法（时间复杂度 O(|R|² × |Q|)）

##### 算法伪代码（算法导论格式）

```
FIND-DNA-REPEATS(query, reference, min_length)
    max_length = 120
    repeats = []
    ref_len = LENGTH(reference)
    query_len = LENGTH(query)
    
    // 在参考序列中枚举所有可能的起点位置
    for pos = 0 to ref_len - min_length - 1
        // 枚举所有可能的片段长度
        for length = min_length to min(max_length, ref_len - pos)
            // 从参考序列中提取片段
            segment = reference[pos..(pos+length-1)]
            
            // 检查正向重复
            repeat_info = CHECK-CONSECUTIVE-REPEATS(segment, query, pos, FALSE, reference)
            if repeat_info != NIL and repeat_info.count > 1
                repeats.APPEND({
                    position: pos,
                    length: length,
                    repeat_count: repeat_info.count,
                    is_reverse: FALSE,
                    original_sequence: segment,
                    query_position: repeat_info.position
                })
            
            // 检查反向互补重复
            rev_comp = GET-REVERSE-COMPLEMENT(segment)
            repeat_info = CHECK-CONSECUTIVE-REPEATS(rev_comp, query, pos, TRUE, reference)
            if repeat_info != NIL and repeat_info.count > 1
                repeats.APPEND({
                    position: pos,
                    length: length,
                    repeat_count: repeat_info.count,
                    is_reverse: TRUE,
                    original_sequence: segment,
                    query_position: repeat_info.position
                })
    
    // 过滤嵌套重复
    filtered_repeats = FILTER-NESTED-REPEATS(repeats)
    
    // 按重复长度和次数排序
    SORT(filtered_repeats) by (length * repeat_count) in descending order
    
    return filtered_repeats

CHECK-CONSECUTIVE-REPEATS(segment, query, ref_pos, is_reverse, reference)
    query_len = LENGTH(query)
    seg_len = LENGTH(segment)
    
    // 在查询序列中找第一个匹配
    start_idx = 0
    while TRUE
        idx = FIND-SUBSTRING(query, segment, start_idx)
        if idx == -1
            return NIL
            
        // 检查是否有连续重复
        count = 1
        current_pos = idx + seg_len
        
        while current_pos + seg_len <= query_len
            if query[current_pos..(current_pos+seg_len-1)] == segment
                count = count + 1
                current_pos = current_pos + seg_len
            else
                break
        
        // 验证在参考序列中只出现一次
        if count > 1
            if is_reverse
                orig_seq = GET-REVERSE-COMPLEMENT(segment)
                if COUNT-OCCURRENCES(reference, orig_seq) <= 1
                    return {position: idx, count: count}
            else
                if COUNT-OCCURRENCES(reference, segment) <= 1
                    return {position: idx, count: count}
        
        // 继续搜索下一个位置
        start_idx = idx + 1
        
        if start_idx >= query_len
            return NIL

FILTER-NESTED-REPEATS(repeats)
    if repeats is empty
        return []
    
    // 按位置和反向标志分组
    position_groups = {}
    for repeat in repeats
        key = (repeat.position, repeat.is_reverse)
        if key not in position_groups
            position_groups[key] = []
        position_groups[key].APPEND(repeat)
    
    // 从每组中选择最长的重复
    filtered_repeats = []
    for group in position_groups.values()
        best_repeat = MAX(group, key = length)
        filtered_repeats.APPEND(best_repeat)
    
    return filtered_repeats

GET-REVERSE-COMPLEMENT(dna)
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    reversed_dna = REVERSE(dna)
    result = ""
    for base in reversed_dna
        result = result + complement[base]
    return result
```

##### 基础算法复杂度分析

**时间复杂度**：O(|R|² × |Q|)
- 枚举所有可能的参考序列片段：O(|R|²)
- 对每个片段，在查询序列中检查连续重复：O(|Q|)

**空间复杂度**：O(N + max_length)，其中N是找到的重复片段数量

#### 二、优化算法（时间复杂度 O(|R| × |Q|)）

##### 优化思路

基础算法的主要性能瓶颈在于双重循环枚举所有可能的参考序列片段。优化的关键是使用索引和种子扩展技术，避免穷举所有可能的片段长度。

##### 优化算法伪代码（算法导论格式）

```
FIND-DNA-REPEATS-OPTIMIZED(query, reference, min_length)
    // 步骤1: 为参考序列创建索引 - O(|R|)
    segment_positions = {}
    segment_counts = {}
    
    // 为每个长度为min_length的子串建立索引
    for pos = 0 to LENGTH(reference) - min_length
        segment = reference[pos..(pos+min_length-1)]
        if segment not in segment_positions
            segment_positions[segment] = []
        segment_positions[segment].APPEND(pos)
    
    // 计算每个片段的出现次数
    for segment, positions in segment_positions
        segment_counts[segment] = LENGTH(positions)
    
    // 步骤2: 使用种子扩展技术查找潜在重复片段 - O(|Q|)
    potential_repeats = []
    
    for i = 0 to LENGTH(query) - min_length
        // 正向检查
        query_segment = query[i..(i+min_length-1)]
        if query_segment in segment_counts and segment_counts[query_segment] = 1
            ref_pos = segment_positions[query_segment][0]
            // 向两边扩展找到最大匹配
            extended_length = EXTEND-MATCH(query, reference, i, ref_pos, min_length)
            if extended_length ≥ min_length
                potential_repeats.APPEND({
                    position: ref_pos,
                    length: extended_length,
                    query_position: i
                })
        
        // 反向互补检查
        rev_comp = GET-REVERSE-COMPLEMENT(query_segment)
        if rev_comp in segment_counts and segment_counts[rev_comp] = 1
            ref_pos = segment_positions[rev_comp][0]
            extended_length = EXTEND-MATCH-REVERSE-COMPLEMENT(query, reference, i, ref_pos, min_length)
            if extended_length ≥ min_length
                potential_repeats.APPEND({
                    position: ref_pos,
                    length: extended_length,
                    query_position: i,
                    is_reverse: TRUE
                })
    
    // 步骤3: 检查连续重复 - O(|Q|)
    repeats = []
    
    for repeat in potential_repeats
        position = repeat.position
        length = repeat.length
        query_position = repeat.query_position
        is_reverse = repeat.is_reverse or FALSE
        
        if is_reverse
            segment = GET-REVERSE-COMPLEMENT(reference[position..(position+length-1)])
        else
            segment = reference[position..(position+length-1)]
        
        // 检查连续重复
        count = 1
        current_pos = query_position + length
        
        while current_pos + length ≤ LENGTH(query)
            if query[current_pos..(current_pos+length-1)] = segment
                count = count + 1
                current_pos = current_pos + length
            else
                break
        
        if count > 1
            repeats.APPEND({
                position: position,
                length: length,
                repeat_count: count,
                is_reverse: is_reverse,
                original_sequence: reference[position..(position+length-1)],
                query_position: query_position
            })
    
    // 步骤4: 过滤和排序 - O(N log N)
    filtered_repeats = FILTER-NESTED-REPEATS(repeats)
    SORT(filtered_repeats) by (length * repeat_count) in descending order
    
    return filtered_repeats

EXTEND-MATCH(query, reference, query_pos, ref_pos, seed_length)
    // 初始长度为种子长度
    length = seed_length
    
    // 向右扩展
    q_right = query_pos + seed_length
    r_right = ref_pos + seed_length
    
    while q_right < LENGTH(query) and r_right < LENGTH(reference) and query[q_right] = reference[r_right]
        length = length + 1
        q_right = q_right + 1
        r_right = r_right + 1
    
    // 向左扩展
    q_left = query_pos - 1
    r_left = ref_pos - 1
    
    while q_left ≥ 0 and r_left ≥ 0 and query[q_left] = reference[r_left]
        length = length + 1
        q_left = q_left - 1
        r_left = r_left - 1
    
    return length
```

##### 优化算法复杂度分析

**时间复杂度**：O(|R| × |Q|)
- 步骤1 - 为参考序列创建索引：O(|R|)
- 步骤2 - 使用种子扩展查找潜在重复：O(|Q|)
- 步骤3 - 检查连续重复：O(|Q|)
- 步骤4 - 过滤和排序：O(N log N)，其中N是找到的重复片段数量（通常N << |R|和|Q|）

总体上，优化算法的时间复杂度从O(|R|² × |Q|)降低到O(|R| × |Q|)，这在处理长序列时能显著提高性能。

**空间复杂度**：O(|R| + |Q| + N)
- 索引数据结构：O(|R|)
- 潜在重复列表：O(|Q|)，最坏情况下每个查询位置都有一个潜在重复
- 结果存储：O(N)

#### 优化创新点

1. **索引预处理**：构建参考序列的索引，使得只需O(1)时间即可检查片段的唯一性，避免了重复计算。

2. **种子扩展法**：从最小长度的匹配开始，通过双向扩展寻找最长匹配，避免了枚举所有可能长度的开销。

3. **哈希优化**：使用哈希表存储片段位置和出现次数，使查询操作的时间复杂度保持在O(1)级别。

4. **减少冗余检查**：只对参考序列中唯一出现的片段执行重复检测，大大减少了不必要的计算。

#### 算法实现对比

两种算法在实现上的主要区别：

1. **基础算法**：通过双重循环枚举参考序列中所有可能的片段，然后检查每个片段是否在查询序列中连续重复。

2. **优化算法**：先为参考序列建立索引，然后扫描查询序列，利用索引快速定位可能的匹配，再通过扩展种子找到最长匹配，最后检查连续重复。

优化算法特别适合处理长DNA序列，当序列长度增加时，性能优势更加显著。在处理几十万或几百万碱基的序列时，优化算法可以将运行时间从数小时缩短到几分钟。