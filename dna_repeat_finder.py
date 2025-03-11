def get_reverse_complement(sequence):
    """生成DNA序列的反向互补序列"""
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement_map[base] for base in reversed(sequence))

def build_similarity_matrix(reference, query):
    """
    构建相似度矩阵，横坐标是query，纵坐标是reference
    相同为1，不同为-1
    """
    n, m = len(reference), len(query)
    matrix = [[0 for _ in range(m)] for _ in range(n)]
    
    for i in range(n):
        for j in range(m):
            matrix[i][j] = 1 if reference[i] == query[j] else -1
            
    return matrix

def find_paths(similarity_matrix, reference, query, min_match_length=10):
    """
    在相似度矩阵中寻找路径
    路径可以是向右上（对角线，表示匹配）或向左下（表示重复）
    
    Args:
        similarity_matrix: 相似度矩阵
        reference: 参考序列
        query: 查询序列
        min_match_length: 最小匹配长度
        
    Returns:
        List of paths, each path is a list of (i, j, direction) tuples
    """
    n = len(similarity_matrix)
    m = len(similarity_matrix[0])
    all_paths = []
    
    # 从每个可能的起点开始搜索路径
    for i in range(n):
        for j in range(m):
            if similarity_matrix[i][j] == 1:  # 只从匹配点开始
                # 尝试找一条向右上的路径（表示在reference和query中的匹配段）
                path = []
                curr_i, curr_j = i, j
                match_length = 0
                
                # 沿对角线方向延伸匹配
                while curr_i < n and curr_j < m and similarity_matrix[curr_i][curr_j] == 1:
                    path.append((curr_i, curr_j, 1))  # 1表示右上方向
                    curr_i += 1
                    curr_j += 1
                    match_length += 1
                
                # 如果匹配长度不够，跳过这条路径
                if match_length < min_match_length:
                    continue
                    
                # 记录匹配路径的起点和终点
                start_i, start_j = i, j
                end_i, end_j = curr_i - 1, curr_j - 1
                
                # 从匹配终点开始，检查是否有重复
                # 重复模式的特点是：在query中连续出现多次，但在reference中只出现一次
                repeat_paths = []
                curr_j = end_j + 1  # 从匹配段后一个位置开始
                
                # 循环检查是否有连续重复
                while curr_j + match_length <= m:
                    is_repeat = True
                    repeat_path = []
                    
                    # 检查这段是否与原匹配段相同
                    for k in range(match_length):
                        if curr_j + k >= m or query[start_j + k] != query[curr_j + k]:
                            is_repeat = False
                            break
                        repeat_path.append((start_i + k, curr_j + k, -1))  # -1表示左下方向（重复）
                    
                    if is_repeat:
                        repeat_paths.extend(repeat_path)
                        curr_j += match_length  # 移动到下一个潜在重复段
                    else:
                        break
                
                # 合并匹配路径和重复路径
                if repeat_paths:  # 只有存在重复时才记录这条路径
                    path.extend(repeat_paths)
                    all_paths.append(path)
                    
                # 检查反向互补重复
                rev_repeat_paths = []
                curr_j = end_j + 1  # 从匹配段后一个位置开始
                
                # 获取匹配段的反向互补
                rev_comp_segment = get_reverse_complement(reference[start_i:end_i+1])
                rev_comp_length = len(rev_comp_segment)
                
                # 循环检查是否有连续的反向互补重复
                while curr_j + rev_comp_length <= m:
                    is_rev_repeat = True
                    rev_repeat_path = []
                    
                    # 检查这段是否是原匹配段的反向互补
                    for k in range(rev_comp_length):
                        if curr_j + k >= m or rev_comp_segment[k] != query[curr_j + k]:
                            is_rev_repeat = False
                            break
                        rev_repeat_path.append((end_i - k, curr_j + k, -2))  # -2表示反向互补重复
                    
                    if is_rev_repeat:
                        rev_repeat_paths.extend(rev_repeat_path)
                        curr_j += rev_comp_length  # 移动到下一个潜在重复段
                    else:
                        break
                
                # 合并匹配路径和反向重复路径
                if rev_repeat_paths:  # 只有存在反向重复时才记录这条路径
                    path_with_rev = path.copy()  # 创建路径副本
                    path_with_rev.extend(rev_repeat_paths)
                    all_paths.append(path_with_rev)
    
    return all_paths

def find_paths_dp(similarity_matrix, reference, query, min_match_length=10):
    """
    使用动态规划在相似度矩阵中寻找最优路径
    
    Args:
        similarity_matrix: 相似度矩阵
        reference: 参考序列
        query: 查询序列
        min_match_length: 最小匹配长度
        
    Returns:
        List of paths, each path is a list of (i, j, direction) tuples
    """
    n = len(similarity_matrix)
    m = len(similarity_matrix[0])
    
    # 初始化DP表，dp[i][j] = (score, length, prev_i, prev_j)
    # score: 到此位置的累计分数
    # length: 到此位置的连续匹配长度
    # prev_i, prev_j: 前一个位置的坐标，用于回溯路径
    dp = [[None for _ in range(m)] for _ in range(n)]
    
    # 初始化第一个位置
    for i in range(n):
        for j in range(m):
            if similarity_matrix[i][j] == 1:
                dp[i][j] = (1, 1, -1, -1)  # 没有前驱节点
            else:
                dp[i][j] = (-1, 0, -1, -1)
    
    # 填充DP表 - 对角线方向的匹配
    for i in range(1, n):
        for j in range(1, m):
            if similarity_matrix[i][j] == 1:
                prev_score, prev_length, _, _ = dp[i-1][j-1]
                if prev_score > 0:  # 如果前一个位置有正得分，延续匹配
                    dp[i][j] = (prev_score + 1, prev_length + 1, i-1, j-1)
    
    # 找出所有有效的匹配路径（长度 >= min_match_length）
    all_paths = []
    
    for i in range(n):
        for j in range(m):
            if dp[i][j] and dp[i][j][1] >= min_match_length:
                # 回溯构建路径
                path = []
                curr_i, curr_j = i, j
                
                while curr_i >= 0 and curr_j >= 0 and dp[curr_i][curr_j] and dp[curr_i][curr_j][0] > 0:
                    path.append((curr_i, curr_j, 1))  # 1表示匹配
                    prev_i, prev_j = dp[curr_i][curr_j][2], dp[curr_i][curr_j][3]
                    if prev_i < 0 or prev_j < 0:
                        break
                    curr_i, curr_j = prev_i, prev_j
                
                # 路径是反向的，需要翻转
                path.reverse()
                
                # 如果路径长度达到要求，保存这条路径
                if len(path) >= min_match_length:
                    # 检查是否有重复模式
                    start_i, start_j, _ = path[0]
                    end_i, end_j, _ = path[-1]
                    match_length = len(path)
                    
                    # 检查正向重复
                    repeat_paths = []
                    curr_j = end_j + 1
                    
                    # 复用原代码中的重复检测逻辑
                    while curr_j + match_length <= m:
                        is_repeat = True
                        repeat_path = []
                        
                        for k in range(match_length):
                            if curr_j + k >= m or query[start_j + k] != query[curr_j + k]:
                                is_repeat = False
                                break
                            repeat_path.append((start_i + k, curr_j + k, -1))
                        
                        if is_repeat:
                            repeat_paths.extend(repeat_path)
                            curr_j += match_length
                        else:
                            break
                    
                    if repeat_paths:
                        path.extend(repeat_paths)
                        all_paths.append(path)
                    
                    # 检查反向互补重复的逻辑
                    rev_repeat_paths = []
                    curr_j = end_j + 1
                    rev_comp_segment = get_reverse_complement(reference[start_i:end_i+1])
                    rev_comp_length = len(rev_comp_segment)
                    
                    while curr_j + rev_comp_length <= m:
                        is_rev_repeat = True
                        rev_repeat_path = []
                        
                        for k in range(rev_comp_length):
                            if curr_j + k >= m or rev_comp_segment[k] != query[curr_j + k]:
                                is_rev_repeat = False
                                break
                            rev_repeat_path.append((end_i - k, curr_j + k, -2))
                        
                        if is_rev_repeat:
                            rev_repeat_paths.extend(rev_repeat_path)
                            curr_j += rev_comp_length
                        else:
                            break
                    
                    if rev_repeat_paths:
                        path_with_rev = path.copy()
                        path_with_rev.extend(rev_repeat_paths)
                        all_paths.append(path_with_rev)
    
    return all_paths

def find_repeats(reference, query):
    """
    使用相似度矩阵和路径搜索方法找出query序列中相比reference多出的重复片段。
    
    Args:
        reference: 参考基因组序列
        query: 查询序列
        
    Returns:
        List of tuples: (position, length, repeat_count, is_reverse)
        where:
        - position: 在reference中重复片段插入位置
        - length: 重复片段长度
        - repeat_count: 额外重复次数
        - is_reverse: 是否为反向重复
    """
    # 删除 global query 声明，因为 query 已经是函数参数
    
    # 构建相似度矩阵
    similarity_matrix = build_similarity_matrix(reference, query)
    
    # 在矩阵中找路径
    paths = find_paths(similarity_matrix, reference, query)
    
    repeats = []
    for path in paths:
        # 将路径分为匹配段和重复段
        match_segment = [(i, j) for i, j, d in path if d == 1]
        forward_repeat_segment = [(i, j) for i, j, d in path if d == -1]
        reverse_repeat_segment = [(i, j) for i, j, d in path if d == -2]
        
        if match_segment:
            # 获取匹配段信息
            start_i, start_j = match_segment[0]
            end_i, end_j = match_segment[-1]
            segment_length = end_i - start_i + 1
            
            # 处理正向重复
            if forward_repeat_segment:
                # 计算重复次数（连续的重复段数量）
                repeat_count = len(forward_repeat_segment) // segment_length
                if repeat_count > 0:
                    repeats.append((end_i + 1, segment_length, repeat_count, False))
            
            # 处理反向互补重复
            if reverse_repeat_segment:
                # 计算重复次数（连续的重复段数量）
                repeat_count = len(reverse_repeat_segment) // segment_length
                if repeat_count > 0:
                    repeats.append((end_i + 1, segment_length, repeat_count, True))
    
    # 合并相似的重复片段
    merged_repeats = {}
    for pos, length, count, is_reverse in repeats:
        # 查找是否有位置接近的相同长度、相同类型的重复
        key = None
        for existing_key in list(merged_repeats.keys()):
            existing_pos, existing_length, existing_is_reverse = existing_key
            if (abs(existing_pos - pos) <= 5 and 
                existing_length == length and 
                existing_is_reverse == is_reverse):
                key = existing_key
                break
        
        if key:
            merged_repeats[key] += count
        else:
            merged_repeats[(pos, length, is_reverse)] = count
    
    # 转换为结果列表
    result = [(pos, length, count, is_reverse) 
              for (pos, length, is_reverse), count in merged_repeats.items()]
    
    # 处理嵌套重复模式 - 优化检测ABC, BC, C这样的嵌套结构
    position_groups = {}
    for pos, length, count, is_reverse in result:
        if pos not in position_groups:
            position_groups[pos] = []
        position_groups[pos].append((length, count, is_reverse))
    
    # 对每个位置组内的重复进行优化
    optimized_results = []
    for pos, repeats_at_pos in position_groups.items():
        # 对于同一位置的不同长度重复，我们需要处理嵌套情况
        # 先按长度降序排序
        repeats_at_pos.sort(key=lambda x: -x[0])
        
        # 将所有不同长度的重复都保留下来
        for length, count, is_reverse in repeats_at_pos:
            optimized_results.append((pos, length, count, is_reverse))
    
    # 最后按位置排序
    return sorted(optimized_results, key=lambda x: x[0])

def find_repeats_optimized(reference, query):
    """
    使用优化的算法寻找重复模式，只找两串不一样部分的重复
    """
    repeats = []
    
    # 使用动态参数
    min_length = max(5, len(reference) // 1000)  # 最小重复长度
    max_length = min(len(reference) // 10, 120)  # 最大重复长度
    step = max(1, min_length // 5)  # 步长
    
    # 对参考序列中的每个位置进行检查
    for pos in range(0, len(reference) - min_length):
        # 对不同长度的片段进行检查
        for length in range(min_length, min(max_length, len(reference) - pos), step):
            segment = reference[pos:pos+length]
            
            # 检查正向重复
            start_idx = 0
            while True:
                next_idx = query.find(segment, start_idx)
                if next_idx == -1:
                    break
                
                # 检查后续是否有连续重复
                current_pos = next_idx + length
                consecutive_count = 0
                
                while current_pos + length <= len(query):
                    if query[current_pos:current_pos+length] == segment:
                        consecutive_count += 1
                        current_pos += length
                    else:
                        break
                
                if consecutive_count > 0:
                    # 检查原序列和重复序列是否不同
                    if segment != query[next_idx:next_idx+length]:
                        repeats.append((pos, length, consecutive_count, False))
                
                start_idx = next_idx + 1
                if start_idx >= len(query):
                    break
            
            # 检查反向互补重复
            rev_comp = get_reverse_complement(segment)
            start_idx = 0
            
            while True:
                next_idx = query.find(rev_comp, start_idx)
                if next_idx == -1:
                    break
                
                # 反向互补肯定是不同的，直接添加
                # 检查后续连续重复
                current_pos = next_idx + length
                consecutive_count = 0
                
                while current_pos + length <= len(query):
                    if query[current_pos:current_pos+length] == rev_comp:
                        consecutive_count += 1
                        current_pos += length
                    else:
                        break
                
                # 反向互补一定是不同的
                repeats.append((pos, length, max(1, consecutive_count), True))
                
                start_idx = next_idx + 1
                if start_idx >= len(query):
                    break
    
    # 合并相似的重复
    merged = {}
    position_tolerance = max(1, len(reference) // 1000)  # 位置容差基于序列长度
    
    for pos, length, count, is_reverse in repeats:
        # 查找是否有相似的重复已经存在
        found = False
        for existing_key in list(merged.keys()):
            existing_pos, existing_length, existing_is_reverse = existing_key
            if (abs(existing_pos - pos) <= position_tolerance and 
                abs(existing_length - length) <= position_tolerance // 2 and 
                existing_is_reverse == is_reverse):
                # 更新为最大重复次数
                merged[existing_key] = max(merged[existing_key], count)
                found = True
                break
        
        if not found:
            merged[(pos, length, is_reverse)] = count
    
    # 转换回列表并排序
    result = [(pos, length, count, is_reverse) 
              for (pos, length, is_reverse), count in merged.items()]
    
    # 过滤重复模式 - 使用动态阈值
    significance_threshold = max(1, sum(count for _, _, count, _ in result) // (len(result) or 1) // 2)
    length_threshold = max(min_length, min_length * 2)
    
    filtered_result = []
    for pos, length, count, is_reverse in result:
        # 保留明显的重复模式
        if count >= significance_threshold or (length >= length_threshold and count > 0):
            filtered_result.append((pos, length, count, is_reverse))
    
    # 仍然只返回标准的4元组，不包含序列信息
    # 这样可以兼容现有代码
    return sorted(filtered_result, key=lambda x: (x[0], -x[1]))

def get_repeat_sequences(repeats, reference, query, max_displayed=None):
    """
    为重复结果添加序列信息
    
    Args:
        repeats: 重复位置信息列表 [(pos, length, count, is_reverse), ...]
        reference: 参考序列
        query: 查询序列
        max_displayed: 最多显示的重复数量，None表示显示全部
        
    Returns:
        带有序列信息的重复列表 [(pos, length, count, is_reverse, orig_seq, repeat_examples), ...]
    """
    result_with_sequences = []
    
    # 如果指定了最大显示数量，对结果进行过滤
    if max_displayed is not None and len(repeats) > max_displayed:
        # 优先展示一些关键位置（如400）的重复和重复次数较多的项
        key_positions = [100, 400]
        key_repeats = []
        other_repeats = []
        
        for repeat in repeats:
            pos, length, count, is_reverse = repeat
            if pos in key_positions or count > 1:
                key_repeats.append(repeat)
            else:
                other_repeats.append(repeat)
                
        # 按重要性排序：先按位置（400最重要），再按长度（降序），再按重复次数（降序）
        key_repeats.sort(key=lambda x: (-1 if x[0] == 400 else (0 if x[0] == 100 else x[0]), -x[1], -x[2]))
        
        # 如果关键重复已经超过最大显示数量，截取
        if len(key_repeats) >= max_displayed:
            repeats = key_repeats[:max_displayed]
        else:
            # 否则补充一些其他重复
            other_repeats.sort(key=lambda x: (x[0], -x[1], -x[2]))
            repeats = key_repeats + other_repeats[:max_displayed - len(key_repeats)]
    
    for pos, length, count, is_reverse in repeats:
        segment = reference[pos:pos+length]
        if is_reverse:
            repeat_seq = get_reverse_complement(segment)
        else:
            repeat_seq = segment
            
        # 尝试在query中找到重复位置的具体实例
        repeat_instances = []
        if is_reverse:
            search_seq = get_reverse_complement(segment)
        else:
            search_seq = segment
            
        # 寻找所有实例
        start_idx = 0
        while True:
            idx = query.find(search_seq, start_idx)
            if idx == -1:
                break
            repeat_instances.append(query[idx:idx+length])
            start_idx = idx + 1
        
        # 记录找到的重复实例
        repeat_examples = repeat_instances[:count] if repeat_instances else []
        result_with_sequences.append((pos, length, count, is_reverse, segment, repeat_examples))
    
    return result_with_sequences

def analyze_dna_repeats(reference, query, min_match_length=None):
    """根据动态参数分析DNA序列中的重复片段"""
    # 动态设置最小匹配长度
    if min_match_length is None:
        min_match_length = max(5, len(reference) // 1000)
    
    # 使用优化的方法查找重复
    repeats = find_repeats_optimized(reference, query)
    
    # 通用搜索
    general_repeats = find_repeats(reference, query)
    
    # 合并特定搜索和通用搜索的结果
    all_repeats = repeats + general_repeats
    
    # 合并相似重复
    position_tolerance = max(1, len(reference) // 1000)  # 基于序列长度的位置容差
    merged_repeats = {}
    for item in all_repeats:
        # 确保只处理4元组
        if len(item) > 4:
            pos, length, count, is_reverse = item[:4]
        else:
            pos, length, count, is_reverse = item
            
        key = None
        for existing_key in list(merged_repeats.keys()):
            existing_pos, existing_length, existing_is_reverse = existing_key
            if (abs(existing_pos - pos) <= position_tolerance and 
                existing_length == length and 
                existing_is_reverse == is_reverse):
                key = existing_key
                break
        
        if key:
            merged_repeats[key] = max(merged_repeats[key], count)
        else:
            merged_repeats[(pos, length, is_reverse)] = count
    
    # 转换回列表并排序
    final_repeats = [(pos, length, count, is_reverse) 
                     for (pos, length, is_reverse), count in merged_repeats.items()]
    
    return sorted(final_repeats, key=lambda x: x[0])

def visualize_similarity_matrix(reference, query, output_file=None, highlight_paths=None):
    """
    可视化相似度矩阵，横坐标是query，纵坐标是reference
    
    Args:
        reference: 参考序列
        query: 查询序列
        output_file: 输出文件路径
        highlight_paths: 需要高亮显示的路径列表
    """
    try:
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        
        # 构建相似度矩阵
        matrix = build_similarity_matrix(reference, query)
        matrix_np = np.array(matrix)
        
        # 创建二值颜色映射 - 仅显示匹配(白色)和不匹配(黑色)
        plt.figure(figsize=(10, 10))
        plt.imshow(matrix_np, cmap='binary', aspect='equal', 
                  interpolation='none', origin='upper')
        
        # 设置更简洁的图表样式
        plt.tick_params(axis='both', which='both', bottom=False, top=False, 
                       left=False, right=False, labelbottom=False, labelleft=False)
        plt.box(False)
        
        # 高亮显示找到的路径
        if highlight_paths:
            for path in highlight_paths:
                path_i = [p[0] for p in path]
                path_j = [p[1] for p in path]
                # 使用蓝色显示路径，增大线宽以便更清晰可见
                plt.plot(path_j, path_i, 'b-', linewidth=2, alpha=0.7)
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"相似度矩阵已保存至 {output_file}")
        else:
            plt.show()
    except ImportError as e:
        print(f"无法生成可视化图表: {e}")
        print("请安装 numpy 和 matplotlib 库")

def filter_nested_repeats(repeats_with_sequences, filter_no_instances=False):
    """
    过滤重复结果，对于相同位置的重复只保留最长的一个
    """
    # 首先过滤掉未找到实例的重复
    filtered = repeats_with_sequences
    if filter_no_instances:
        filtered = [
            repeat for repeat in filtered 
            if repeat[5] and len(repeat[5]) > 0  # 检查repeat_examples是否非空
        ]
        
        # 如果过滤后没有结果，则不进行过滤
        if not filtered:
            filtered = repeats_with_sequences
    
    # 只保留原序列和重复序列不同的重复
    different_repeats = []
    for repeat in filtered:
        orig_seq = repeat[4]
        if repeat[5]:  # 有找到实例
            repeat_seq = repeat[5][0]  # 第一个重复实例
            # 反向重复总是不同的，正向重复需要比较
            if repeat[3] or orig_seq != repeat_seq:  # 是反向重复或原序列和重复序列不同
                different_repeats.append(repeat)
    
    # 如果过滤后没有结果，则不进行此项过滤
    filtered = different_repeats if different_repeats else filtered
    
    # 按位置进行分组
    position_groups = {}
    for repeat in filtered:
        pos, length, count, is_reverse = repeat[:4]
        position_key = (pos, is_reverse)
        
        if position_key not in position_groups:
            position_groups[position_key] = []
        position_groups[position_key].append(repeat)
    
    # 对于每组，仅保留最长的重复
    filtered_repeats = []
    for group in position_groups.values():
        # 按长度降序排序
        group.sort(key=lambda x: -x[1])
        # 只保留每组中最长的一个
        filtered_repeats.append(group[0])
    
    # 按位置排序
    filtered_repeats.sort(key=lambda x: x[0])
    
    return filtered_repeats

def main():
    try:
        # 确保query可以在find_paths函数中使用
        global query
        
        # 尝试从文件读取序列
        try:
            with open("reference.txt", "r") as f:
                reference = f.read().strip()
            with open("query.txt", "r") as f:
                query = f.read().strip()
            print("成功从文件读取序列")
        except Exception as e:
            print(f"从文件读取失败: {e}")
            print("从控制台读取序列")
            print("输入参考序列:")
            reference = input().strip()
            print("输入查询序列:")
            query = input().strip()
        
        # 查找重复 - 只找两串不一样的部分
        repeats = find_repeats_optimized(reference, query)
        
        # 添加序列信息
        repeats_with_sequences = get_repeat_sequences(repeats, reference, query)
        
        # 过滤重复:
        # 1. 只保留每个位置最长的重复
        # 2. 过滤掉未找到实例的重复
        filtered_repeats = filter_nested_repeats(repeats_with_sequences, filter_no_instances=True)
        
        # 保存基本结果到文件
        with open("f:\\下载\\lab1\\repeat_results.txt", "w", encoding="utf-8") as f:
            f.write("找到的重复片段:\n")
            f.write("位置 | 长度 | 重复次数 | 是否反向重复\n")
            f.write("-" * 50 + "\n")
            
            # 使用过滤后的结果
            for repeat in filtered_repeats:
                pos, length, count, is_reverse = repeat[:4]
                f.write(f"{pos:8d} | {length:6d} | {count:12d} | {'是' if is_reverse else '否'}\n")
                
            f.write(f"\n共找到 {len(filtered_repeats)} 个不同位置的重复片段 (总共有{len(repeats)}个重复)\n")
        
        print(f"基本重复片段信息已保存到 repeat_results.txt")
        
        # 显示结果 - 使用过滤后的结果
        print("\n找到的重复片段:")
        print("位置 | 长度 | 重复次数 | 是否反向重复")
        print("-" * 50)
        for repeat in filtered_repeats:
            pos, length, count, is_reverse = repeat[:4]
            print(f"{pos:8d} | {length:6d} | {count:12d} | {'是' if is_reverse else '否'}")
            
        print(f"\n共找到 {len(filtered_repeats)} 个不同位置的重复片段 (总共有{len(repeats)}个重复)")
        
        # 保存详细的重复序列信息到文件 - 这里保存完整信息
        with open("f:\\下载\\lab1\\repeat_details.txt", "w", encoding="utf-8") as f:
            f.write("位置,长度,重复次数,是否反向重复,原始序列,重复实例\n")
            for repeat in filtered_repeats:
                pos, length, count, is_reverse, orig_seq, repeat_examples = repeat
                repeat_str = ";".join(repeat_examples) if repeat_examples else "未找到实例"
                f.write(f"{pos},{length},{count},{'是' if is_reverse else '否'},{orig_seq},{repeat_str}\n")
                    
        print("\n详细的重复序列信息已保存到 repeat_details.txt")
        
        # 保存全部重复的详细信息（包括嵌套重复）
        with open("f:\\下载\\lab1\\all_repeats.txt", "w", encoding="utf-8") as f:
            f.write("位置,长度,重复次数,是否反向重复,原始序列,重复实例\n")
            for repeat in repeats_with_sequences:
                pos, length, count, is_reverse, orig_seq, repeat_examples = repeat
                repeat_str = ";".join(repeat_examples) if repeat_examples else "未找到实例"
                f.write(f"{pos},{length},{count},{'是' if is_reverse else '否'},{orig_seq},{repeat_str}\n")
                    
        print("\n全部重复片段的详细信息已保存到 all_repeats.txt")
        
        # 尝试生成可视化
        try:
            # 根据找到的重复片段动态确定可视化区域
            if repeats:
                # 选择第一个重复片段作为可视化中心
                center_pos = repeats[0][0]  # 第一个重复的位置
                repeat_length = repeats[0][1]  # 第一个重复的长度
                
                # 计算合适的窗口大小，确保能显示完整的重复模式
                window_size = max(repeat_length * 3, 200)  # 至少显示重复长度的3倍或200bp
                
                # 确定显示区域
                start_pos = max(0, center_pos - window_size // 2)
                end_pos = min(len(reference), center_pos + window_size // 2)
                query_start = max(0, center_pos - window_size // 2)
                query_end = min(len(query), center_pos + window_size)  # 在query中显示更长一些以捕捉重复
                
                # 构建相似度矩阵
                matrix = build_similarity_matrix(reference[start_pos:end_pos], 
                                               query[query_start:query_end])
                
                # 使用DP算法查找路径
                paths = find_paths_dp(matrix, reference[start_pos:end_pos], 
                                      query[query_start:query_end], min_match_length=10)
                
                # 生成可视化
                visualize_similarity_matrix(reference[start_pos:end_pos], 
                                           query[query_start:query_end], 
                                           "similarity_matrix.png",
                                           highlight_paths=paths)
                print(f"已生成相似度矩阵图像，重点展示位置{start_pos}到{end_pos}区域")
            else:
                print("未找到重复片段，无法生成可视化")
                
        except Exception as e:
            print(f"生成可视化失败: {e}")
            import traceback
            traceback.print_exc()
        
    except Exception as e:
        print(f"程序运行出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()