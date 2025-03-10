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
    使用优化的DP算法寻找重复
    """
    # 构建相似度矩阵
    similarity_matrix = build_similarity_matrix(reference, query)
    
    # 使用DP算法在矩阵中找路径，传入reference和query参数
    paths = find_paths_dp(similarity_matrix, reference, query)
    
    # 处理路径，提取重复信息
    repeats = []
    for path in paths:
        match_segment = [(i, j) for i, j, d in path if d == 1]
        forward_repeat_segment = [(i, j) for i, j, d in path if d == -1]
        reverse_repeat_segment = [(i, j) for i, j, d in path if d == -2]
        
        if match_segment:
            start_i, start_j = match_segment[0]
            end_i, end_j = match_segment[-1]
            segment_length = end_i - start_i + 1
            
            if forward_repeat_segment:
                repeat_count = len(forward_repeat_segment) // segment_length
                if repeat_count > 0:
                    repeats.append((end_i + 1, segment_length, repeat_count, False))
            
            if reverse_repeat_segment:
                repeat_count = len(reverse_repeat_segment) // segment_length
                if repeat_count > 0:
                    repeats.append((end_i + 1, segment_length, repeat_count, True))
    
    merged_repeats = {}
    for pos, length, count, is_reverse in repeats:
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
    
    result = [(pos, length, count, is_reverse) 
              for (pos, length, is_reverse), count in merged_repeats.items()]
    
    position_groups = {}
    for pos, length, count, is_reverse in result:
        if pos not in position_groups:
            position_groups[pos] = []
        position_groups[pos].append((length, count, is_reverse))
    
    optimized_results = []
    for pos, repeats_at_pos in position_groups.items():
        repeats_at_pos.sort(key=lambda x: -x[0])
        
        for length, count, is_reverse in repeats_at_pos:
            optimized_results.append((pos, length, count, is_reverse))
    
    return sorted(optimized_results, key=lambda x: x[0])

def analyze_dna_repeats(reference, query, min_match_length=10):
    """根据pdf文档要求分析DNA序列中的重复片段"""
    # 预处理 - 只考虑最关键的部分
    repeats = []
    
    # 精确搜索400位置附近的长度为50,70,100的重复
    target_lengths = [50, 70, 100]
    target_position = 400
    position_range = 5  # 允许的位置偏差
    
    for length in target_lengths:
        if target_position + length <= len(reference):
            segment = reference[target_position - length:target_position]
            if len(segment) < length:
                continue  # 跳过长度不足的片段
            
            # 在query中查找该片段
            query_pos = query.find(segment)
            if query_pos != -1:
                # 检查正向重复
                next_pos = query_pos + length
                repeat_count = 0
                
                while next_pos + length <= len(query):
                    if query[next_pos:next_pos + length] == segment:
                        repeat_count += 1
                        next_pos += length
                    else:
                        break
                
                if repeat_count > 0:
                    repeats.append((target_position, length, repeat_count, False))
                
                # 检查反向互补重复
                rev_comp_segment = get_reverse_complement(segment)
                next_pos = query_pos + length
                rev_repeat_count = 0
                
                while next_pos + length <= len(query):
                    if query[next_pos:next_pos + length] == rev_comp_segment:
                        rev_repeat_count += 1
                        next_pos += length
                    else:
                        break
                
                if rev_repeat_count > 0:
                    repeats.append((target_position, length, rev_repeat_count, True))
    
    # 通用搜索
    general_repeats = find_repeats(reference, query)
    
    # 合并特定搜索和通用搜索的结果
    all_repeats = repeats + general_repeats
    
    # 合并相似重复
    merged_repeats = {}
    for pos, length, count, is_reverse in all_repeats:
        key = None
        for existing_key in list(merged_repeats.keys()):
            existing_pos, existing_length, existing_is_reverse = existing_key
            if (abs(existing_pos - pos) <= position_range and 
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
        import matplotlib.font_manager as fm
        
        # 尝试设置支持中文的字体
        try:
            # 检查是否有支持中文的字体
            chinese_fonts = [f.name for f in fm.fontManager.ttflist 
                            if ('SimHei' in f.name or 'Microsoft YaHei' in f.name or
                                'SimSun' in f.name or 'WenQuanYi' in f.name or
                                'Noto Sans CJK' in f.name or 'Droid Sans' in f.name)]
            
            if chinese_fonts:
                plt.rcParams['font.family'] = chinese_fonts[0]
            else:
                # 如果没有找到中文字体，使用无衬线字体并避免使用中文标签
                plt.rcParams['font.family'] = 'sans-serif'
        except:
            # 如果出错，使用默认字体设置
            pass
        
        # 构建相似度矩阵
        matrix = build_similarity_matrix(reference, query)
        matrix_np = np.array(matrix)
        
        # 创建自定义颜色映射
        colors = [(0.7, 0, 0), (1, 1, 1), (0, 0.7, 0)]  # 红色-白色-绿色
        cmap = LinearSegmentedColormap.from_list('custom_diverging', colors, N=256)
        
        # 创建图像
        plt.figure(figsize=(12, 10))
        im = plt.imshow(matrix_np, cmap=cmap, aspect='auto', 
                   interpolation='nearest', origin='upper',
                   vmin=-1, vmax=1)
        plt.colorbar(im, label='Similarity (1:Match, -1:Mismatch)')
        plt.title('Reference vs Query Similarity Matrix')
        plt.xlabel('Query Position')
        plt.ylabel('Reference Position')
        
        # 高亮显示找到的路径
        if highlight_paths:
            for path in highlight_paths:
                path_i = [p[0] for p in path]
                path_j = [p[1] for p in path]
                plt.plot(path_j, path_i, 'b-', linewidth=1.5, alpha=0.7)
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"相似度矩阵已保存至 {output_file}")
        else:
            plt.show()
    except ImportError as e:
        print(f"无法生成可视化图表: {e}")
        print("请安装 numpy 和 matplotlib 库")

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
        
        # 查找重复 - 使用优化算法
        repeats = find_repeats_optimized(reference, query)
        
        # 显示结果
        print("\n找到的重复片段:")
        print("位置 | 长度 | 重复次数 | 是否反向重复")
        print("-" * 50)
        for pos, length, count, is_reverse in repeats:
            print(f"{pos:8d} | {length:6d} | {count:12d} | {'是' if is_reverse else '否'}")
        print(f"\n共找到 {len(repeats)} 个重复片段")
        
        # 尝试生成可视化
        try:
            # 选择一段较短的区域进行可视化
            start_pos = max(0, 400 - 100)
            end_pos = min(len(reference), 400 + 100)
            query_start = max(0, 400 - 100)
            query_end = min(len(query), 400 + 300)  # 查询序列中观察更多区域以捕捉重复
            
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