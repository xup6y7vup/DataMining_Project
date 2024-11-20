# 計算圖的數量
def count_graphs(file_path):
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            # 檢查每行是否以 "t #" 開頭
            if line.startswith("t #"):
                count += 1
    return count

# 設定 UK_graph.txt 的路徑
file_path = "C:/Users/Lockdream/Desktop/data_mining/Pang/data/UK/UK_graph.txt"

# 計算並顯示圖的總數
graph_count = count_graphs(file_path)
print("圖的總數量為:", graph_count)


# 計算 patterns 的數量
def count_patterns(file_path):
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            # 檢查每行是否以 "t #" 開頭
            if line.startswith("t #"):
                count += 1
    return count

# 設定 UK_pattern.txt 的路徑
pattern_file_path = "C:/Users/Lockdream/Desktop/data_mining/Pang/data/UK/UK_pattern.txt"

# 計算並顯示 patterns 的總數量
pattern_count = count_patterns(pattern_file_path)
print("patterns 的總數量為:", pattern_count)
