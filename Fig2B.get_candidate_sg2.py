import itertools

# 原始序列
seq = "GTTCGAGAGCTATGCTGGAAACAGCATAGCAAGTTCGAA"

# 定义配对位置 (1-based 到 0-based 转换)
pairs = [(2,39),(3,38),(4,37),(5,36),(6,35),
         (9,30),(10,29),(11,28),(12,27),
         (13,26),(14,25),(15,24),(16,23),(17,22)]
pairs = [(i-1, j-1) for i,j in pairs]

# 合法碱基对
valid_pairs = ["AT","TA","GC","CG"]

def mutate_pair(a, b):
    """返回某一对碱基的所有新配对（排除原样）"""
    original = a+b
    return [p for p in valid_pairs if p != original]

def gc_content(s):
    gc = s.count("G") + s.count("C")
    return gc / len(s)

def has_four_consecutive(s):
    for i in range(len(s)-3):
        if s[i] == s[i+1] == s[i+2] == s[i+3]:
            return True
    return False

def generate_mutants(seq, pairs, k):
    """
    在 k 对碱基上突变，返回所有符合要求的突变序列
    """
    from copy import deepcopy
    seq_list = list(seq)
    
    for chosen_pairs in itertools.combinations(range(len(pairs)), k):
        # 每对的可能突变集合
        mutation_options = []
        for idx in chosen_pairs:
            i, j = pairs[idx]
            a, b = seq_list[i], seq_list[j]
            mutation_options.append(mutate_pair(a,b))
        
        # 笛卡尔积：组合所有可能突变
        for combo in itertools.product(*mutation_options):
            new_seq = seq_list[:]
            for idx, new_pair in zip(chosen_pairs, combo):
                i, j = pairs[idx]
                new_seq[i], new_seq[j] = new_pair[0], new_pair[1]
            s = "".join(new_seq)
            if not has_four_consecutive(s) and 0.3 <= gc_content(s) <= 0.7:
                yield s

if __name__ == "__main__":
    for k in range(1, 15):
        filename = f"MM{k}.txt"
        with open(filename, "w") as f:
            for mutant in generate_mutants(seq, pairs, k):
                f.write(mutant + "\n")
        print(f"突变 {k} 对完成，结果写入 {filename}")

