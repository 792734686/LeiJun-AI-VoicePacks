#coding=utf-8

# 读取txt文件中以制表符分隔的每一行数据的第二个数据并按照数字大小排列，输出到一个txt文件中
with open('NBD.txt', 'r') as file:
    lines = file.readlines()
    data = [line.split()[1] for line in lines if len(line) > 1]

sorted_data = sorted(data, key=lambda x: float(x))  # 按照数字大小排序数据

with open('NBD.rpt', 'w') as file:
    for item in sorted_data:
        file.write("%s\n" % item)  # 将排序后的数据写入新的txt文件中
        

# 读取txt文件中以制表符分隔的每一行数据的第二个数据并按照数字大小排列，输出到一个txt文件中
with open('EELM.txt', 'r') as file:
    lines = file.readlines()
    data = [line.split()[1] for line in lines if len(line) > 1]

sorted_data = sorted(data, key=lambda x: float(x))  # 按照数字大小排序数据

with open('EELM.rpt', 'w') as file:
    for item in sorted_data:
        file.write("%s\n" % item)  # 将排序后的数据写入新的txt文件中
