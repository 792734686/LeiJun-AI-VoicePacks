#coding=utf-8
import numpy as np# 可输出txt文件
from odbAccess import*# 打开odb的库
import os
import sys

myodb=openOdb(path='Job-1.odb')

sets = ['SET-8']
boundtype = ['fix']

for i in range(len(sets)):
    setname = sets[i]
    boundname = boundtype[i]+'.txt'
    nodelabel = []
    for setnode in myodb.rootAssembly.nodeSets[setname].nodes[0]:
        nodelabel.append(setnode.label)

        output_file_path = boundname
        # 使用with语句来打开文件，保证文件正常关闭
        with open(output_file_path, 'w') as file:
            # 遍历列表中的每一个元素
            for node in nodelabel:
                # 将列表中的每个元素转化为字符串
                # 使用str.join将数值转换为以逗号分隔的字符串
                line = str(node) + '\n'
                # 写入文件
                file.write(line)
