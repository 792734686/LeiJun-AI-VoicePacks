#coding=utf-8
import numpy as np# 可输出txt文件
from odbAccess import*# 打开odb的库
import os
import sys

myodb=openOdb(path='Job-3.odb')
instancename = 'YUWEI-COPY-1'
Stepname = 'Step-1'
SETName = 'SET-3'
myinstance = myodb.rootAssembly.instances[instancename]
origin_coordinate = []
totalFramenum = len(myodb.steps[Stepname].frames)
region1 = myodb.rootAssembly.instances[instancename].nodeSets[SETName]

workpath = os.getcwd()
NUfile = "NU"
path = workpath + "./" + NUfile
if not os.path.exists(path):
    os.makedirs(path)
NUpath = workpath + ".\\" + NUfile + "\\"

for node in myinstance.nodes:
    node_origin_coordinate = node.coordinates
    origin_coordinate.append(node_origin_coordinate)

with open('Nxy0.txt', 'w') as file:
    # 遍历列表中的每一个元素
    for coordinates in origin_coordinate:
        # 将列表中的每个元素转化为字符串
        # 使用str.join将数值转换为以逗号分隔的字符串
        line = '  '.join(map(str, coordinates)) + '\n'
        # 写入文件
        file.write(line)

for framenum in range(totalFramenum):

    displacement = myodb.steps[Stepname].frames[framenum].fieldOutputs['U'].getSubset(region = region1).values
    NU = []

    for node_disp_value in displacement:
        nodedisp = node_disp_value.data
        NU.append(nodedisp)

    with open(NUpath + 'NU_{}.txt'.format(framenum), 'w') as file:
        # 遍历列表中的每一个元素
        for nu in NU:
            # 将列表中的每个元素转化为字符串
            # 使用str.join将数值转换为以逗号分隔的字符串
            line = '  '.join(map(str, nu)) + '\n'
            # 写入文件
            file.write(line)

    def plus(a,b):
        c = []
        for i in range(len(a)):
            c1 = a[i][0]+b[i][0]
            c2 = a[i][1]+b[i][1]
            c3 = a[i][2]+b[i][2]
            c.append([c1,c2,c3])
        return c

    deformed_coordinate = plus(origin_coordinate,NU)
    with open(NUpath + 'Nxy_{}.txt'.format(framenum), 'w') as file:
        # 遍历列表中的每一个元素
        for dcoordinates in deformed_coordinate:
            # 将列表中的每个元素转化为字符串
            # 使用str.join将数值转换为以逗号分隔的字符串
            line = '  '.join(map(str, dcoordinates)) + '\n'
            # 写入文件
            file.write(line)

