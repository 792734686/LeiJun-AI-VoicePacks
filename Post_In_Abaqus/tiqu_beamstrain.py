#coding=utf-8
import numpy as np# 可输出txt文件
from odbAccess import*# 打开odb的库
import os
import sys
myodb=openOdb(path='Job-2.odb')
instancename = 'PART-1-1'
framenum = -1
Stepname = 'Step-1'
myinstance = myodb.rootAssembly.instances[instancename]

ELM_SE = myodb.steps[Stepname].frames[framenum].fieldOutputs['SE'].values
beamstrain1 = []
for SE in ELM_SE:
    section_strain = SE.data
    beamstrain1.append(section_strain)

ELM_SK = myodb.steps[Stepname].frames[framenum].fieldOutputs['SK'].values
beamstrain2 = []
for SK in ELM_SK:
    section_kappa = SK.data
    beamstrain2.append(section_kappa)

beamstrain = []
for i in range(len(beamstrain1)):
    e1 = beamstrain1[i][0]
    e2 = beamstrain2[i][0]
    e3 = - beamstrain2[i][1]
    e4 = beamstrain1[i][1]
    e5 = beamstrain1[i][2]
    e6 = beamstrain2[i][2]
    beamstrain.append([e1,e2,e3,e4,e5,e6])

with open('SE.txt','w') as file:
    # 遍历列表中的每一个元素
    for E in beamstrain:
        # 将列表中的每个元素转化为字符串
        # 使用str.join将数值转换为以逗号分隔的字符串
        line = '  '.join(map(str,E)) + '\n'
        # 写入文件
        file.write(line)