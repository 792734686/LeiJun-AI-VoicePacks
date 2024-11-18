#coding=utf-8
from abaqus import *
from abaqusConstants import *
import __main__
import regionToolset

# 获取当前模型、装配和步骤
myModel = mdb.models['Model-2']
myAssembly = myModel.rootAssembly
myStep = myModel.StaticStep(name='Step-1', previous='Initial')

# 要施加位移边界条件的节点集合
#nodes = myAssembly.instances['Part-1-mesh-1-1'].nodes
nodes = mdb.models['Model-2'].rootAssembly.instances['Part-1-mesh-1-1'].nodes

# 打开文件
with open('Ndsp.txt', 'r') as file:
    # 逐行读取数据并存储为列表
    Ndsp = []
    for line in file:
        # 使用map函数将每行数据转换为整数或浮点数，并存储为列表
        Ndsp.append(list(map(float, line.split('\t'))))

# 循环遍历每个节点并添加位移边界条件
for node in nodes:
    nodelabel = node.label
    nodedisp = Ndsp[nodelabel-1]
    U1 = nodedisp[1]
    U2 = nodedisp[2]
    U3 = nodedisp[3]
    nodei = nodes[nodelabel-1:nodelabel]
    region = regionToolset.Region(nodes=nodei)
    myModel.DisplacementBC(name='BC-'+str(node.label), createStepName='Step-1', region=region, u1=U1, u2=U2, u3=U3)
