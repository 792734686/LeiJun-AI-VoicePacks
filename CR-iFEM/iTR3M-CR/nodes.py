#coding=utf-8
from abaqus import *
from abaqusConstants import *
import __main__
import regionToolset

modelname = 'Model-2'
instancename = 'Part-1-1'
# 获取当前模型、装配和步骤
myModel = mdb.models[modelname]
myAssembly = myModel.rootAssembly

Elements = myAssembly.instances[instancename].elements
nodes = []

for elm in Elements:
    elementlable = elm.label
    parentnode = elm.connectivity
    nodes.append([parentnode[0],parentnode[1],parentnode[2]])

output_file_path = 'Enod.txt'
# 使用with语句来打开文件，保证文件正常关闭
with open(output_file_path, 'w') as file:
    # 遍历列表中的每一个元素
    for element_nodes in nodes:
        # 将列表中的每个元素转化为字符串
        # 使用str.join将数值转换为以逗号分隔的字符串
        line = '  '.join(map(str, element_nodes)) + '\n'
        # 写入文件
        file.write(line)