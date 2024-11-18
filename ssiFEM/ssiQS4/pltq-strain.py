#coding=utf-8
#import pandas as pd
from odbAccess import*# 打开odb的库
data_list = []
data_list1=[]
#单元对应的节点号
with open('node.txt') as file:
    lines = file.readlines()
    # 分割每行，并转换为整数
    data_list = [line.strip().split() for line in lines]  # split() 用于分割行内容

# 将字符串列表转换为整数列表
int_data_list = []
for sublist in data_list:
    int_sublist = []
    for item in sublist:
        try:
            number = int(item)
            int_sublist.append(number)
        except ValueError:
            print("1")
    int_data_list.append(int_sublist)
#每个节点的坐标
with open('nodezuobiao1.txt', 'r') as file:
    lines = file.readlines()

# lines 现在包含了文件的每一行，格式未改变。

# 如果要将数据存储到一个数组或列表中，可以这样做：
xyzdata = []  # 创建一个空列表用于存储数据
for line in lines:
    # 假定每行的数据通过空格分隔
    row = line.strip().split()  # 去掉开头和结尾的空白符，然后分割
    xyzdata.append(row)  # 将分割后的数据作为列表添加到data列表中
float_data = [[float(val) for val in row] for row in xyzdata]
#
odb=openOdb(path='Job-hk.odb')# 读取odb
scratchOdb = session.ScratchOdb(odb)
a=len(int_data_list)
for item in range(a):
    #print(data_list[0][0])
    node1=int_data_list[item][0]
    node2 = int_data_list[item][1]
    node3 = int_data_list[item][3]
    scratchOdb.rootAssembly.DatumCsysByThreePoints(name='CSYS-%s' %item,
                                                   coordSysType=CARTESIAN,
                                                   origin=(float_data[node1-1][0], float_data[node1-1][1],
                                                           float_data[node1-1][2]),
                                                   point1=(float_data[node2-1][0], float_data[node2-1][1],
                                                           float_data[node2-1][2]),
                                                   point2=(float_data[node3-1][0], float_data[node3-1][1],
                                                           float_data[node3-1][2]))

step = odb.steps['Step-1']
strain_bottom=[]
strain_top=[]
# 得到最后一帧的结果数据，也可以是其他帧
lastFrame = step.frames[1]
# 获取应变场
strain_field = lastFrame.fieldOutputs['E']
for item in range(a):
    myCsys = scratchOdb.rootAssembly.datumCsyses['CSYS-%s' %item]
    StrainTransformed = strain_field.getTransformedField(datumCsys=myCsys)
    selected_strains = []  # 对于E11, E22, E12
    #print(odb.rootAssembly.instances['PART-1-1'])
    region = odb.rootAssembly.instances['PART-1-1'].elementSets['SET-1']
    # 遍历所有实例和相应的单元
        # 获取该实例中单元集合的应变数据
    strains = StrainTransformed.getSubset(region=region)
    # 遍历每个单元的应变数据


    for value in strains.values:
        # 提取特定的应变分量 (E11, E22, E12)
        e11 = value.data[0]  # E11
        e22 = value.data[1]  # E22
        e12 = value.data[3]  # E12

        # 创建一个包含单元标号和特定应变分量的列表
        element_strains = [value.elementLabel, e11, e22, e12]

        # 将该列表添加到 selected_strains 中
        selected_strains.append(element_strains)
    strain_bottom.append(selected_strains[item])
    strain_top.append(selected_strains[item+a])
#输出底部应变
output_file_path = 'bottom-strains.txt'
# 使用with语句来打开文件，保证文件正常关闭
with open(output_file_path, 'w') as file:
    # 遍历列表中的每一个元素
    for element_strain in strain_bottom:
        # 将列表中的每个元素转化为字符串
        # 使用str.join将数值转换为以逗号分隔的字符串
        line = '  '.join(map(str, element_strain)) + '\n'
        # 写入文件
        file.write(line)
#输出顶部应变
output_file_path = 'top-strains.txt'
# 使用with语句来打开文件，保证文件正常关闭
with open(output_file_path, 'w') as file:
    # 遍历列表中的每一个元素
    for element_strain in strain_top:
        # 将列表中的每个元素转化为字符串
        # 使用str.join将数值转换为以逗号分隔的字符串
        line = '  '.join(map(str, element_strain)) + '\n'
        # 写入文件
        file.write(line)
#output_file_path = 'all-strains2.txt'
# 使用with语句来打开文件，保证文件正常关闭
#with open(output_file_path, 'w') as file:
    # 遍历列表中的每一个元素
   # for element_strain in selected_strains:
        # 将列表中的每个元素转化为字符串
        # 使用str.join将数值转换为以逗号分隔的字符串
      #  line = '  '.join(map(str, element_strain)) + '\n'
        # 写入文件
       # file.write(line)