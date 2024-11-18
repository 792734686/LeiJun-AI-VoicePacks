#coding=utf-8
import numpy as np
from odbAccess import*# 打开odb的库
import os
import sys


workpath = os.getcwd()
strainfile = "LE"
path = workpath + "./" + strainfile
if not os.path.exists(path):
    os.makedirs(path)
strainpath = workpath + ".\\" + strainfile + "\\"

def frame(nodeCoordinateList):
    X1 = np.array(nodeCoordinateList[0])
    X2 = np.array(nodeCoordinateList[1])
    X3 = np.array(nodeCoordinateList[2])
    X12 = X2-X1
    X13 = X3-X1
    n1 = X12/np.linalg.norm(X12)
    n3 = np.cross(n1,X13)
    n3 = n3/np.linalg.norm(n3)
    n2 = np.cross(n3,n1)
    aFrame_mat = np.matrix([n1,n2,n3])
    return aFrame_mat

myodb=openOdb(path='Job-3.odb')
instancename = 'PART-1-1'
Stepname = 'Step-1'
SETName = 'SET-3'
myinstance = myodb.rootAssembly.instances[instancename]
totalFramenum = len(myodb.steps[Stepname].frames)
region1 = myodb.rootAssembly.instances[instancename].elementSets[SETName]
region2 = myodb.rootAssembly.instances[instancename].nodeSets[SETName]

initial_coordinate = []
for Nodes in myinstance.nodeSets[SETName].nodes:
    initial_coordinate.append(np.array(Nodes.coordinates))

with open ('Enod.txt','r') as file:
    Enod = []
    for line in file:
        Enod.append(list(line.strip('\n').split('  ')))
    for i in range(len(Enod)):
        Enod[i] = [int(x) for x in Enod[i]]

for framenum in range(totalFramenum):
    LEoutput = myodb.steps[Stepname].frames[framenum].fieldOutputs['LE'].getSubset(region = region1).values
    ELMnum = int(len(LEoutput) / 2)
    LEBT_init = []
    LETP_init = []
    frame_init = []
    for ELMS in range(ELMnum):
        ELM_lable = LEoutput[ELMS].elementLabel
        ELM_LEBT_init = [LEoutput[ELMS].data[0], LEoutput[ELMS].data[1], LEoutput[ELMS].data[3]]
        ELM_LEBT_init = np.array(ELM_LEBT_init)
        ELM_LETP_init = [LEoutput[ELMS+ELMnum].data[0], LEoutput[ELMS+ELMnum].data[1], LEoutput[ELMS+ELMnum].data[3]]
        ELM_LETP_init = np.array(ELM_LETP_init)
        ELM_frame_init = np.matrix(LEoutput[ELMS].localCoordSystem)
        LEBT_init.append(ELM_LEBT_init)
        LETP_init.append(ELM_LETP_init)
        frame_init.append(ELM_frame_init)

    Uoutput = myodb.steps[Stepname].frames[framenum].fieldOutputs['U'].getSubset(region = region2).values
    NU = []

    for Nodes in Uoutput:
        nodelable = Nodes.nodeLabel
        NU.append(np.array(Nodes.data))

    Nxy_mat = np.matrix(initial_coordinate)+np.matrix(NU)
    Nxy = np.matrix.tolist(Nxy_mat)

    frame_recent = []
    for ELMs in range(len(Enod)):
        node1 = int(Enod[ELMs][0])
        node2 = int(Enod[ELMs][1])
        node3 = int(Enod[ELMs][2])
        node1coordinate = Nxy[node1-1]
        node2coordinate = Nxy[node2-1]
        node3coordinate = Nxy[node3-1]
        frame_recent.append(frame([node1coordinate,node2coordinate,node3coordinate]))

    LEBT_recent = []
    LETP_recent = []
    for ELMs in range(len(Enod)):
        T0 = frame_init[ELMs]
        T = frame_recent[ELMs]
        R = T * T0.T
        ELMstrainBT = LEBT_init[ELMs].tolist()
        strain_tensor = np.matrix([[ELMstrainBT[0], 0.5*ELMstrainBT[2], 0], [0.5*ELMstrainBT[2], ELMstrainBT[1], 0], [0, 0, 0]])
        strain_tensor_mat = 2*R*strain_tensor*R.T
        for i in range(3):
            strain_tensor_mat[i, i] = 0.5 * strain_tensor_mat[i, i]
        strain_tensor = np.matrix.tolist(strain_tensor_mat)
        ELMstrainBT = [ELMs+1,strain_tensor[0][0], strain_tensor[1][1], strain_tensor[0][1]]
        LEBT_recent.append(ELMstrainBT)


        ELMstrainTP = LETP_init[ELMs].tolist()
        strain_tensor = np.matrix([[ELMstrainTP[0], 0.5*ELMstrainTP[2], 0], [0.5*ELMstrainTP[2], ELMstrainTP[1], 0], [0, 0, 0]])
        strain_tensor_mat = 2*R*strain_tensor*R.T
        for i in range(3):
            strain_tensor_mat[i, i] = 0.5 * strain_tensor_mat[i, i]
        strain_tensor = np.matrix.tolist(strain_tensor_mat)
        ELMstrainTP = [ELMs+1,strain_tensor[0][0], strain_tensor[1][1], strain_tensor[0][1]]
        LETP_recent.append(ELMstrainTP)

    with open(strainpath+'NEBT_{}.txt'.format(framenum), 'w') as file:
        # 遍历列表中的每一个元素
        for LEBT in LEBT_recent:
            # 将列表中的每个元素转化为字符串
            # 使用str.join将数值转换为以逗号分隔的字符串
            line = '  '.join(map(str, LEBT)) + '\n'
            # 写入文件
            file.write(line)

    with open(strainpath+'NETP_{}.txt'.format(framenum), 'w') as file:
        # 遍历列表中的每一个元素
        for LETP in LETP_recent:
            # 将列表中的每个元素转化为字符串
            # 使用str.join将数值转换为以逗号分隔的字符串
            line = '  '.join(map(str, LETP)) + '\n'
            # 写入文件
            file.write(line)


