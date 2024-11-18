clear;
%输入矩阵Nxy,Enod:
%Enod:[单元号|节点1|节点2|节点3|节点4]
%Nxy:[节点号|x坐标|y坐标|z坐标]
%--------------------------------------------
Enod = importdata('Enod.txt');
Enod = [(1:size(Enod,1))',Enod];
%--------------------------------------------
Nxy = importdata('Nxy0.txt');
Nxy = [(1:size(Nxy,1))',Nxy];
%--------------------------------------------
t = 30;      %板厚
EELM = Enod(:,1);   %输入测量应变的单元
Eep = EptranPLN_ss('NETP.txt');    %转换为8应变分量Eep(单元数1,8)
Eep2 = zeros(size(Eep));
%-----------------------------------------------
[GK,B] = gstiffm_ssiQS4(Nxy,Enod,EELM,t);
GF = gforcem_ssiQS4(Nxy, Enod, Eep,Eep2,EELM,t);      %等效载荷向量
%----------------------------------------------
FIX = importdata('fix.txt');
[GK,GF] = boundary(FIX,[1,2,3,4,5,6],GK,GF);
%----------------------------------------------
Ndsp = SloveiQS4n(GK,GF);   %返回节点位移值与约束反力*