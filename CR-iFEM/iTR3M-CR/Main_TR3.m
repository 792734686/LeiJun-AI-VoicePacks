clear;
%输入矩阵Nxy,Enod:
%Enod:[单元号|节点1|节点2|节点3|节点4]
%Nxy:[节点号|x坐标|y坐标|z坐标]
%--------------------------------------------
Enod = importdata('Enod.rpt').data;
Enod = [(1:size(Enod,1))',Enod];
%--------------------------------------------
Nxy = importdata('Nxy.rpt').data;
Nxy = Nxy(:,1:4);
%--------------------------------------------
t = 0.04;      %板厚
% EELM = importdata('EELM1.rpt');
EELM = Enod(:,1);   %输入测量应变的单元
Eep = zeros(size(Enod,1),1);    %转换为8应变分量Eep(单元数1,8)
Eep2 = zeros(size(Eep));
%-----------------------------------------------
tic;
GK = gstiffm_TR3(Nxy,Enod,EELM,t);
GF = zeros(size(Nxy,1)*6,1);
GF(6*(1-1)+3) = 20;
GF(6*(601-1)+1) = -20;
%载荷向量
%----------------------------------------------
XSY = importdata('xsymm.rpt');
[GK,GF] = boundary(XSY,[1,5,6],GK,GF);
ZSY = importdata('zsymm.rpt');
[GK,GF] = boundary(ZSY,[3,4,5],GK,GF);
U2 = importdata('u2.rpt');
[GK,GF] = boundary(U2,2,GK,GF);
%----------------------------------------------
Ndsp = SloveiQS4n(GK,GF);   %返回节点位移值与约束反力*
toc;