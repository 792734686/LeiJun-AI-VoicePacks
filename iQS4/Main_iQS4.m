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
t = 1;      %板厚
% EELM = importdata('EELM1.rpt');
EELM = Enod(:,1);   %输入测量应变的单元
Eep = EptranPLN(t);    %转换为8应变分量Eep(单元数1,8)
Eep2 = zeros(size(Eep));
%-----------------------------------------------
tic;
[GK,B] = gstiffm_iQS4(Nxy,Enod,EELM,t);
GF = gforcem_iQS4(Nxy, Enod, Eep,Eep2,EELM,t);      %等效载荷向量
%----------------------------------------------
FIX = importdata('fix.rpt');
[GK,GF] = boundary(FIX,[1,2,3,4,5,6],GK,GF);
%----------------------------------------------
Ndsp = SloveiQS4n(GK,GF);   %返回节点位移值与约束反力*
toc;