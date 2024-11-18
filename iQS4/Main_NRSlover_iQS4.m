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
t = 0.1;      %板厚
EELM = Enod(:,1);   %输入测量应变的单元
Eep = EptranPLN(t);    %转换为8应变分量Eep(单元数1,8)
Eep2 = EptranPLN(t);
% -------------------------------------------
NSTP = 200;
%--------------------------------------------
[Ndsp, iterations] = NR_Solver_iQS4(Nxy,Enod,EELM,t,Eep,Eep2,1e-6,100,NSTP);