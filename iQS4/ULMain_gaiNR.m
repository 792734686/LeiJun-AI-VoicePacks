clear
%输入矩阵Nxy,Enod:
%Enod:[单元号|节点1|节点2|节点3|节点4]
%Nxy:[节点号|x坐标|y坐标|z坐标]
%--------------------------------------------
Enod = importdata('Enod.rpt').data;
Enod = [(1:size(Enod,1))',Enod];
%--------------------------------------------
Nxy = importdata('Nxy.rpt').data;
Nxy0 = Nxy(:,1:4);
Nxy = Nxy0;
Rxy = zeros(size(Nxy,1),3); %节点转角矩阵，RX|RY|RZ
FIX = importdata('fix.rpt');%边界条件矩阵，节点号|方向|位移量

EELM = Enod(:,1);   %输入测量应变的单元
%-----------------------------------------------
t = 20;      %板厚
NSTP = 1;   %应变数据的分析步个数
%-----------------------------------------------
tic;
for i = 1:NSTP
    Eep = NSTP\EptranPLN(t);    %转换为8应变分量Eep(单元数1,8)
    Eep2 = NSTP\EptranPLN(t);
    %----------------------------------------------

    [Ndsp, iterations] = NR_Solver_iQS4(Nxy,Enod,EELM,t,Eep,Eep2,1e-6,20,i); %返回节点位移值与约束反力*
    Nxy = [(1:size(Nxy,1))',Nxy(:,2:4)+Ndsp(:,2:4)];
    Rxy = Rxy+Ndsp(:,5:7);
end
toc;
Ndsp = [(1:size(Nxy,1))',Nxy(:,2:4)-Nxy0(:,2:4),Rxy];