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
FIX = importdata('fix.rpt');
XSYMM = importdata('xsymm.rpt');
YSYMM = importdata('ysymm.rpt');
ZSYMM = importdata('zsymm.rpt');
EELM = importdata('EELM.rpt');   %输入测量应变的单元
%-----------------------------------------------
t = 20;      %板厚
NSTP = 10;   %应变数据的分析步个数
%-----------------------------------------------
tic;
for i = 1:NSTP
    Eep = NSTP\EptranPLN(t);    %转换为8应变分量Eep(单元数1,8)
    Eep2 = NSTP\Eptran2PLN(t);
    [GK,B] = gstiffm_iQS4(Nxy,Enod,EELM,t);
    GF = gforcem_iQS4(Nxy, Enod,Eep,Eep2,EELM,t);      %等效载荷向量
    %----------------------------------------------
    [GK,GF] = boundary(FIX,[1,2,3,4,5,6],GK,GF);
    [GK,GF] = boundary(XSYMM,[1,5,6],GK,GF);
    [GK,GF] = boundary(YSYMM,[2,4,6],GK,GF);
    [GK,GF] = boundary(ZSYMM,[3,4,5],GK,GF);
    Ndsp = SloveiQS4n(GK,GF);   %返回节点位移值与约束反力*
    Nxy = [(1:size(Nxy,1))',Nxy(:,2:4)+Ndsp(:,2:4)];
    Rxy = Rxy+Ndsp(:,5:7);
end
toc;
Ndsp = [(1:size(Nxy,1))',Nxy(:,2:4)-Nxy0(:,2:4),Rxy];