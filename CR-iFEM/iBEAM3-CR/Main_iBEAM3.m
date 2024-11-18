clear;
%输入矩阵Nxy,Enod:
%Enod:[单元号|节点1|节点2]
%Nxy:[节点号|x坐标|y坐标|z坐标]
%--------------------------------------------
Enod = importdata('Enod.txt');
Enod = [(1:size(Enod,1))',Enod];
%--------------------------------------------
Nxy = importdata('Nxy0.txt');
Nxy = [(1:size(Nxy,1))',Nxy];
%--------------------------------------------
% EELM = importdata('EELM.txt');
EELM = Enod(:,1);   %输入测量应变的单元
Eep = beamstrain([1,1,1,1,1,1]);
%-----------------------------------------------
tic;
[GK,GF] = gstiff_force_iBEAM3(Enod,Nxy,EELM,Eep);
%----------------------------------------------
FIX = importdata('fix.txt');
[GK,GF] = boundary(FIX,[1,2,3,4,5,6],GK,GF);
%----------------------------------------------
Ndsp = Slove_6DOF(GK,GF);   %返回节点位移值与约束反力*
toc;

function SE = beamstrain(switcher)
full_strain = importdata('SE.txt');
SE = zeros(size(full_strain));
for i = 1:size(switcher,2)
    SE(:,i) = switcher(i)*full_strain(:,i);
end
end