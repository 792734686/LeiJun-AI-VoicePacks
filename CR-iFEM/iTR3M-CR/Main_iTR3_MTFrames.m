clear;
%输入矩阵Nxy,Enod:
%Enod:[单元号|节点1|节点2|节点3|节点4]
%Nxy:[节点号|x坐标|y坐标|z坐标]
% ---------------输出文件目录------------------
OTPfilename = 'file4_Ndsp';
if exist ([cd,'\',OTPfilename,'\'],'dir') == 7
    rmdir(OTPfilename,'s')
end

mkdir(OTPfilename)
%--------------------------------------------
Enod = importdata('Enod.txt');
Enod = [(1:size(Enod,1))',Enod];
%--------------------------------------------
Nxy = importdata('Nxy0.txt');
Nxy = [(1:size(Nxy,1))',Nxy];
%--------------------------------------------
t = 2;      %板厚
strainTfile = [cd,'\LE\'];
totalframe = length(dir(strcat(strainTfile,'\','*.txt')))/2-1;   %应变数据的分析步个数
% -------------------------------------------
tic;
for framenum = 1:totalframe
% EELM = importdata('EELM.txt');
EELM = Enod(:,1);   %输入测量应变的单元
Eep = EptranPLN_MTFrames(strainTfile,framenum,t);    %转换为8应变分量Eep(单元数1,8)
Eep2 = zeros(size(Eep));
%-----------------------------------------------

GK = gstiffm_iTR3(Nxy,Enod,EELM,t);
GF = gforcem_iTR3(Nxy,Enod, Eep,Eep2,EELM,t);      %等效载荷向量
%----------------------------------------------
FIX = importdata('fix.txt');
% BOUND = importdata('bound.txt');
% XSYMM = importdata('xsymm.rpt');
% YSYMM = importdata('ysymm.rpt');
% ZSYMM = importdata('zsymm.rpt');
% U2 = importdata('u2.rpt');
[GK,GF] = boundary(FIX,[1,2,3,4,5,6],GK,GF);
% [GK,GF] = boundary(BOUND,[1,2,4,5],GK,GF);
% [GK,GF] = boundary(XSYMM,[1,5,6],GK,GF);
% [GK,GF] = boundary(YSYMM,[2,4,6],GK,GF);
% [GK,GF] = boundary(ZSYMM,[3,4,5],GK,GF);
% [GK,GF] = boundary(U2,2,GK,GF);
%----------------------------------------------
Ndsp = SloveiQS4n(GK,GF);   %返回节点位移值与约束反力*
filename = [cd,'\',OTPfilename,'\','Ndsp_',num2str(framenum),'.txt'];
writematrix(Ndsp(:,2:4),filename,'delimiter','\t');
filename = [cd,'\',OTPfilename,'\','NU_',num2str(framenum),'.txt'];
writematrix(Ndsp(:,2:7),filename,'delimiter','\t');
disp(['Frame_', num2str(framenum),' Solve Succeed!']);
end
toc;