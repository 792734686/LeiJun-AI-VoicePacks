clear
%输入矩阵Nxy,Enod:
%Enod:[单元号|节点1|节点2|节点3|节点4]
%Nxy:[节点号|x坐标|y坐标|z坐标]
% ---------------输出文件目录------------------
OTPfilename = 'file3_Ndsp';
if exist ([cd,'\',OTPfilename,'\'],'dir') == 7
    rmdir(OTPfilename,'s')
end

mkdir(OTPfilename)
%--------------------------------------------
Enod = importdata('Enod.txt');
Enod = [(1:size(Enod,1))',Enod];
%--------------------------------------------
Nxy = importdata('Nxy0.txt');
Nxy0 = [(1:size(Nxy,1))',Nxy];
Nxy = Nxy0;
NDSP = zeros(size(Nxy,1),3);
Rxy = zeros(size(Nxy,1),3); %节点转角矩阵，RX|RY|RZ
FIX = importdata('fix.txt');
% BOUND = importdata('bound.txt');
% XSYMM = importdata('xsymm.rpt');
% YSYMM = importdata('ysymm.rpt');
% ZSYMM = importdata('zsymm.rpt');

% EELM = importdata('EELM.rpt');   %输入测量应变的单元
EELM = Enod(:,1);
%-----------------------------------------------
t = 2;      %板厚
strainTfile = [cd,'\LE\'];
totalframe = length(dir(strcat(strainTfile,'\','*.txt')))/2-1;   %应变数据的分析步个数

%-----------------------------------------------
tic;
for framenum = 1:totalframe
        Eep = EptranPLN_MTFrames(strainTfile,framenum,t) - ...
            EptranPLN_MTFrames(strainTfile,framenum - 1,t);    %转换为8应变分量Eep(单元数1,8)
        Eep2 = zeros(size(Eep));
%     % ------------------------------------------------
%     NEBT_file = ['NEBT_,',num2str(i),'.txt'];
%     NETP_file = ['NETP_,',num2str(i),'.txt'];
%     if i == 1
%         Eep = EptranPLN_fileInput(['NEBT_',num2str(1),'.txt'],['NETP_',num2str(1),'.txt'],t);
%         Eep2 = zeros(size(Eep));
%     else
%         Eep = EptranPLN_fileInput(['NEBT_',num2str(i),'.txt'],['NETP_',num2str(i),'.txt'],t)-...
%             EptranPLN_fileInput(['NEBT_',num2str(i-1),'.txt'],['NETP_',num2str(i-1),'.txt'],t);
%         Eep2 = zeros(size(Eep));
%     end
    % ------------------------------------------------
    
    GK = gstiffm_iTR3(Nxy,Enod,EELM,t);
    GF = gforcem_iTR3(Nxy, Enod,Eep,Eep2,EELM,t);      %等效载荷向量
    %----------------------------------------------
    [GK,GF] = boundary(FIX,[1,2,3,4,5,6],GK,GF);
%     [GK,GF] = boundary(BOUND,[1,2,4,5],GK,GF);
%     [GK,GF] = boundary(XSYMM,[1,5,6],GK,GF);
%     [GK,GF] = boundary(YSYMM,[2,4,6],GK,GF);
%     [GK,GF] = boundary(ZSYMM,[3,4,5],GK,GF);
    Ndsp = SloveiQS4n(GK,GF);   %返回节点位移值与约束反力*
    NDSP = NDSP + Ndsp(:,2:4);
    Rxy = RotMat_node(Rxy',Ndsp(:,5:7)');
    Rxy = Rxy';
    Nxy = [(1:size(Nxy,1))',Nxy(:,2:4)+Ndsp(:,2:4)];
    filename = [cd,'\',OTPfilename,'\','Ndsp_',num2str(framenum),'.txt'];
    writematrix(NDSP,filename,'delimiter','\t');
    filename = [cd,'\',OTPfilename,'\','NU_',num2str(framenum),'.txt'];
    writematrix([NDSP,Rxy],filename,'delimiter','\t');
    disp(['Frame_', num2str(framenum),' Solve Succeed!']);
end
toc;
Ndsp = [(1:size(Nxy,1))',Nxy(:,2:4)-Nxy0(:,2:4),Rxy];