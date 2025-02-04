function GF = gforcem_ssiQS4(NXY, ELEM,Eep,Eep2,EELM,t)
% function to assemble global stiffness matrix 
% NXY --- node coordinate matrix
% ELEM --- element node number matrix
% EMAT --- elemnt material matrix
% GK --- global stiffness matrix
% Eep --- 总应变矩阵，行号为单元号，膜应变（3列）|弯曲应变（3列）|法向剪切应变（2列）
n=size(NXY,1);      %返回结构总节点个数
GF = zeros(6*n,1);
n=size(ELEM,1);     % 返回结构总单元个数
for i=1:n
%     E = EMAT(i,2);
%     v = EMAT(i,3);
%     t = EMAT(i,4);
%     D = [1,v,0;v,1,0;0,0,(1-v)/2];
%     D = E/(1-v*v)*D;
%     Eepi = Eep;
    if ismember(i,EELM)
        we = 1;
        wg = 1e-4;
        Eepi = Eep(i,:)';       %返回单元载荷向量
    else
        we = 1;
        wg = 1e-4;
        Eepi = Eep2(i,:)'; 
    end

    NN = ELEM(i,2:end);     %返回节点编号
    ENC = [NXY(NN,2),NXY(NN,3),NXY(NN,4)];    %4*2，返回单元的节点坐标（行），4列表示4节点
    EF = eforcem_ssiQS4( ENC,4,Eepi,we,wg,t );  % 4--ngs
    n1 = NN(1);
    n2 = NN(2);
    n3 = NN(3);
    n4 = NN(4);
    F = zeros(6*size(NXY,1),1);
    F(6*n1-5:6*n1,1) = EF(1:6);
    F(6*n2-5:6*n2,1) = EF(7:12);
    F(6*n3-5:6*n3,1) = EF(13:18);
    F(6*n4-5:6*n4,1) = EF(19:24);
    GF = GF+F;
end
end

