function GK = gstiffm_iTR3(NXY, ELEM,EELM,t)
% function to assemble global stiffness matrix
% NXY --- node coordinate matrix    节点坐标矩阵，节点号|x坐标|y坐标
% ELEM --- element node number matrix   单元节点号矩阵，单元号|节点1|节点2|节点3|节点4
% EMAT --- elemnt material matrix
% GK --- global stiffness matrix    总刚度阵，6n*6n
% we --- 膜应变权重系数
% wk --- 弯曲应变权重系数
% wg --- 法向剪切应变权重系数
% t --- 板厚
% EELM --- 贴片的单元号
n=size(NXY,1);
GK = zeros(6*n,6*n);
n=size(ELEM,1);
for i=1:n
    %     E = EMAT(i,2);
    %     v = EMAT(i,3);
    %     t = EMAT(i,4);
    %     D = [1,v,0;v,1,0;0,0,(1-v)/2];
    %     D = E/(1-v*v)*D;
    if ismember(i,EELM)
        we = 1;
        wk = 1;
        wg = 1e-4;
    else
        we = 1e-4;
        wk = 1e-4;
        wg = 1e-4;
    end
    NN = ELEM(i,2:end); %第i个单元的节点号向量，节点1|2|3
    ENC = [NXY(NN,2),NXY(NN,3),NXY(NN,4)];    %返回第i个单元的3个节点坐标，[xi|yi|zi]
    EK = estiffm_iTR3( ENC,we,wk,wg,t );  % 4--ngs
    for j=1:3   % 3 - number of node in an element
        for k=1:3   % 3 - number of node in an element
            jj = NN(j);
            kk = NN(k);
            GK(6*jj-5,6*kk-5) = GK(6*jj-5,6*kk-5) + EK(6*j-5,6*k-5);
            GK(6*jj-5,6*kk-4) = GK(6*jj-5,6*kk-4) + EK(6*j-5,6*k-4);
            GK(6*jj-5,6*kk-3) = GK(6*jj-5,6*kk-3) + EK(6*j-5,6*k-3);
            GK(6*jj-5,6*kk-2) = GK(6*jj-5,6*kk-2) + EK(6*j-5,6*k-2);
            GK(6*jj-5,6*kk-1) = GK(6*jj-5,6*kk-1) + EK(6*j-5,6*k-1);
            GK(6*jj-5,6*kk) = GK(6*jj-5,6*kk) + EK(6*j-5,6*k);
            %-------------------------------------------------------
            GK(6*jj-4,6*kk-5) = GK(6*jj-4,6*kk-5) + EK(6*j-4,6*k-5);
            GK(6*jj-4,6*kk-4) = GK(6*jj-4,6*kk-4) + EK(6*j-4,6*k-4);
            GK(6*jj-4,6*kk-3) = GK(6*jj-4,6*kk-3) + EK(6*j-4,6*k-3);
            GK(6*jj-4,6*kk-2) = GK(6*jj-4,6*kk-2) + EK(6*j-4,6*k-2);
            GK(6*jj-4,6*kk-1) = GK(6*jj-4,6*kk-1) + EK(6*j-4,6*k-1);
            GK(6*jj-4,6*kk) = GK(6*jj-4,6*kk) + EK(6*j-4,6*k);
            %-------------------------------------------------------
            GK(6*jj-3,6*kk-5) = GK(6*jj-3,6*kk-5) + EK(6*j-3,6*k-5);
            GK(6*jj-3,6*kk-4) = GK(6*jj-3,6*kk-4) + EK(6*j-3,6*k-4);
            GK(6*jj-3,6*kk-3) = GK(6*jj-3,6*kk-3) + EK(6*j-3,6*k-3);
            GK(6*jj-3,6*kk-2) = GK(6*jj-3,6*kk-2) + EK(6*j-3,6*k-2);
            GK(6*jj-3,6*kk-1) = GK(6*jj-3,6*kk-1) + EK(6*j-3,6*k-1);
            GK(6*jj-3,6*kk) = GK(6*jj-3,6*kk) + EK(6*j-3,6*k);
            %-------------------------------------------------------
            GK(6*jj-2,6*kk-5) = GK(6*jj-2,6*kk-5) + EK(6*j-2,6*k-5);
            GK(6*jj-2,6*kk-4) = GK(6*jj-2,6*kk-4) + EK(6*j-2,6*k-4);
            GK(6*jj-2,6*kk-3) = GK(6*jj-2,6*kk-3) + EK(6*j-2,6*k-3);
            GK(6*jj-2,6*kk-2) = GK(6*jj-2,6*kk-2) + EK(6*j-2,6*k-2);
            GK(6*jj-2,6*kk-1) = GK(6*jj-2,6*kk-1) + EK(6*j-2,6*k-1);
            GK(6*jj-2,6*kk) = GK(6*jj-2,6*kk) + EK(6*j-2,6*k);
            %-------------------------------------------------------
            GK(6*jj-1,6*kk-5) = GK(6*jj-1,6*kk-5) + EK(6*j-1,6*k-5);
            GK(6*jj-1,6*kk-4) = GK(6*jj-1,6*kk-4) + EK(6*j-1,6*k-4);
            GK(6*jj-1,6*kk-3) = GK(6*jj-1,6*kk-3) + EK(6*j-1,6*k-3);
            GK(6*jj-1,6*kk-2) = GK(6*jj-1,6*kk-2) + EK(6*j-1,6*k-2);
            GK(6*jj-1,6*kk-1) = GK(6*jj-1,6*kk-1) + EK(6*j-1,6*k-1);
            GK(6*jj-1,6*kk) = GK(6*jj-1,6*kk) + EK(6*j-1,6*k);
            %-------------------------------------------------------
            GK(6*jj,6*kk-5) = GK(6*jj,6*kk-5) + EK(6*j,6*k-5);
            GK(6*jj,6*kk-4) = GK(6*jj,6*kk-4) + EK(6*j,6*k-4);
            GK(6*jj,6*kk-3) = GK(6*jj,6*kk-3) + EK(6*j,6*k-3);
            GK(6*jj,6*kk-2) = GK(6*jj,6*kk-2) + EK(6*j,6*k-2);
            GK(6*jj,6*kk-1) = GK(6*jj,6*kk-1) + EK(6*j,6*k-1);
            GK(6*jj,6*kk) = GK(6*jj,6*kk) + EK(6*j,6*k);
            %-------------------------------------------------------
        end
    end
end
end

