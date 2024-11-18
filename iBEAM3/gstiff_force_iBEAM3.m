function [GK,GF] = gstiff_force_iBEAM3(Enod,NXY,EELM,Eep)
n=size(NXY,1);
GK = zeros(6*n,6*n);
GF = zeros(6*n,1);
n = size(Enod,1);
for i=1:n
    %     E = EMAT(i,2);
    %     v = EMAT(i,3);
    %     t = EMAT(i,4);
    %     D = [1,v,0;v,1,0;0,0,(1-v)/2];
    %     D = E/(1-v*v)*D;
    if ismember(i,EELM)
        w = [1,1,1,1e-4,1e-4,1];
        e = Eep(i,:);
    else
        w = 1e-4*[1,1,1,1,1,1];
        e = zeros(1,6);
    end
    NN = Enod(i,2:end); %第i个单元的节点号向量，节点1|2
    ENC = [NXY(NN,2),NXY(NN,3),NXY(NN,4)]';    %返回第i个单元的2个节点坐标，[X1|X2]
% ---------------------------------------------------
% 静力凝聚
    EK = estiff_iBEAM3(w,ENC,2);  % 2--ngs
    EF = eforce_iBEAM3(w,e,ENC,2);
    KRR = EK(1:12,1:12);
    KRO = EK(1:12,13:14);
    KOO = EK(13:14,13:14);
    FR = EF(1:12);
    FO = EF(13:14);
    EK = KRR-KRO*KOO^(-1)*KRO';
    EF = FR-KRO*KOO^(-1)*FO;
% 旋转
    Te = transmat_iBEAM3(ENC);
    EK = Te'*EK*Te;
    EF = Te'*EF;
% ----------------------------------------------------
    for j=1:2   % 4 - number of node in an element
        for k=1:2   % 4 - number of node in an element
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
    n1 = NN(1);
    n2 = NN(2);
    F = zeros(6*size(NXY,1),1);
    F(6*n1-5:6*n1,1) = EF(1:6);
    F(6*n2-5:6*n2,1) = EF(7:12);
    GF = GF + F;
end
end