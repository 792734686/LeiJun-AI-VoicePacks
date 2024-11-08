function [GK,GF] = cr_globalKF_iBEAM3(NXY,NXY0,Enod,E,E2,EELM,Ndsp,NOD_ROT)
n=size(NXY,1);
GK = zeros(6*n,6*n);
GF = zeros(6*n,1);
n=size(Enod,1);
NXY(:,2:4) = NXY(:,2:4)+Ndsp(:,2:4);
NXY0 = NXY0(:,2:4);
NXY00 = importdata('Nxy0.txt');
RXY = Ndsp(:,5:7);
for i=1:n
    if ismember(i,EELM)
        w = [1,1,1,1,1,1];
        e = E(i,:);
    else
        w = 1e-4*[1,1,1,1,1,1];
        e = E2(i,:);
    end
    NN = Enod(i,2:end); %第i个单元的节点号向量，节点1|2|3
    ENC = [NXY(NN,2),NXY(NN,3),NXY(NN,4)]';    %返回第i个单元的3个节点坐标，[xi|yi|zi]
    ENC0 = [NXY0(NN,1),NXY0(NN,2),NXY0(NN,3)]';
    ENC00 = [NXY00(NN,1),NXY00(NN,2),NXY00(NN,3)]';
    Rxy = [RXY(NN,1),RXY(NN,2),RXY(NN,3)]';%结点转动（全局坐标系,列）
    Rxy_iSTEP = [NOD_ROT(NN,1),NOD_ROT(NN,2),NOD_ROT(NN,3)]';
    %Thetaxy 结点转动（单元坐标系,列）
    %NUxy 结点位移（全局坐标系，列）
    [T0,len0] = transmat0_iBEAM3_CR(ENC0,ENC00,Rxy_iSTEP);
    [EK,EF] = cr_elementKF_iBEAM3(ENC,e,w,Rxy,T0,len0,2);
    % ------------------------------------------------------------------
    n1 = NN(1);
    n2 = NN(2);
    F = zeros(6*size(NXY,1),1);
    F(6*n1-5:6*n1,1) = EF(1:6);
    F(6*n2-5:6*n2,1) = EF(7:12);
    GF = GF + F;
    % ------------------------------------------------------------------
    for j=1:2   % 2 - number of node in an element
        for k=1:2   % 2 - number of node in an element
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