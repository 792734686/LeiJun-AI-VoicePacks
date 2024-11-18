function [GK,GF_in] = cr_globalKF_TR3(NXY,NXY0,ELEM,Eep,Eep2,EELM,t,Ndsp,Theta)
n=size(NXY,1);
GK = zeros(6*n,6*n);
GF_in = zeros(6*n,1);
n=size(ELEM,1);
NXY(:,2:4) = NXY(:,2:4)+Ndsp(:,2:4);
NXY0 = NXY0(:,2:4);
RXY = Ndsp(:,5:7);
ThetaXY = Theta';
for i=1:n
    if ismember(i,EELM)
        we = 1;
        wk = 1;
        wg = 1e-4;
        Eepi = Eep(i,:)';
    else
        we = 1e-4;
        wk = 1e-4;
        wg = 1e-4;
        Eepi = Eep2(i,:)';
    end
    NN = ELEM(i,2:end); %第i个单元的节点号向量，节点1|2|3
    ENC = [NXY(NN,2),NXY(NN,3),NXY(NN,4)];    %返回第i个单元的3个节点坐标，[xi|yi|zi]
    ENC0 = [NXY0(NN,1),NXY0(NN,2),NXY0(NN,3)];
    Rxy = [RXY(NN,1),RXY(NN,2),RXY(NN,3)]';%结点转动（全局坐标系,列）
    %Thetaxy 结点转动（单元坐标系,列）
    Thetaxy = [ThetaXY(NN,1),ThetaXY(NN,2),ThetaXY(NN,3)]';
    %NUxy 结点位移（全局坐标系，列）
    NUxy = [Ndsp(NN,2),Ndsp(NN,3),Ndsp(NN,4)]';
    [EK,f_in] = cr_elementKF_TR3(ENC,ENC0,Eepi,we,wk,wg,t,Rxy,Thetaxy,NUxy,i);
    % ------------------------------------------------------------------
    n1 = NN(1);
    n2 = NN(2);
    n3 = NN(3);
    F_in = zeros(6*size(NXY,1),1);
    F_in(6*n1-5:6*n1) = f_in(1:6);
    F_in(6*n2-5:6*n2) = f_in(7:12);
    F_in(6*n3-5:6*n3) = f_in(13:18);
    GF_in = GF_in+F_in;
    % ------------------------------------------------------------------
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
% [GF(6*(21-1)+3),GF(6*(21*2-1)+3),GF(6*(21*3-1)+3),GF(6*(21*4-1)+3),GF(6*(21*5-1)+3),GF(6*(21*6-1)+3),GF(6*(21*7-1)+3)]
end