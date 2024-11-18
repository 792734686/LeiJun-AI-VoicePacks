function [SHP,DSHP] = shape_iQS4 (s,t,ENC)
% shape function for an element of 2-dimension and 4-node——返回形函数矩阵与形函数微分矩阵
% xi,eta --- local coordinate——局部坐标
% SHP(4) --- shape functions——形函数矩阵(未修改为QS4形函数)
% DSHP(6,4) --- derivative of shape function——形函数的微分矩阵
% ENC(4,2) --- 单元节点坐标，x|y
% ---------------------------
a1 = ENC(1,1) - ENC(2,1);
a2 = ENC(2,1) - ENC(3,1);
a3 = ENC(3,1) - ENC(4,1);
a4 = ENC(4,1) - ENC(1,1);
b1 = ENC(2,2) - ENC(1,2);
b2 = ENC(3,2) - ENC(2,2);
b3 = ENC(4,2) - ENC(3,2);
b4 = ENC(1,2) - ENC(4,2);
s1 = -1;
s2 = 1;
s3 = 1;
s4 = -1;
t1 = -1;
t2 = -1;
t3 = 1;
t4 = 1;
% ---------------------------
N1s = -0.25*(1-t);
N2s = 0.25*(1-t);
N3s = 0.25*(1+t);
N4s = -0.25*(1+t);

N1t = -0.25*(1-s);
N2t = -0.25*(1+s);
N3t = 0.25*(1+s);
N4t = 0.25*(1-s);
% ---------------------------
N5s = -1/8*s*(1+t1*t);
N6s = 1/16*s2*(1-t^2);
N7s = -1/8*s*(1+t3*t);
N8s = 1/16*s4*(1-t^2);

N5t = 1/16*t1*(1-s^2);
N6t = -1/8*t*(1+s2*s);
N7t = 1/16*t3*(1-s^2);
N8t = -1/8*t*(1+s4*s);
% ---------------------------
M1s = a4*N8s-a1*N5s;
M2s = a1*N5s-a2*N6s;
M3s = a2*N6s-a3*N7s;
M4s = a3*N7s-a4*N8s;

M1t = a4*N8t-a1*N5t;
M2t = a1*N5t-a2*N6t;
M3t = a2*N6t-a3*N7t;
M4t = a3*N7t-a4*N8t;
% ---------------------------
L1s = b4*N8s-b1*N5s;
L2s = b1*N5s-b2*N6s;
L3s = b2*N6s-b3*N7s;
L4s = b3*N7s-b4*N8s;

L1t = b4*N8t-b1*N5t;
L2t = b1*N5t-b2*N6t;
L3t = b2*N6t-b3*N7t;
L4t = b3*N7t-b4*N8t;
% ---------------------------
SHP(1) = 0.25*(1-s)*(1-t);
SHP(2) = 0.25*(1+s)*(1-t);
SHP(3) = 0.25*(1+s)*(1+t);
SHP(4) = 0.25*(1-s)*(1+t);
%----------------------------
DSHP = [N1s,N2s,N3s,N4s;
        N1t,N2t,N3t,N4t;
        M1s,M2s,M3s,M4s;
        M1t,M2t,M3t,M4t;
        L1s,L2s,L3s,L4s;
        L1t,L2t,L3t,L4t];
%----------------------------
end
