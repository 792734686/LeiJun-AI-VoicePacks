function [KTP,EF] = cr_elementKF_iBEAM3(ENC,e,w,Rxy,T0,len0,ngs)

len = norm(ENC(:,1)-ENC(:,2));

    function R = funR(array)
        a = norm(array);
        if a == 0
            R = eye(3);
        else
            R = eye(3)+sin(a)/a*skew(array)+...
                0.5*(sin(0.5*a)/(0.5*a))^2*skew(array)^2;
        end
    end
% Ri ---- 结点转动矩阵，全局
R1 = funR(Rxy(:,1));
R2 = funR(Rxy(:,2));
% T0,T ---- 坐标变换矩阵
[T,b_ary,b1_ary,b2_ary] = transmat_iBEAM3_CR(ENC,R1,R2,T0);
oo = zeros(3);
Te = [T,oo,oo,oo;
      oo,T,oo,oo;
      oo,oo,T,oo;
      oo,oo,oo,T];
% ------------------
R1e = T*R1*T0';
R2e = T*R2*T0';
Theta1 = Array_by_Quaternions(R1e);
Theta2 = Array_by_Quaternions(R2e);
NU1 = [0,0,0]';
NU2 = [len-len0,0,0]';

%% 计算参考系
B_ary = T*b_ary;
B1_ary = T*b1_ary;
B2_ary = T*b2_ary;
b1 = B_ary(1);
b2 = B_ary(2);
b11 = B1_ary(1);
b12 = B1_ary(2);
b21 = B2_ary(1);
b22 = B2_ary(2);
eta = b1/b2;
eta11 = b11/b2;
eta12 = b12/b2;
eta21 = b21/b2;
eta22 = b22/b2;

%% 计算G,S
S = [0,0,0,1,0,0,0,0,0,1,0,0;
     0,0,0,0,1,0,0,0,-len,0,1,0;
     0,0,0,0,0,1,0,len,0,0,0,1]';
G1 = [0,0,eta/len;
      0,0,1/len;
      0,-1/len,0];
G2 = [eta12/2,-eta11/2,0;
      0,0,0;
      0,0,0];
G3 = [0,0,-eta/len;
      0,0,-1/len;
      0,1/len,0];
G4 = [eta22/2,-eta21/2,0;
      0,0,0;
      0,0,0];
G = [G1,G2,G3,G4];
%% 计算P
P = eye(12)-S*G;
%% 计算H,Lambda
    function H_theta = funH(theta)
        a = norm(theta);
        if a <= 0.05
            eta = 1/12 + 1/720*a^2 + 1/30240*a^4 + 1/1209600*a^6;
        else
            eta = (1-0.5*a*cot(0.5*a))/(a^2);
        end
        H_theta = eye(3)-0.5*skew(theta)+eta*skew(theta)^2;
    end
H1 = funH(Theta1);
H2 = funH(Theta2);
H(1:3,1:3) = eye(3);
H(4:6,4:6) = H1;
H(7:9,7:9) = eye(3);
H(10:12,10:12) = H2;
Lambda = H*P*Te;
%% 切线刚度阵与内力矢量
EK = estiff_iBEAM3(w,ENC,ngs);
ef = eforce_iBEAM3(w,e,ENC,ngs);
%凝聚
KRR = EK(1:12,1:12);
KRO = EK(1:12,13:14);
KOO = EK(13:14,13:14);
FR = ef(1:12);
FO = ef(13:14);
EK = KRR-KRO*KOO^(-1)*KRO';
ef = FR-KRO*KOO^(-1)*FO;
% 内力矢量
f_in = EK*[NU1;Theta1;NU2;Theta2];
f_p = P'*H'*f_in;
np1 = f_p(1:3);
mp1 = f_p(4:6);
np2 = f_p(7:9);
mp2 = f_p(10:12);
m1 = f_in(4:6);
m2 = f_in(10:12);
Fn = [skew(np1)',zeros(3),skew(np2)',zeros(3)]';
Fnm = [skew(np1)',skew(mp1)',skew(np2)',skew(mp2)']';

    function Ma = funM(theta,m)
        a = norm(theta);
        Ha = funH(theta);
        if norm(theta) <= 0.1
            eta = 1/12 + 1/720*a^2 + 1/30240*a^4 + 1/1209600*a^6;
            mu = 1/360 + 1/7560*a^2 + 1/201600*a^4 + 1/5987520*a^6;
        else
            eta = (1-0.5*a*cot(0.5*a))/(a^2);
            mu = (a^2+4*cos(a)+a*sin(a)-4)/(4*a^4*sin(0.5*a)^2);
        end
        Ma = (eta*((theta'*m)*eye(3)+theta*m'-2*m*theta')+...
            mu*skew(theta)^2*m*theta'-0.5*skew(m))*Ha;
    end

M1 = funM(Theta1,m1);
M2 = funM(Theta2,m2);
M = zeros(12);
M(4:6,4:6) = M1;
M(10:12,10:12) = M2;
% 材料刚度阵
KM = Lambda'*EK*Lambda;
% 几何刚度阵
KGR = -Te'*Fnm*G*Te;
KGP = -Te'*G'*Fn'*P*Te;
KGM = Te'*P'*M*P*Te;
% 载荷刚度阵
n1e = ef(1:3);
m1e = ef(4:6);
n2e = ef(7:9);
m2e = ef(10:12);
Fp = [skew(n1e);skew(m1e);skew(n2e);skew(m2e)];
KL = Te'*Fp*G*Te;
% 刚度阵
KTP = KM+KGR+KGP+KGM+KL;
EF = Te'*(ef-f_p);

end