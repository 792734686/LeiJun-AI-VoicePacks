function [KT,f_in] = cr_elementKF_TR3(ENC,ENC0,Eepi,we,wk,wg,t,Rxy,Thetaxy,NUxy,i)
%ENC ----- 节点坐标，全局，行
%NUxy ----- 结点平动位移，全局，列
[Te,enc] = transmat_iTR3(ENC);
% -----------------------------------
[Te0,enc0] = transmat_iTR3(ENC0);
T = Te(1:3,1:3);
T0 = Te0(1:3,1:3);
% -----------------------------------
R1 = expm(skew(Rxy(:,1)));
R2 = expm(skew(Rxy(:,2)));
R3 = expm(skew(Rxy(:,3)));
Rd1 = T*R1*T0';
Rd2 = T*R2*T0';
Rd3 = T*R3*T0';
Thetaxy1 = invskew(logm(Rd1));
Thetaxy2 = invskew(logm(Rd2));
Thetaxy3 = invskew(logm(Rd3));
Thetaxy = [Thetaxy1,Thetaxy2,Thetaxy3];
% -----------------------------------
Nuxy(:,1) = [(enc(1,:)-enc0(1,:))';0];
Nuxy(:,2) = [(enc(2,:)-enc0(2,:))';0];
Nuxy(:,3) = [(enc(3,:)-enc0(3,:))';0];
% Nuxy 结点平动位移，局部，列
% enc 单元结点坐标，局部，行，xy
%% 计算G
A = triArea(enc(1,:),enc(2,:),enc(3,:));
x32 = enc(3,1)-enc(2,1);
y32 = enc(3,2)-enc(2,2);
x13 = enc(1,1)-enc(3,1);
y13 = enc(1,2)-enc(3,2);
x21 = enc(2,1)-enc(1,1);
y21 = enc(2,2)-enc(1,2);
l12 = x21;
G1u = 0.5/A*[0,0,x32;0,0,y32;0,-2*A/l12,0];
G2u = 0.5/A*[0,0,x13;0,0,y13;0,2*A/l12,0];
G3u = 0.5/A*[0,0,x21;0,0,y21;0,0,0];
G1 = [G1u,zeros(3)];
G2 = [G2u,zeros(3)];
G3 = [G3u,zeros(3)];
G = [G1,G2,G3];
%% 计算U,S
delta = eye(3);
U11 = (delta(1,1)-1/3)*eye(3);
U12 = (delta(1,2)-1/3)*eye(3);
U13 = (delta(1,3)-1/3)*eye(3);
U21 = (delta(2,1)-1/3)*eye(3);
U22 = (delta(2,2)-1/3)*eye(3);
U23 = (delta(2,3)-1/3)*eye(3);
U31 = (delta(3,1)-1/3)*eye(3);
U32 = (delta(3,2)-1/3)*eye(3);
U33 = (delta(3,3)-1/3)*eye(3);
X1 = [enc(1,:),0];
X2 = [enc(2,:),0];
X3 = [enc(3,:),0];
S1 = skew(X1);
S2 = skew(X2);
S3 = skew(X3);
%% 计算P
P11 = [U11+S1*G1u,zeros(3);-G1u,delta(1,1)*eye(3)];
P12 = [U12+S1*G2u,zeros(3);-G2u,delta(1,2)*eye(3)];
P13 = [U13+S1*G3u,zeros(3);-G3u,delta(1,3)*eye(3)];
P21 = [U21+S2*G1u,zeros(3);-G1u,delta(2,1)*eye(3)];
P22 = [U22+S2*G2u,zeros(3);-G2u,delta(2,2)*eye(3)];
P23 = [U23+S2*G3u,zeros(3);-G3u,delta(2,3)*eye(3)];
P31 = [U31+S3*G1u,zeros(3);-G1u,delta(3,1)*eye(3)];
P32 = [U32+S3*G2u,zeros(3);-G2u,delta(3,2)*eye(3)];
P33 = [U33+S3*G3u,zeros(3);-G3u,delta(3,3)*eye(3)];
P = [P11,P12,P13;P21,P22,P23;P31,P32,P33];
%% 计算H
% Thetaxy 结点转动向量，局部，列
    function H_theta = funH(theta)
        if norm(theta) == 0
            eta = 1/12;
        else
            eta = (1-0.5*norm(theta)*cot(0.5*norm(theta)))/(norm(theta)^2);
        end
        H_theta = eye(3)-0.5*skew(theta)+eta*skew(theta)^2;
    end
H1 = funH(Thetaxy(:,1));
H2 = funH(Thetaxy(:,2));
H3 = funH(Thetaxy(:,3));
H = eye(18);
H(4:6,4:6) = H1;
H(10:12,10:12) = H2;
H(16:18,16:18) = H3;
%% 计算Lambda
% Lambda11 = [(U11+S1*G1u)*T,zeros(3);-H1*G1*T,delta(1,1)*H1*T];
% Lambda12 = [(U12+S1*G2u)*T,zeros(3);-H1*G2*T,delta(1,2)*H1*T];
% Lambda13 = [(U13+S1*G3u)*T,zeros(3);-H1*G3*T,delta(1,3)*H1*T];
% Lambda21 = [(U21+S2*G1u)*T,zeros(3);-H2*G1*T,delta(2,1)*H2*T];
% Lambda22 = [(U22+S2*G2u)*T,zeros(3);-H2*G2*T,delta(2,2)*H2*T];
% Lambda23 = [(U23+S2*G3u)*T,zeros(3);-H2*G3*T,delta(2,3)*H2*T];
% Lambda31 = [(U31+S3*G1u)*T,zeros(3);-H3*G1*T,delta(3,1)*H3*T];
% Lambda32 = [(U32+S3*G2u)*T,zeros(3);-H3*G2*T,delta(3,2)*H3*T];
% Lambda33 = [(U33+S3*G3u)*T,zeros(3);-H3*G3*T,delta(3,3)*H3*T];
% Lambda = [Lambda11,Lambda12,Lambda13;
%           Lambda21,Lambda22,Lambda23;
%           Lambda31,Lambda32,Lambda33];
Lambda = H*P*Te;
%% 计算切线刚度阵与内力向量
% 局部坐标系刚度阵与内力向量
% Nuxy 单元结点平动位移，局部，列
Npxy = [Nuxy;Thetaxy];
ke = estiffmLC_TR3(ENC,we,wk,wg,t);
f = ke*Umat_to_lin(Npxy');
m1 = f(4:6);
m2 = f(10:12);
m3 = f(16:18);
f_p = P'*H'*f;
np1 = f_p(1:3);
mp1 = f_p(4:6);
np2 = f_p(7:9);
mp2 = f_p(10:12);
np3 = f_p(13:15);
mp3 = f_p(16:18);
Fnm = [skew(np1)',skew(mp1)',skew(np2)',skew(mp2)',skew(np3)',skew(mp3)']';
Fn = [skew(np1)',zeros(3),skew(np2)',zeros(3),skew(np3)',zeros(3)]';
% --------------------------------------------
    function Ma = funM(theta,m)
        Ha = funH(theta);
        if norm(theta) == 0
            eta = 1/12;
            mu = 1/360;
        else
            eta = (1-0.5*norm(theta)*cot(0.5*norm(theta)))/(norm(theta)^2);
            x = norm(theta);
            mu = (x^2+4*cos(x)+x*sin(x)-4)/(4*x^4*sin(0.5*x)^2);
        end
        Ma = (eta*((theta'*m)*eye(3)+theta*m'-2*m*theta')+...
            mu*skew(theta)^2*m*theta'-0.5*skew(m))*Ha;
    end
M1 = funM(Thetaxy(:,1),m1);
M2 = funM(Thetaxy(:,2),m2);
M3 = funM(Thetaxy(:,3),m3);
M = zeros(18);
M(4:6,4:6) = M1;
M(10:12,10:12) = M2;
M(16:18,16:18) = M3;
% ----------------------------------------------
KM = Lambda'*ke*Lambda;
KGR = -Te'*Fnm*G*Te;
KGP = -Te'*G'*Fn'*P*Te;
KGM = Te'*P'*M*P*Te;
KT = KM+KGR+KGP+KGM;
f_in = Te'*P'*H'*f;
% if i == 119
%     ke
%     Npxy
%     
% elseif i == 159
%     ENC
%     
% end
end