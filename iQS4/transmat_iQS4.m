function [Te,enc] = transmat_iQS4(ENC)
%Te --- 单元局部坐标到全局坐标的转换矩阵
%ENC(4,3) --- 单元节点全局坐标矩阵
%enc(4,2) --- 单元节点局部坐标矩阵
if size(ENC,2) == 2
    ENC(:,3) = 0;
end

X1 = ENC(1,:)';
X2 = ENC(2,:)';
X3 = ENC(3,:)';
X4 = ENC(4,:)';
%-----------------------------
A = X3-X1;
B = X4-X2;
c1 = 0.5*(X2+X1);
c2 = 0.5*(X3+X2);
c3 = 0.5*(X4+X3);
c4 = 0.5*(X1+X4);
d1 = norm(X2-X1);
d2 = norm(X3-X2);
d3 = norm(X4-X3);
d4 = norm(X1-X4);
%-----------------------------
n = cross(A,B)/norm(cross(A,B));
p = (A+B)/norm(A+B);
l = cross(p,n);
T = [l,p,n]';
Te(1:3,1:3) = T;
Te(4:6,4:6) = T;
Te(7:9,7:9) = T;
Te(10:12,10:12) = T;
Te(13:15,13:15) = T;
Te(16:18,16:18) = T;
Te(19:21,19:21) = T;
Te(22:24,22:24) = T;
%------------------------------------------------
C = (c1*d1+c2*d2+c3*d3+c4*d4)/(d1+d2+d3+d4);
x1 = dot(X1-C,l);
x2 = dot(X2-C,l);
x3 = dot(X3-C,l);
x4 = dot(X4-C,l);
y1 = dot(X1-C,p);
y2 = dot(X2-C,p);
y3 = dot(X3-C,p);
y4 = dot(X4-C,p);
enc = [x1,y1;
       x2,y2;
       x3,y3;
       x4,y4];
end