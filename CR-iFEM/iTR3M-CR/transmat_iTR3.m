function [Te,enc] = transmat_iTR3(ENC)
%Te --- 单元局部坐标到全局坐标的转换矩阵
%ENC(3,3) --- 单元节点全局坐标矩阵
%enc(3,2) --- 单元节点局部坐标矩阵
if size(ENC,2) == 2
    ENC(:,3) = 0;
end

X1 = ENC(1,:)';
X2 = ENC(2,:)';
X3 = ENC(3,:)';
X0 = 1/3*(X1+X2+X3);
X10 = X1-X0;
X20 = X2-X0;
X30 = X3-X0;
%-----------------------------
L12 = X2-X1;
L13 = X3-X1;
L23 = X3-X2;
%-----------------------------
l = L12/norm(L12);
n = cross(L12,L13)/norm(cross(L12,L23));
p = cross(n,l);
% n = cross(A,B)/norm(cross(A,B));
% p = (A+B)/norm(A+B);
% l = cross(p,n);
T = [l,p,n]';
Te(1:3,1:3) = T;
Te(4:6,4:6) = T;
Te(7:9,7:9) = T;
Te(10:12,10:12) = T;
Te(13:15,13:15) = T;
Te(16:18,16:18) = T;
%------------------------------------------------
% x1 = 0;
% y1 = 0;
% x2 = norm(L12);
% y2 = 0;
% x3 = dot(L13,l);
% y3 = dot(L13,p);
x1 = dot(X10,l);
y1 = dot(X10,p);
x2 = dot(X20,l);
y2 = dot(X20,p);
x3 = dot(X30,l);
y3 = dot(X30,p);
enc = [x1,y1;
       x2,y2;
       x3,y3];
end