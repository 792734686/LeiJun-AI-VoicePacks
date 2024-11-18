function [T,b,b1,b2] = transmat_iBEAM3_CR(ENC,R1,R2,T0)
b1 = R1*T0'*[0,1,0]';
b2 = R2*T0'*[0,1,0]';
b = 0.5*(b1+b2);
X1 = ENC(:,1);
X2 = ENC(:,2);
n1 = (X2-X1)/norm(X2-X1);
n3 = cross(n1,b)/norm(cross(n1,b));
n2 = cross(n3,n1);
T = [n1,n2,n3]';
end