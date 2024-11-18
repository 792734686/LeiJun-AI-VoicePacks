function Te = transmat_iBEAM3(ENC)
X1 = ENC(:,1);
X2 = ENC(:,2);
n1 = (X2-X1)/norm(X2-X1);
if cross(n1,[0,0,-1]') == 0
    n2 = [1,0,0]';
else
    n2 = [0,0,-1]';
end
n3 = cross(n1,n2)/norm(cross(n1,n2));
n2 = cross(n3,n1);
T = [n1,n2,n3]';
Te(1:3,1:3) = T;
Te(4:6,4:6) = T;
Te(7:9,7:9) = T;
Te(10:12,10:12) = T;
end