function Ndsp = Slove_6DOF(KS,SP)
% structural stiffness equation solution function
% KS(6N, 6N) -- structural stiffness matrix     总刚6n*6n
% SP(2N, 1) -- nodal load vector    总载荷向量6n*1
% cons(nc, 3) -- node constraint matrix     边界条件矩阵，节点号|自由度方向|位移量
% cons(i, 3) -- [node, 1(x) 2(y) 3(z) 4(rx) 5(ry) 6(rz), value]
% Ndsp(n, 7) --- 节点位移矩阵，节点号|x|y|z|rx|ry|rz
KS1 = KS;   %总刚
SP1 = SP;   %载荷
KS1 = sparse(KS1);
NDSP1 = KS1\SP1;
nn = size(NDSP1,1)/6;
Ndsp = zeros(nn,7);
for i = 1:nn
    Ndsp(i,1) = i;
    Ndsp(i,2) = NDSP1(6*i-5, 1);
    Ndsp(i,3) = NDSP1(6*i-4, 1);
    Ndsp(i,4) = NDSP1(6*i-3, 1);
    Ndsp(i,5) = NDSP1(6*i-2, 1);
    Ndsp(i,6) = NDSP1(6*i-1, 1);
    Ndsp(i,7) = NDSP1(6*i, 1);
end
end
