clear;
tol = 1e-4;
maxiter = 1000;
iterations = 0;
%输入矩阵Nxy,Enod:
%Enod:[单元号|节点1|节点2|节点3|节点4]
%Nxy:[节点号|x坐标|y坐标|z坐标]
%--------------------------------------------
Enod = importdata('Enod.rpt').data;
Enod = [(1:size(Enod,1))',Enod];
%--------------------------------------------
Nxy = importdata('Nxy.rpt').data;
Nxy = Nxy(:,1:4);
%--------------------------------------------
G13 = zeros(size(Enod,1),1);
G23 = G13;
Ndsp = zeros(size(Nxy,1),7);

while iterations < maxiter
    t = 20;      %板厚
    [G13,G23] = StrainG(Enod, Ndsp, Nxy);
    EELM = Enod(:,1);   %输入测量应变的单元
    Eep = EptranPLN(t);    %转换为8应变分量Eep(单元数1,8)
    Eep(:,7) = G13;
    Eep(:,8) = G23;
    Eep2 = EptranPLN(t);
    %-----------------------------------------------
    [GK,B] = gstiffm_iQS4(Nxy,Enod,EELM,t);
    GF = gforcem_iQS4(Nxy, Enod, Eep,Eep2,EELM,t);      %等效载荷向量
    %----------------------------------------------
    FIX = importdata('fix.rpt');
    [GK,GF] = boundary(FIX,[1,2,3,4,5,6],GK,GF);
    %----------------------------------------------
    Ndsp0 = Ndsp;
    Ndsp = SloveiQS4n(GK,GF);   %返回节点位移值与约束反力*
    
    tolerance = tol*max(max(Ndsp(:,2:7)));
    if norm(Ndsp0(:,2:7)-Ndsp(:,2:7)) <= tolerance
        disp('solver succeed')
        break
    end
    iterations = iterations+1;
end

if iterations == maxiter
    disp('too many steps for solver')
else
    disp('solve succeed')
end