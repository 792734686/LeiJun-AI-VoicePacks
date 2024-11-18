function EF = eforcemLC_iTR3( ENC,Eepi,we,wk,wg,t )
% element stiffness matrix of 2-dimension and 4-node 
% ENC(4,2)--- element node coordinate   单元节点坐标矩阵，[xi|yi]
% t ----- thickness of element  
% D ----- elastic matrix for 2-dimension    
% ngs ---- number of gausss intigrate point     gauss积分点个数
% EK(6*4,6*4) ---  element stiffness matrix     单元刚度阵
eei = Eepi(1:3,1);
kei = Eepi(4:6,1);
gei = Eepi(7:8,1);
% -----------------------------------------------------------------
[~,enc] = transmat_iTR3(ENC);
% -----------------------------------------------------------------
[B1,JCB] = geom_iTR3( enc,-1/3,-1/3 );
[B2,~] = geom_iTR3( enc,2/3,-1/3 );
[B3,~] = geom_iTR3( enc,-1/3,2/3 );
[B4,~] = geom_iTR3( enc,0,0 );
Bm1 = B1(1:3,:);
Bb1 = B1(4:6,:);
Bs1 = B1(7:8,:);
EF1 = (we*(Bm1'*eei)+wk*t^2*(Bb1'*kei)+wg*2/3*(Bs1'*gei))*det(JCB)/24;
Bm2 = B2(1:3,:);
Bb2 = B2(4:6,:);
Bs2 = B2(7:8,:);
EF2 = (we*(Bm2'*eei)+wk*t^2*(Bb2'*kei)+wg*2/3*(Bs2'*gei))*det(JCB)/24;
Bm3 = B3(1:3,:);
Bb3 = B3(4:6,:);
Bs3 = B3(7:8,:);
EF3 = (we*(Bm3'*eei)+wk*t^2*(Bb3'*kei)+wg*2/3*(Bs3'*gei))*det(JCB)/24;
Bm4 = B4(1:3,:);
Bb4 = B4(4:6,:);
Bs4 = B4(7:8,:);
EF4 = (we*(Bm4'*eei)+wk*t^2*(Bb4'*kei)+wg*2/3*(Bs4'*gei))*det(JCB)*3/8;
EF = EF1+EF2+EF3+EF4;
% -----------------------------------------------------------------
% EF = Te'*EF;
end
