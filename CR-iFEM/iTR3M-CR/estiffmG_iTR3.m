 function EK = estiffmG_iTR3( ENC,we,wk,wg,t )
% element stiffness matrix of 2-dimension and 4-node 
% ENC(4,2)--- element node coordinate
% t ----- thickness of element  板厚
% D ----- elastic matrix for 2-dimension
% ngs ---- number of gausss intigrate point
% EK(6*4,6*4) ---  element stiffness matrix
[Te,enc] = transmat_iTR3(ENC);
% -----------------------------------------------------------------
[B1,JCB] = geom_iTR3( enc,0,0 );
[B2,~] = geom_iTR3( enc,1,0 );
[B3,~] = geom_iTR3( enc,0,1 );
[B4,~] = geom_iTR3( enc,1/3,1/3 );
Bm1 = B1(1:3,:);
Bb1 = B1(4:6,:);
Bs1 = B1(7:8,:);
EK1 = (we*(Bm1'*Bm1)+wk*t^2*(Bb1'*Bb1)+wg*2/3*(Bs1'*Bs1))*det(JCB)/24;
Bm2 = B2(1:3,:);
Bb2 = B2(4:6,:);
Bs2 = B2(7:8,:);
EK2 = (we*(Bm2'*Bm2)+wk*t^2*(Bb2'*Bb2)+wg*2/3*(Bs2'*Bs2))*det(JCB)/24;
Bm3 = B3(1:3,:);
Bb3 = B3(4:6,:);
Bs3 = B3(7:8,:);
EK3 = (we*(Bm3'*Bm3)+wk*t^2*(Bb3'*Bb3)+wg*2/3*(Bs3'*Bs3))*det(JCB)/24;
Bm4 = B4(1:3,:);
Bb4 = B4(4:6,:);
Bs4 = B4(7:8,:);
EK4 = (we*(Bm4'*Bm4)+wk*t^2*(Bb4'*Bb4)+wg*2/3*(Bs4'*Bs4))*det(JCB)*3/8;
EK = EK1+EK2+EK3+EK4;
% -----------------------------------------------------------------
% EK = Te'*EK*Te;
% EK = t*EK;
end
