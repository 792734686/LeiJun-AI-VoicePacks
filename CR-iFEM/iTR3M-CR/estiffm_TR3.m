 function EK = estiffm_TR3(ENC,we,wk,wg,t)
% element stiffness matrix of 2-dimension and 4-node 
% ENC(4,2)--- element node coordinate
% t ----- thickness of element  板厚
% D ----- elastic matrix for 2-dimension
% ngs ---- number of gausss intigrate point
% EK(6*4,6*4) ---  element stiffness matrix
[Te,enc] = transmat_iTR3(ENC);
E = 6.825e7;
mu = 0.3;
Dm = E/(1-mu^2)*[1,mu,0;mu,1,0;0,0,0.5*(1-mu)];
Db = Dm;
Ds = E/(1+mu)*[1,0;0,1];
% -----------------------------------------------------------------
we = 1;
wk = 1;
wg = 1;
% -------------------------------------
[B1,JCB] = geom_iTR3( enc,0,0 );
[B2,~] = geom_iTR3( enc,1/5-1/3,1/5-1/3 );
[B3,~] = geom_iTR3( enc,3/5-1/3,1/5-1/3 );
[B4,~] = geom_iTR3( enc,1/5-1/3,3/5-1/3 );
Bm1 = B1(1:3,:);
Bb1 = B1(4:6,:);
Bs1 = B1(7:8,:);
EK1 = (we*t*(Bm1'*Dm*Bm1)+wk*t^3/12*(Bb1'*Db*Bb1)+wg*t*5/12*(Bs1'*Ds*Bs1))*det(JCB)*-27/96;
Bm2 = B2(1:3,:);
Bb2 = B2(4:6,:);
Bs2 = B2(7:8,:);
EK2 = (we*t*(Bm2'*Dm*Bm2)+wk*t^3/12*(Bb2'*Db*Bb2)+wg*t*5/12*(Bs2'*Ds*Bs2))*det(JCB)*25/96;
Bm3 = B3(1:3,:);
Bb3 = B3(4:6,:);
Bs3 = B3(7:8,:);
EK3 = (we*t*(Bm3'*Dm*Bm3)+wk*t^3/12*(Bb3'*Db*Bb3)+wg*t*5/12*(Bs3'*Ds*Bs3))*det(JCB)*25/96;
Bm4 = B4(1:3,:);
Bb4 = B4(4:6,:);
Bs4 = B4(7:8,:);
EK4 = (we*t*(Bm4'*Dm*Bm4)+wk*t^3/12*(Bb4'*Db*Bb4)+wg*t*5/12*(Bs4'*Ds*Bs4))*det(JCB)*25/96;
EK = EK1+EK2+EK3+EK4;
% -----------------------------------------------------------------
EK = Te'*EK*Te;
% EK = t*EK;
end
