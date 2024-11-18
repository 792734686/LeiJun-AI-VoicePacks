function [EK,B] = FEM_estiffm_iQS4( ENC,ngs,we,wk,wg,t )
% element stiffness matrix of 2-dimension and 4-node 
% ENC(4,2)--- element node coordinate
% t ----- thickness of element  板厚
% D ----- elastic matrix for 2-dimension
% ngs ---- number of gausss intigrate point
% EK(6*4,6*4) ---  element stiffness matrix
EK=zeros(6*4,6*4);
[ GSP,GSW ] = gauspw( ngs );
[Te,enc] = transmat_iQS4(ENC);
%%
%材料属性
E = 210000;
mu = 0.3;
G = E/(2*(1+mu));
%%
De = E/(1-mu^2)*[1,mu,0;mu,1,0;0,0,0.5*(1-mu)];
Dk = E*t^3/(12*(1-mu^2))*[1,mu,0;mu,1,0;0,0,0.5*(1-mu)];
Dg = G*t/1.2*[1,0;0,1];
for i = 1:ngs
    for j = 1:ngs
        xi = GSP(i);
        eta = GSP(j);
        wi = GSW(i);
        wj = GSW(j);
%         [~,DSHP] = shape_2d4n (xi,eta,ENC);
        [B,JCB] = geom_iQS4( enc,xi,eta );
        Be = B(1:3,:);
        Bk = B(4:6,:);
        Bg = B(7:8,:);
        EK = EK + (Be'*De*Be+Bk'*Dk*Bk+Bg'*Dg*Bg)*det(JCB)*wi*wj;
    end
end
EK = Te'*EK*Te;
% EK = t*EK;
end
