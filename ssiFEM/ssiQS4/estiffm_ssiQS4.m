function [EK,B] = estiffm_ssiQS4( ENC,ngs,we,wg,t )
% element stiffness matrix of 2-dimension and 4-node 
% ENC(4,2)--- element node coordinate
% t ----- thickness of element  板厚
% D ----- elastic matrix for 2-dimension
% ngs ---- number of gausss intigrate point
% EK(6*4,6*4) ---  element stiffness matrix
EK=zeros(6*4,6*4);
[ GSP,GSW ] = gauspw( ngs );
[Te,enc] = transmat_iQS4(ENC);
for i = 1:ngs
    for j = 1:ngs
        xi = GSP(i);
        eta = GSP(j);
        wi = GSW(i);
        wj = GSW(j);
%         [~,DSHP] = shape_2d4n (xi,eta,ENC);
        [B,JCB] = geom_ssiQS4( enc,xi,eta,t );
        Bm = B(1:3,:);
        Bs = B(4:5,:);
        EK = EK + (we*(Bm'*Bm)+wg*(Bs'*Bs))*det(JCB)*wi*wj;
    end
end
EK = Te'*EK*Te;
% EK = t*EK;
end
