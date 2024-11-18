function [EK,B] = estiffm_LargeDeform_iQS4( ENC,ELM_Ui,ngs,we,wk,wg,t,Eepi )
% element stiffness matrix of 2-dimension and 4-node 
% ENC(4,2)--- element node coordinate
% t ----- thickness of element  板厚
% D ----- elastic matrix for 2-dimension
% ngs ---- number of gausss intigrate point
% EK(6*4,6*4) ---  element stiffness matrix
EK=zeros(6*4,6*4);
[ GSP,GSW ] = gauspw( ngs );
[Te,enc] = transmat_iQS4(ENC);
T = Te(1:3,1:3);
elm_ui = zeros(4,6);
elm_ui(1:4,1:3) = (T*ELM_Ui(:,1:3)')';
elm_ui(1:4,4:6) = (T*ELM_Ui(:,4:6)')';
for i = 1:ngs
    for j = 1:ngs
        xi = GSP(i);
        eta = GSP(j);
        wi = GSW(i);
        wj = GSW(j);
%         [~,DSHP] = shape_2d4n (xi,eta,ENC);
        [B,JCB] = geom_LargeDeform_iQS4( enc,xi,eta,elm_ui );
        Beta = matrix_Beta_iQS4(enc,xi,eta,Eepi);
        Bm = B(1:3,:);
        Bb = B(4:6,:);
        Bs = B(7:8,:);
        EK = EK + (we*(Bm'*Bm-Beta)+wk*t^2*(Bb'*Bb)+wg*(Bs'*Bs))*det(JCB)*wi*wj;
    end
end
EK = Te'*EK*Te;
% EK = t*EK;
end
