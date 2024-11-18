function EF = eforcem_iQS4( ENC,ngs,Eepi,we,wk,wg,t )
% element stiffness matrix of 2-dimension and 4-node 
% ENC(4,2)--- element node coordinate   单元节点坐标矩阵，[xi|yi]
% t ----- thickness of element  
% D ----- elastic matrix for 2-dimension    
% ngs ---- number of gausss intigrate point     gauss积分点个数
% EK(6*4,6*4) ---  element stiffness matrix     单元刚度阵
EF=zeros(24,1);
[ GSP,GSW ] = gauspw( ngs );
eei = Eepi(1:3,1);
kei = Eepi(4:6,1);
gei = Eepi(7:8,1);
[Te,enc] = transmat_iQS4(ENC);
for i = 1:ngs
    for j = 1:ngs
        xi = GSP(i);
        eta = GSP(j);
        wi = GSW(i);
        wj = GSW(j);
%         [~,DSHP] = shape_iQS4 (xi,eta,ENC);
        [B,JCB] = geom_iQS4( enc,xi,eta );
        Bm = B(1:3,:);
        Bb = B(4:6,:);
        Bs = B(7:8,:);
        EF = EF + (we*Bm'*eei+wk*t^2*Bb'*kei+wg*2/3*Bs'*gei)*det(JCB)*wi*wj;
    end
end
EF = Te'*EF;
end
