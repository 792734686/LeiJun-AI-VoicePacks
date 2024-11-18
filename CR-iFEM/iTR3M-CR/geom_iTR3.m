function [B,JCB] = geom_iTR3( ENC,xi,eta )
% geometric matrix for an element of iQS4——输出单元应变矩阵与雅可比矩阵
% ENC(4,2) --- element node coordinate——单元节点坐标，4行2列，x坐标|y坐标
% DSHP(6,4) --- derivative matrix of shape function——形函数微分矩阵
% xi, eta --- local coordinate——局部坐标
% B(8,24) --- element geometric matrix,8*24
% JCB(2,2) --- Jacobian matrix
[SHP,DSHP] = shape_iTR3 (xi,eta,ENC);
LDSHP = DSHP(1:2,:);
NuDSHP = DSHP(3:4,:);
NvDSHP = DSHP(5:6,:);
JCB = LDSHP*ENC;
% IJCB = inv(JCB);
LDSHPG = JCB\LDSHP;
NuDSHPG = JCB\NuDSHP;
NvDSHPG = JCB\NvDSHP;
B1 = [LDSHPG(1,1),0,0,0,0,NuDSHPG(1,1);
      0,LDSHPG(2,1),0,0,0,NvDSHPG(2,1);
      LDSHPG(2,1),LDSHPG(1,1),0,0,0,NuDSHPG(2,1)+NvDSHPG(1,1);
      0,0,0,0,LDSHPG(1,1),0;
      0,0,0,-LDSHPG(2,1),0,0;
      0,0,0,-LDSHPG(1,1),LDSHPG(2,1),0;
      0,0,LDSHPG(1,1),-NuDSHPG(1,1),SHP(1)-NvDSHPG(1,1),0;
      0,0,LDSHPG(2,1),-(NuDSHPG(2,1)+SHP(1)),-NvDSHPG(2,1),0];
B2 = [LDSHPG(1,2),0,0,0,0,NuDSHPG(1,2);
      0,LDSHPG(2,2),0,0,0,NvDSHPG(2,2);
      LDSHPG(2,2),LDSHPG(1,2),0,0,0,NuDSHPG(2,2)+NvDSHPG(1,2);
      0,0,0,0,LDSHPG(1,2),0;
      0,0,0,-LDSHPG(2,2),0,0;
      0,0,0,-LDSHPG(1,2),LDSHPG(2,2),0;
      0,0,LDSHPG(1,2),-NuDSHPG(1,2),SHP(2)-NvDSHPG(1,2),0;
      0,0,LDSHPG(2,2),-(NuDSHPG(2,2)+SHP(2)),-NvDSHPG(2,2),0];
B3 = [LDSHPG(1,3),0,0,0,0,NuDSHPG(1,3);
      0,LDSHPG(2,3),0,0,0,NvDSHPG(2,3);
      LDSHPG(2,3),LDSHPG(1,3),0,0,0,NuDSHPG(2,3)+NvDSHPG(1,3);
      0,0,0,0,LDSHPG(1,3),0;
      0,0,0,-LDSHPG(2,3),0,0;
      0,0,0,-LDSHPG(1,3),LDSHPG(2,3),0;
      0,0,LDSHPG(1,3),-NuDSHPG(1,3),SHP(3)-NvDSHPG(1,3),0;
      0,0,LDSHPG(2,3),-(NuDSHPG(2,3)+SHP(3)),-NvDSHPG(2,3),0];
B = [B1, B2, B3];
B(3,:) = B(3,:);
end

