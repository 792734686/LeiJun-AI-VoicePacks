function [B,JCB] = geom_LargeDeform_iQS4( ENC,xi,eta,elm_ui )
% geometric matrix for an element of iQS4——输出单元应变矩阵与雅可比矩阵
% ENC(4,2) --- element node coordinate——单元节点坐标，4行2列，x坐标|y坐标
% DSHP(6,4) --- derivative matrix of shape function——形函数微分矩阵
% xi, eta --- local coordinate——局部坐标
% B(8,24) --- element geometric matrix,8*24
% JCB(2,2) --- Jacobian matrix
% elm_ui(4,6) --- 
[SHP,DSHP] = shape_iQS4 (xi,eta,ENC);
NDSHP = DSHP(1:2,:);
MDSHP = DSHP(3:4,:);
LDSHP = DSHP(5:6,:);
JCB = NDSHP*ENC;
% IJCB = inv(JCB);
NDSHPG = JCB\NDSHP;
MDSHPG = JCB\MDSHP;
LDSHPG = JCB\LDSHP;
%%--------------------------------------
theta_xi = elm_ui(:,4);
theta_yi = elm_ui(:,5);
theta_1 = SHP*theta_yi;
theta_2 = -SHP*theta_xi;
%%--------------------------------------
B1 = [NDSHPG(1,1),0,0,0,0.5*theta_1*SHP(1),LDSHPG(1,1);
      0,NDSHPG(2,1),0,-0.5*theta_2*SHP(1),0,MDSHPG(2,1);
      NDSHPG(2,1),NDSHPG(1,1),0,-0.5*theta_1*SHP(1),0.5*theta_2*SHP(1),LDSHPG(2,1)+MDSHPG(1,1);
      0,0,0,0,NDSHPG(1,1),0;
      0,0,0,-NDSHPG(2,1),0,0;
      0,0,0,-NDSHPG(1,1),NDSHPG(2,1),0;
      0,0,NDSHPG(1,1),-LDSHPG(1,1),SHP(1)-MDSHPG(1,1),0;
      0,0,NDSHPG(2,1),-(LDSHPG(2,1)+SHP(1)),-MDSHPG(2,1),0];
B2 = [NDSHPG(1,2),0,0,0,0.5*theta_1*SHP(2),LDSHPG(1,2);
      0,NDSHPG(2,2),0,-0.5*theta_2*SHP(2),0,MDSHPG(2,2);
      NDSHPG(2,2),NDSHPG(1,2),0,-0.5*theta_1*SHP(2),0.5*theta_2*SHP(2),LDSHPG(2,2)+MDSHPG(1,2);
      0,0,0,0,NDSHPG(1,2),0;
      0,0,0,-NDSHPG(2,2),0,0;
      0,0,0,-NDSHPG(1,2),NDSHPG(2,2),0;
      0,0,NDSHPG(1,2),-LDSHPG(1,2),SHP(2)-MDSHPG(1,2),0;
      0,0,NDSHPG(2,2),-(LDSHPG(2,2)+SHP(2)),-MDSHPG(2,2),0];
B3 = [NDSHPG(1,3),0,0,0,0.5*theta_1*SHP(3),LDSHPG(1,3);
      0,NDSHPG(2,3),0,-0.5*theta_2*SHP(3),0,MDSHPG(2,3);
      NDSHPG(2,3),NDSHPG(1,3),0,-0.5*theta_1*SHP(3),0.5*theta_2*SHP(3),LDSHPG(2,3)+MDSHPG(1,3);
      0,0,0,0,NDSHPG(1,3),0;
      0,0,0,-NDSHPG(2,3),0,0;
      0,0,0,-NDSHPG(1,3),NDSHPG(2,3),0;
      0,0,NDSHPG(1,3),-LDSHPG(1,3),SHP(3)-MDSHPG(1,3),0;
      0,0,NDSHPG(2,3),-(LDSHPG(2,3)+SHP(3)),-MDSHPG(2,3),0];
B4 = [NDSHPG(1,4),0,0,0,0.5*theta_1*SHP(4),LDSHPG(1,4);
      0,NDSHPG(2,4),0,-0.5*theta_2*SHP(4),0,MDSHPG(2,4);
      NDSHPG(2,4),NDSHPG(1,4),0,-0.5*theta_1*SHP(4),0.5*theta_2*SHP(4),LDSHPG(2,4)+MDSHPG(1,4);
      0,0,0,0,NDSHPG(1,4),0;
      0,0,0,-NDSHPG(2,4),0,0;
      0,0,0,-NDSHPG(1,4),NDSHPG(2,4),0;
      0,0,NDSHPG(1,4),-LDSHPG(1,4),SHP(4)-MDSHPG(1,4),0;
      0,0,NDSHPG(2,4),-(LDSHPG(2,4)+SHP(4)),-MDSHPG(2,4),0];
B = [B1, B2, B3, B4];


end

