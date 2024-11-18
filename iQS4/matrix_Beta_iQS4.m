function B = matrix_Beta_iQS4(ENC,xi,eta,Eepi)
[N,SHP] = basic_shape_iQS4 (xi,eta,ENC);
epsilon11 = Eepi(1);
epsilon22 = Eepi(2);
epsilon12 = Eepi(3);
Lamb1 = [zeros(3,5);
         0,0,0,-0.5*SHP(1)*epsilon12,-0.5*SHP(1)*epsilon22;
         0,0,0,0.5*SHP(1)*epsilon11,0.5*SHP(1)*epsilon12;
         0,0,0,0,0];
Lamb2 = [zeros(3,5);
         0,0,0,-0.5*SHP(2)*epsilon12,-0.5*SHP(2)*epsilon22;
         0,0,0,0.5*SHP(2)*epsilon11,0.5*SHP(2)*epsilon12;
         0,0,0,0,0];
Lamb3 = [zeros(3,5);
         0,0,0,-0.5*SHP(3)*epsilon12,-0.5*SHP(3)*epsilon22;
         0,0,0,0.5*SHP(3)*epsilon11,0.5*SHP(3)*epsilon12;
         0,0,0,0,0];
Lamb4 = [zeros(3,5);
         0,0,0,-0.5*SHP(4)*epsilon12,-0.5*SHP(4)*epsilon22;
         0,0,0,0.5*SHP(4)*epsilon11,0.5*SHP(4)*epsilon12;
         0,0,0,0,0];
matrix_Lambda = [Lamb1;Lamb2;Lamb3;Lamb4];
B = matrix_Lambda*N;
end