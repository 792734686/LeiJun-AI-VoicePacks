% NSTP = 16;
% for i = 1:NSTP
%     
%     filename = [cd,'\file_Ndsp\NU_',num2str(i),'.txt'];
%     Ndsp = importdata(filename);
%     [Eep_FEM,NEBT,NETP] = strain(Enod,Nxy,Ndsp,t);
%     writematrix(NEBT,['NEBT_',num2str(i),'.txt'])
%     writematrix(NETP,['NETP_',num2str(i),'.txt'])
%     
% end

[Eep_FEM,NEBT,NETP] = strain(Enod,Nxy,Ndsp,t);
writematrix(NEBT,'NEBT.txt')
writematrix(NETP,'NETP.txt')

function [Eep_FEM,NEBT,NETP] = strain(Enod,Nxy,Ndsp,t)
Nxy0 = importdata('Nxy0.txt');
Nxy0 = [(1:size(Nxy0,1))',Nxy0];
for i = 1:size(Enod,1)
    NN = Enod(i,2:end);
    ENCD = [Nxy(NN,2),Nxy(NN,3),Nxy(NN,4)];
    ENC0 = [Nxy0(NN,2),Nxy0(NN,3),Nxy0(NN,4)];
    [T,encd] = transmat_iTR3(ENCD);
    [T0,enc0] = transmat_iTR3(ENC0);
    [B,~] = geom_iTR3(enc0,0,0 );
    ELM_u_LC = [(encd-enc0)';[0,0,0]];
    ELM_Omega_GC = [Ndsp(NN,5),Ndsp(NN,6),Ndsp(NN,7)]';
    for j = 1:3
        NOD_Omega_GC = ELM_Omega_GC(:,j);
        NOD_R_LC = expm(skew(NOD_Omega_GC));
        NOD_Omega_LC = Array_by_Quaternions(T(1:3,1:3)*NOD_R_LC*T0(1:3,1:3)');
        ELM_Omega_LC(:,j) = NOD_Omega_LC;
    end
    ELM_P_LC = Umat_to_lin([ELM_u_LC;ELM_Omega_LC]');
    ELM_Strain = B*ELM_P_LC;
    Eep_FEM(i,:) = ELM_Strain';
end
NEBT = [Enod(:,1),Eep_FEM(:,1:3) - 0.5*t*Eep_FEM(:,4:6)];
NETP = [Enod(:,1),Eep_FEM(:,1:3) + 0.5*t*Eep_FEM(:,4:6)];
end


