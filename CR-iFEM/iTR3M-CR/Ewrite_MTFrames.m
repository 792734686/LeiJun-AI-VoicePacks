filepath_U = [cd,'\file4_Ndsp\'];
strainTfile = [cd,'\LE\'];
% ===============
t = 2;
% ===============
totalframe = length(dir(strcat(strainTfile,'\','*.txt')))/2-1;
Nxy0 = importdata('Nxy0.txt');

for framenum = 1:totalframe
    file_U = [filepath_U,'NU_',num2str(framenum),'.txt'];
    Ndsp = importdata(file_U);
    Nxy = Nxy0 + Ndsp(:,1:3);
    [Eep_FEM,NEBT,NETP] = strain(Enod,Nxy,Ndsp,t);

    filename = [filepath_U,'NEBT_',num2str(framenum),'.txt'];
    writematrix(NEBT,filename,'delimiter','\t');

    filename = [filepath_U,'NETP_',num2str(framenum),'.txt'];
    writematrix(NETP,filename,'delimiter','\t');
end

function [Eep_FEM,NEBT,NETP] = strain(Enod,Nxy,Ndsp,t)
Nxy0 = importdata('Nxy0.txt');
Nxy0 = [(1:size(Nxy0,1))',Nxy0];
for i = 1:size(Enod,1)
    NN = Enod(i,2:end);
    ENCD = [Nxy(NN,1),Nxy(NN,2),Nxy(NN,3)];
    ENC0 = [Nxy0(NN,2),Nxy0(NN,3),Nxy0(NN,4)];
    [T,encd] = transmat_iTR3(ENCD);
    [T0,enc0] = transmat_iTR3(ENC0);
    [B,~] = geom_iTR3(enc0,0,0);
    ELM_u_LC = [(encd-enc0)';[0,0,0]];
    ELM_Omega_GC = [Ndsp(NN,4),Ndsp(NN,5),Ndsp(NN,6)]';
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
NEBT = Eep_FEM(:,1:3) - 0.5*t*Eep_FEM(:,4:6);
NETP = Eep_FEM(:,1:3) + 0.5*t*Eep_FEM(:,4:6);
end
