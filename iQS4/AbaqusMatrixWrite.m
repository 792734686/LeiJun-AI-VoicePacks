%%%%%%%%%%%%%%%

%% Write ABAQUS generated Matrix into Matlab

ABAQUS_Matrix_rawdata=load('Job-1_MTX_STIF2.mtx');

n = size(Enod,1);

GK_Abaqus=zeros(6*n,6*n);

for i=1:1:max(size(ABAQUS_Matrix_rawdata))

GK_Abaqus(ABAQUS_Matrix_rawdata(i,1),ABAQUS_Matrix_rawdata(i,2))=ABAQUS_Matrix_rawdata(i,3);

end

[GK_Abaqus,GF] = boundary(FIX,[1,2,3,4,5,6],GK_Abaqus,GF);

NF = GK_Abaqus*Umat_to_lin(Ndsp(:,2:7));
NF = Ulin_to_mat(NF,6);

%%%%%%%%%%%%%%%
%%

GK_FEM=FEM_gstiffm_iQS4(Nxy,Enod,t);

[GK_FEM,~] = boundary(FIX,[1,2,3,4,5,6],GK_FEM,GF);

Ndsplin = Ndsp(:,2:7);
NF = GK_FEM*Umat_to_lin(Ndsplin);
NF = Ulin_to_mat(NF,6);

%%
GF_FEM = zeros(size(Ndsp(:,2:7)));
GF_FEM(:,3) = -100;
GF_FEM = Umat_to_lin(GF_FEM);
[GK_FEM,GF_FEM] = boundary(FIX,[1,2,3,4,5,6],GK_FEM,GF_FEM);
ans = GK_FEM\GF_FEM;
ans = Ulin_to_mat(ans,6);
