filepath_UT = [cd,'\NU\'];
filepath_U = [cd,'\file2_Ndsp\'];
strainTfile = [cd,'\LE\'];
% ===============
NodeLabel = 39;
OutPut_Ui = [1,2,3];
% ===============
totalframe = length(dir(strcat(strainTfile,'\','*.txt')))/2-1;

UT = zeros(totalframe,width(OutPut_Ui));
U = zeros(totalframe,width(OutPut_Ui));

for framenum = 1:totalframe
    file_UT = [filepath_UT,'NU_',num2str(framenum),'.txt'];
    file_U = [filepath_U,'Ndsp_',num2str(framenum),'.txt'];
    mat_UT = importdata(file_UT);
    mat_U = importdata(file_U);
    UT(framenum,:) = mat_UT(NodeLabel,OutPut_Ui);
    U(framenum,:) = mat_U(NodeLabel,OutPut_Ui);
end


plot(U(:,1),'*b')
hold on
plot(UT(:,1))