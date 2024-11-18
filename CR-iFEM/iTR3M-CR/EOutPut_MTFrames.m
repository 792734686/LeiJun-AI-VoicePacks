filepath_U = [cd,'\file4_Ndsp\'];
strainTfile = [cd,'\LE\'];
% ===============
ElementLabel = 39;
E = 21000;
mu = 0.3;
% ===============
D = E/(1-mu^2)*[1,mu,0;mu,1,0;0,0,0.5*(1-mu)];
totalframe = length(dir(strcat(strainTfile,'\','*.txt')))/2-1;

for framenum = 1:totalframe
    NEfile = [filepath_U,'NEBT_',num2str(framenum),'.txt'];
    NEBT = importdata(NEfile);
    SBT = zeros(size(NEBT));
    for enum = 1:size(SBT,1)
        SBT(enum,:) = (D*NEBT(enum,:)')';
    end
    filename = [filepath_U,'SBT_',num2str(framenum),'.txt'];
    writematrix(SBT,filename,'delimiter','\t');
end
