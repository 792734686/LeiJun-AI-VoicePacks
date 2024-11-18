function NUwrite(Ndsp,N)
working_disp = cd;
filename = [working_disp,'\file_Ndsp\','NU_',num2str(N),'.txt'];
writematrix(Ndsp,filename,'delimiter','\t');
end