function Ndspwrite(Ndsp,N)
working_disp = cd;
filename = [working_disp,'\file_Ndsp\','Ndsp_',num2str(N),'.txt'];
writematrix(Ndsp(:,2:4),filename,'delimiter','\t');
end