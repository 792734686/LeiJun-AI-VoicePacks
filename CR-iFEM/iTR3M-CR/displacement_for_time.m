time_step = 20;
target_node = 17;
displacement_u = zeros(time_step,3);
for i = 1:time_step
    filename = [cd , '\file2_Ndsp\' , 'Ndsp_' , num2str(i) , '.txt'];
    displacement = importdata(filename);
    displacement_u(i,:) = displacement(target_node,:);
end