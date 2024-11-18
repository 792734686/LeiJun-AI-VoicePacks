function pltdfm2d4n(Enod,Nxy,Ndsp,enlarge)
% function to plot displacement contour figure
% Nxy(N,3) -- ; Enod(M,5) -- ; Evar(M,2); Nvar(N,2) 
figure
hold on
axis equal
axis off
m = size(Enod,1);
for i=1:m
    k = Enod(i,2:5);
    x = Nxy(k,2);
    y = Nxy(k,3);
    x = [x;x(1)];
    y = [y;y(1)];
    plot(x,y,'k')
end
for i=1:m
    k = Enod(i,2:5);
    x = Nxy(k,2)+enlarge*Ndsp(k,2);
    y = Nxy(k,3)+enlarge*Ndsp(k,3);
    x = [x;x(1)];
    y = [y;y(1)];
    plot(x,y,'r')
end
end

