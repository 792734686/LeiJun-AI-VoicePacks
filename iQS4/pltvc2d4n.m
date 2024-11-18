function pltvc2d4n(Enod,Nxy,Nvar)
% function to plot field variable contour figure for element with 4 nodes
figure
hold on
axis equal
axis off
m = size(Enod,1);
for i=1:m
    k = Enod(i,2:5);
    x = Nxy(k,2);
    y = Nxy(k,3);
    c = Nvar(k,2);  % 
    fill(x,y,c)
end
colorbar('location','easto utside')
end

