function [Nxy,Enod] = mesh2d4n(P4,xen,yen)
% function to mesh 2-D plane to 4-node elements
% P4(4,2) -- coordinate matrix of corner points
% xen -- element number in direction x
% yen -- element number in direxction y 
xn = xen + 1;
yn = yen + 1;
XI = linspace(-1,1,xn);
ETA = linspace(-1,1,yn);
[X,Y] = meshgrid(XI,ETA);
%-----coordinate of node-----
for i=1:yn
    for j=1:xn
        k = (i-1)*xn+j;
        NXE(k,1)=X(i,j);
        NXE(k,2)=Y(i,j);
    end
end
n = xn*yn;
Nxy(:,1) = 1:n;
for i = 1:n
    xi = NXE(i,1);
    eta = NXE(i,2);
    n1 = (1-xi)*(1-eta)/4.0;
    n2 = (1+xi)*(1-eta)/4.0;
    n3 = (1+xi)*(1+eta)/4.0;
    n4 = (1-xi)*(1+eta)/4.0;
    Nxy(i,2) = n1*P4(1,1) + n2*P4(2,1) + n3*P4(3,1) + n4*P4(4,1);
    Nxy(i,3) = n1*P4(1,2) + n2*P4(2,2) + n3*P4(3,2) + n4*P4(4,2);
end
%-----nodes in element-----
for i=1:yn-1
    for j=1:xn-1
        k=(i-1)*(xn-1)+j;
        Enod(k,:)=[k,(i-1)*xn+j,(i-1)*xn+j+1,i*xn+j+1,i*xn+j];
    end
end
%----plot mesh and label node and element ----
figure
axis equal
axis off
hold on
k = size(Enod,1);
for i=1:k
    II = Enod(i,2:5);
    fill(Nxy(II,2),Nxy(II,3),'w')   %绘制单元
    hold on
    xc = mean(Nxy(II,2));
    yc = mean(Nxy(II,3));
    h = text(xc,yc,num2str(i));   %标注单元编号
    set(h,'color','b');
set(h,'fontsize',12)
end
k = size(Nxy,1);
for i=1:k
    h=text(Nxy(i,2),Nxy(i,3),num2str(i));   %标注结点号
    set(h,'color','r')
set(h,'fontsize',12)
end
end
