function [KP,KF] = boundary(Nnod,Dof,GK,GF)
%Nnod --- 施加边界条件的节点号
%Dof --- 施加边界条件的自由度
%-----------------------------------
for i = 1:size(Nnod)
    ii = Nnod(i);
    for j = 1:size(Dof,2)
        jj = Dof(j);
        dof = 6*ii-(6-jj);
        GK(:,dof) = 0;
        GK(dof,:) = 0;
        GK(dof,dof) = 1;
        GF(dof,1) = 0;
    end
end
KP = GK;
KF = GF;
end
