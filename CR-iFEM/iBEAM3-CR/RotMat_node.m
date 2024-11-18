function Rxy_new = RotMat_node(Rxy,DRxy)
    function R = funR(array)
        a = norm(array);
        if a == 0
            R = eye(3);
        else
            R = eye(3)+sin(a)/a*skew(array)+...
                0.5*(sin(0.5*a)/(0.5*a))^2*skew(array)^2;
        end
    end
Rxy_new = zeros(size(Rxy));
for i = 1:size(DRxy,2)
    Rxy_i = Rxy(:,i);
    DRxy_i = DRxy(:,i);
    R_oldi = funR(Rxy_i);
    R = funR(DRxy_i);
    R_newi = R*R_oldi;
    Rxy_newi = Array_by_Quaternions(R_newi);
    Rxy_new(:,i) = Rxy_newi;
end
end