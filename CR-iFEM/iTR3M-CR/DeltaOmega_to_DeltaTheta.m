function Theta = DeltaOmega_to_DeltaTheta(DRxy,Theta0)
DOmega = DRxy;
DTheta = zeros(size(Theta0));
for i = 1:size(DRxy,2)
    Theta0_i = Theta0(:,i);
    DOmega_i = DOmega(:,i);
    if norm(Theta0_i) == 0
        H = eye(3);
    else
        eta = (1-0.5*norm(Theta0_i)*cot(0.5*norm(Theta0_i)))/(norm(Theta0_i)^2);
        H = eye(3)-0.5*skew(Theta0_i)+eta*skew(Theta0_i)^2;
    end
    DTheta_i = H*DOmega_i;
    DTheta(:,i) = DTheta_i;
end
Theta = Theta0+DTheta;
end