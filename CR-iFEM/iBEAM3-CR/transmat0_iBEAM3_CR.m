function [T0,len0] = transmat0_iBEAM3_CR(ENC0,ENC00,Rxy_iSTEP)
len0 = norm(ENC0(:,1)-ENC0(:,2));
Ee0 = transmat_iBEAM3(ENC00)';
E0 = Ee0(1:3,1:3);
Rxyi1 = Rxy_iSTEP(:,1);
Rxyi2 = Rxy_iSTEP(:,2);
    function R = funR(array)
        a = norm(array);
        if a == 0
            R = eye(3);
        else
            R = eye(3)+sin(a)/a*skew(array)+...
                0.5*(sin(0.5*a)/(0.5*a))^2*skew(array)^2;
        end
    end
b1 = funR(Rxyi1)*E0*[0,1,0]';
b2 = funR(Rxyi2)*E0*[0,1,0]';
b = 0.5*(b1+b2);
X1 = ENC0(:,1);
X2 = ENC0(:,2);
n1 = (X2-X1)/norm(X2-X1);
n3 = cross(n1,b)/norm(cross(n1,b));
n2 = cross(n3,n1);
T0 = [n1,n2,n3]';
end