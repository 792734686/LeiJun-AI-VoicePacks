function EK = estiff_iBEAM3(w,ENC,ngs)
EK=zeros(14,14);
[GSP,GSW] = gauspw(ngs);
% [Te,enc] = transmat_iQS4(ENC);
for i = 1:ngs
    s = GSP(i);
    wi = GSW(i);
    B = geom_iBEAM3(s,ENC);
    B1 = B(1,:); B2 = B(2,:); B3 = B(3,:);
    B4 = B(4,:); B5 = B(5,:); B6 = B(6,:);
    EK1 = B1'*B1; EK2 = B2'*B2; EK3 = B3'*B3;
    EK4 = B4'*B4; EK5 = B5'*B5; EK6 = B6'*B6;
    EK = EK + ...
    (w(1)*EK1 + w(2)*EK2 + w(3)*EK3 + w(4)*EK4 + w(5)*EK5 + w(6)*EK6)*wi;
end

end