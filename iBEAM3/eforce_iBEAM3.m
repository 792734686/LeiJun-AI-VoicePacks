function EF = eforce_iBEAM3(w,e,ENC,ngs)
% e ---- element strain(centroid),1*6
EF=zeros(14,1);
[GSP,GSW] = gauspw(ngs);
for i = 1:ngs
    s = GSP(i);
    wi = GSW(i);
    B = geom_iBEAM3(s,ENC);
    B1 = B(1,:); B2 = B(2,:); B3 = B(3,:);
    B4 = B(4,:); B5 = B(5,:); B6 = B(6,:);
    EF1 = B1'*e(1); EF2 = B2'*e(2); EF3 = B3'*e(3);
    EF4 = B4'*e(4); EF5 = B5'*e(5); EF6 = B6'*e(6);
    EF = EF + ...
    (w(1)*EF1 + w(2)*EF2 + w(3)*EF3 + w(4)*EF4 + w(5)*EF5 + w(6)*EF6)*wi;
end

end