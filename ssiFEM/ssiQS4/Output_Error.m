NUmag = sqrt(Ndsp(:,2).^2+Ndsp(:,3).^2+Ndsp(:,4).^2);
NU = importdata('NU.txt');
NUmagT = sqrt(NU(:,1).^2+NU(:,2).^2+NU(:,3).^2);
plot(NUmag)
hold on
plot(NUmagT)
Error = mean(abs(NUmagT-NUmag)/max(abs(NUmagT)))*100
dlmwrite('Ndsp.txt', [Ndsp(:,2),Ndsp(:,3),Ndsp(:,4)], 'delimiter', '\t');
U1_err = [relative_deviation(Ndsp(:,2),NUT(:,2)),zeros(size(NUT,1),1),zeros(size(NUT,1),1)];
dlmwrite('ERRU1.txt',U1_err,' ');
U2_err = [relative_deviation(Ndsp(:,3),NUT(:,3)),zeros(size(NUT,1),1),zeros(size(NUT,1),1)];
dlmwrite('ERRU2.txt',U2_err,' ');
U3_err = [relative_deviation(Ndsp(:,4),NUT(:,4)),zeros(size(NUT,1),1),zeros(size(NUT,1),1)];
dlmwrite('ERRU3.txt',U3_err,' ');
Umag_err = [relative_deviation(NUmag,NUmagT),zeros(size(NUT,1),1),zeros(size(NUT,1),1)];
dlmwrite('ERRUmag.txt',Umag_err,' ');