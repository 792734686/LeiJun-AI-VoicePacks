NUmag = sqrt(Ndsp(:,2).^2+Ndsp(:,3).^2+Ndsp(:,4).^2);
NU = importdata('NU.txt');
NUT = [(1:size(NU,1))',NU];
NUmagT = sqrt(NUT(:,2).^2+NUT(:,3).^2+NUT(:,4).^2);
Nxy0 = Nxy;
figure
plot(NUmag)
hold on
plot(NUmagT)
Error = mean(abs(NUmagT-NUmag)./(max(abs(NUmagT))))*100
% --------------------------------------------------------
dlmwrite('Ndsp.txt', [Ndsp(:,2),Ndsp(:,3),Ndsp(:,4)], 'delimiter', '\t');
% dlmwrite('NDSP0.txt', NUT(:,2:4), 'delimiter', '\t');
% U1_err = [relative_deviation(Ndsp(:,2),NUT(:,2)),zeros(size(NUT,1),1),zeros(size(NUT,1),1)];
% dlmwrite('ERRU1.txt',U1_err,' ');
% U2_err = [relative_deviation(Ndsp(:,3),NUT(:,3)),zeros(size(NUT,1),1),zeros(size(NUT,1),1)];
% dlmwrite('ERRU2.txt',U2_err,' ');
% U3_err = [relative_deviation(Ndsp(:,4),NUT(:,4)),zeros(size(NUT,1),1),zeros(size(NUT,1),1)];
% dlmwrite('ERRU3.txt',U3_err,' ');
% Umag_err = [relative_deviation(NUmag,NUmagT),zeros(size(NUT,1),1),zeros(size(NUT,1),1)];
% dlmwrite('ERRUmag.txt',Umag_err,' ');
% max(Umag_err)

% dlmwrite('Enod.txt',Enod(:,2:5),' ');
% dlmwrite('Ncoordinate.txt',NU(:,5:7),' ');

figure
Ncoord = importdata("Nxy0.txt") + NU;
scatter3(Nxy(:,2),Nxy(:,3),Nxy(:,4))
hold on
scatter3(Ncoord(:,1),Ncoord(:,2),Ncoord(:,3))