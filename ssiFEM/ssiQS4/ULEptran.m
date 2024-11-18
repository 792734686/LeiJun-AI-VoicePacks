function Eep = ULEptran(t,NSTP)
GepBT = importdata('NEBT.rpt').data;
GepTP = importdata('NETP.rpt').data;
GepBT1 = importdata('NEBT1.rpt').data;
GepTP1 = importdata('NETP1.rpt').data;

GepBT = NSTP\GepBT;
GepTP = NSTP\GepTP;
GepBT1 = NSTP\GepBT1;
GepTP1 = NSTP\GepTP1;
n = size(GepBT,1);
% exx1 = GepTP(:,1);
% eyy1 = GepTP(:,2);
% gxy1 = GepTP(:,3);
% exx2 = GepBT(:,1);
% eyy2 = GepBT(:,2);
% gxy2 = GepBT(:,3);

exx1 = GepTP1(:,2);
eyy1 = GepTP1(:,1);
gxy1 = -GepTP1(:,3);
exx2 = GepBT(:,2);
eyy2 = GepBT(:,1);
gxy2 = -GepBT(:,3);

exx = 0.5*(exx1+exx2);
eyy = 0.5*(eyy1+eyy2);
gxy = 0.5*(gxy1+gxy2);
kxx = t\(exx1-exx2);
kyy = t\(eyy1-eyy2);
kxy = t\(gxy1-gxy2);
Eep = [exx,eyy,gxy,kxx,kyy,kxy,zeros(n,1),zeros(n,1)];
end