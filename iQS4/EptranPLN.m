function Eep = EptranPLN(t)
    GepBT = importdata('NEBT.txt');
    GepTP = importdata('NETP.txt');
    GepBT(:,1) = [];
    GepTP(:,1) = [];
    
   %===========================
%    g13 = importdata('G13.rpt').data;
%    g13 = g13/80769;
%    g23 = importdata('G23.rpt').data;
%    g23 = g23/80769;
   %===========================
%     GepBT(:,1) = [];
%     GepTP(:,1) = [];
%     GepBT(:,3) = [];
%     GepTP(:,3) = [];
    n = size(GepBT,1);
    exx1 = GepTP(:,1);
    eyy1 = GepTP(:,2);
    gxy1 = GepTP(:,3);
    exx2 = GepBT(:,1);
    eyy2 = GepBT(:,2);
    gxy2 = GepBT(:,3);
    
    exx = 0.5*(exx1+exx2);
    eyy = 0.5*(eyy1+eyy2);
    gxy = 0.5*(gxy1+gxy2);
    kxx = t\(exx1-exx2);
    kyy = t\(eyy1-eyy2);
    kxy = t\(gxy1-gxy2);
    Eep = [exx,eyy,gxy,kxx,kyy,kxy,zeros(n,1),zeros(n,1)];
%     Eep = [exx,eyy,gxy,kxx,kyy,kxy,g13,g23];
end