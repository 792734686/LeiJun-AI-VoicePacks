function Eep = EptranPLN_MTFrames(strainTfile,framenum,t)
    NEBTfilename = [strainTfile,'NEBT_',num2str(framenum),'.txt'];
    NETPfilename = [strainTfile,'NETP_',num2str(framenum),'.txt'];

    GepBT = importdata(NEBTfilename);
    GepTP = importdata(NETPfilename);
    GepBT(:,1) = [];
    GepTP(:,1) = [];
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
end