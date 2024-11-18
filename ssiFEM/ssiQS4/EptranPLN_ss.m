function Eep = EptranPLN_ss(filename)
    GepTP = importdata(filename);
    GepTP(:,1) = [];
    n = size(GepTP,1);
    exx = GepTP(:,1);
    eyy = GepTP(:,2);
    gxy = GepTP(:,3);
    Eep = [exx,eyy,gxy,zeros(n,1),zeros(n,1)];
end