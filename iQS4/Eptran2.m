function Eep = Eptran2(t)
    GepBT = importdata('NEBTnurbs.txt');
    GepTP = importdata('NETPnurbs.txt');
%     GepBT(:,1) = [];
%     GepTP(:,1) = [];
%     GepBT(:,3) = [];
%     GepTP(:,3) = [];
    %-------------------------------
    GEBT = zeros(size(GepBT));
    GETP = zeros(size(GepTP));
    EPLOC = importdata('NELOC.txt');
    for i = 1:size(EPLOC,1)
        e = GepBT(i,:)';
        E = EPLOC(i,2:4)';
        if E(3)/e(3) > 0
            TK = [1,0,0;
                  0,1,0;
                  0,0,1];
        else
            TK = [0,1,0;
                  1,0,0;
                  0,0,-1];
        end
%         TK = E*pinv(e);
        EBTi = E';
        ETPi = TK*GepTP(i,:)';
        GEBT(i,:) = EBTi;
        GETP(i,:) = ETPi;
    end
    %--------------------------------
    n = size(GEBT,1);
    exx1 = GETP(:,1);
    eyy1 = GETP(:,2);
    gxy1 = GETP(:,3);
    exx2 = GEBT(:,1);
    eyy2 = GEBT(:,2);
    gxy2 = GEBT(:,3);
    exx = 0.5*(exx1+exx2);
    eyy = 0.5*(eyy1+eyy2);
    gxy = 0.5*(gxy1+gxy2);
    kxx = t\(exx1-exx2);
    kyy = t\(eyy1-eyy2);
    kxy = t\(gxy1-gxy2);
    Eep = [exx,eyy,gxy,kxx,kyy,kxy,zeros(n,1),zeros(n,1)];
end