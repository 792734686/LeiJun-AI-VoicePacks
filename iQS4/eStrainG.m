function [g13,g23] = eStrainG(ENC,NUi);
    [Te,enc] = transmat_iQS4(ENC);
    ui = Te*Umat_to_lin(NUi);
    xi = 0;
    eta = 0;
    
    [B,~] = geom_iQS4( enc,xi,eta );
    Bs = B(7:8,:);
    ElementG = Bs*ui;
    g13 = ElementG(1);
    g23 = ElementG(2);
end