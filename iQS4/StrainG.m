function [G13,G23] = StrainG(Enod, Ndsp, Nxy) 
    NELM = size(Enod,1);
    G13 = zeros(NELM,1);
    G23 = zeros(NELM,1);
    
    for i = 1:NELM
        NODi = Enod(i,2:5);
        NUi = [];
        for nn = NODi
            NUi = [NUi;Ndsp(nn,2:7)];
        end
        NN = Enod(i,2:end);     %返回节点编号
        ENC = [Nxy(NN,2),Nxy(NN,3),Nxy(NN,4)];    %4*2，返回单元的节点坐标（行），4列表示4节点
        [g13,g23] = eStrainG(ENC,NUi);
        G13(i) = g13;
        G23(i) = g23;
        
    end
    
end