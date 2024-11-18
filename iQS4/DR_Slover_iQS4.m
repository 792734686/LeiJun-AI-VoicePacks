function [Ndsp, iterations] = DR_Slover_iQS4(Nxy,ELEM,EELM,t,Eep,Eep2, tol, max_iter)
    x = zeros(6*size(Nxy,1),1);
    iterations = 0;
    Ndsp = zeros(size(Nxy,1),7);
%     Ndsp = importdata('NDSP--400.txt');
    [K0,~] = gstiffm_LargeDeform_iQS4(Nxy,ELEM,EELM,t,Eep,Eep2,Ndsp);
    fx = zeros(size(x));
    FIX = importdata('fix.rpt');
    HOLD23 = importdata('hold23.rpt');
    HOLD235 = importdata('hold235.rpt');
    [K0,~] = boundary(FIX,[1,2,3,4,5,6],K0,fx);
    [K0,~] = boundary(HOLD23,[2,3],K0,fx);
    [K0,~] = boundary(HOLD235,[2,3,5],K0,fx);
    
    while iterations < max_iter
        Jx = gstiffm_LargeDeform_iQS4(Nxy,ELEM,EELM,t,Eep,Eep2,Ndsp);
        
        fx = gforcem_iQS4(Nxy,ELEM,Eep,Eep2,EELM,t);
        %-------------------------------------------
        [Jx,fx] = boundary(FIX,[1,2,3,4,5,6],Jx,fx);
        [Jx,fx] = boundary(HOLD23,[2,3],Jx,fx);
        [Jx,fx] = boundary(HOLD235,[2,3,5],Jx,fx);
        %-------------------------------------------
        fx = fx-Jx*x;
        
        K0 = sparse(K0);
        delta_x =K0 \ fx;
        x = x + delta_x
        Ndsp = [Nxy(:,1),Ulin_to_mat(x,6)];
        
        tolerance = tol*max(abs(x));
        if norm(delta_x) < tolerance
            break;
        end
        
        iterations = iterations + 1
    end
    
    if iterations == max_iter
        disp('Too many iterations for Solver.');
    else
        disp('Solve Succeed.');
    end
end
