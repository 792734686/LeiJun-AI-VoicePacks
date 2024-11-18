function [Ndsp, iterations] = NR_Solver_iQS4(Nxy,ELEM,EELM,t,Eep,Eep2, tol, max_iter,NSTP)
    x = zeros(6*size(Nxy,1),1);
    iterations = 0;
    Ndsp = zeros(size(Nxy,1),7);
    bar = waitbar(0,'waiting for solving...');
%     Ndsp = importdata('NDSP0.txt');
    
    while iterations < max_iter
        [Jx,~] = gstiffm_LargeDeform_iQS4(Nxy,ELEM,EELM,t,Eep,Eep2,Ndsp);
        
        fx = gforcem_iQS4(Nxy,ELEM,Eep,Eep2,EELM,t);
        %-------------------------------------------
        FIX = importdata('fix.rpt');
        [Jx,fx] = boundary(FIX,[1,2,3,4,5,6],Jx,fx);
        %-------------------------------------------
        fx = fx-Jx*x;
        
        Jx = sparse(Jx);
        delta_x =Jx \ fx;
        x = x + delta_x;
        Ndsp = [Nxy(:,1),Ulin_to_mat(x,6)];
        
        tolerance = tol*max(abs(x));
        
        if norm(delta_x) < tolerance
            waitbar(1,bar,'Solve Succeed!');
            break;
        end
        
        iterations = iterations + 1;
        waitbar(iterations/max_iter,bar,['Solving...',' Step = ',num2str(NSTP),'Interations = ',num2str(iterations)]);
    end
    
    if iterations == max_iter
        disp(['Too many iterations for ','Step- ',num2str(NSTP)]);
    else
        disp(['Step- ',num2str(NSTP),' Solve Succeed.']);
    end
    close(bar);
end
