function [NDSP, iterations] = NR_Solver_iQS4(Nxy,ELEM,EELM,t,Eep,Eep2, tol, max_iter,NSTP)
    Eep = Eep/NSTP;
    Eep2 = Eep2/NSTP;
    Nxy0 = Nxy;
    NDSP = zeros(size(Nxy,1),7);
    NDSP(:,1) = Nxy(:,1);
    Uxy = zeros(size(Nxy,1),3);
    
    for i = 1:NSTP
    x = zeros(6*size(Nxy,1),1);
    Ndsp = zeros(size(Nxy,1),7);
    iterations = 0;
    bar = waitbar(0,'waiting for solving...');
%     Ndsp = importdata('NDSP0.txt');

        
    while iterations < max_iter
        [Jx,~] = gstiffm_LargeDeform_iQS4(Nxy+[zeros(size(Nxy,1),1),Uxy],ELEM,EELM,t,Eep,Eep2,Ndsp);
        
        fx = gforcem_iQS4(Nxy+[zeros(size(Nxy,1),1),Uxy],ELEM,Eep,Eep2,EELM,t);
        %-------------------------------------------
        FIX = importdata('fix.rpt');
        [Jx,fx] = boundary(FIX,[1,2,3,4,5,6],Jx,fx);
        %-------------------------------------------
        fx = fx-Jx*x;
        
        Jx = sparse(Jx);
        delta_x =Jx \ fx;
        
        Ndsp = Ulin_to_mat(x,6);
        DNdsp = Ulin_to_mat(delta_x,6);
        Uxy = Ndsp(:,1:3);
        Rxy = Ndsp(:,4:6);
        DUxy = DNdsp(:,1:3);
        DRxy = DNdsp(:,4:6);
        Uxy = Uxy+DUxy;
        Rxy = RotMat_node(Rxy',DRxy')';
        x = Umat_to_lin([Uxy,Rxy]);
        Ndsp = [Nxy(:,1),Ulin_to_mat(x,6)];
        
        tolerance = tol*max(abs(x));
        
        if norm(delta_x) < tolerance
            waitbar(1,bar,'Solve Succeed!');
            break;
        end
        
        iterations = iterations + 1;
        waitbar(iterations/max_iter,bar,['Solving...',' Step = ',num2str(i),'Interations = ',num2str(iterations)]);
    end
    
    NDSP(:,2:4) = NDSP(:,2:4) + Ndsp(:,2:4);
    NDSP(:,5:7) = RotMat_node(NDSP(:,5:7)',Ndsp(:,5:7)')';
    Nxy(:,2:4) = Nxy0(:,2:4)+NDSP(:,2:4);
    
    if iterations == max_iter
        disp(['Too many iterations for ','Step- ',num2str(i)]);
    else
        disp(['Step- ',num2str(i),' Solve Succeed.']);
    end
    close(bar);
    end
end
