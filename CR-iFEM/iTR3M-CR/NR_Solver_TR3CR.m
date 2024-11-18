function [Ndsp, iterations,converged] = NR_Solver_TR3CR(Nxy,ELEM,EELM,t,Eep,Eep2,tol,max_iter,NSTP,load_factor)
    x = zeros(6*size(Nxy,1),1);
    iterations = 0;
    Nxy0 = Nxy;
    Ndsp = zeros(size(Nxy,1),7);
    Theta = zeros(size(Nxy,1),3)';
    Rxy = zeros(size(Nxy,1),3)';
    bar = waitbar(0,'waiting for solving...');
    converged = true;
    
    while iterations < max_iter
        [Jx,F_in] = cr_globalKF_TR3(Nxy,Nxy0,ELEM,Eep,Eep2,EELM,t,Ndsp,Theta);
        L = zeros(size(Nxy,1)*6,1);
        L(6*(17-1)+5) = 26.167;
        L(6*(34-1)+5) = 26.167;
        fx = load_factor*L-F_in;
        %-------------------------------------------
        FIX = importdata('fix.txt');
        [Jx,fx] = boundary(FIX,[1,2,3,4,5,6],Jx,fx);
%         XSY = importdata('xsymm.rpt');
%         [Jx,fx] = boundary(XSY,[1,5,6],Jx,fx);
%         ZSY = importdata('zsymm.rpt');
%         [Jx,fx] = boundary(ZSY,[3,4,5],Jx,fx);
%         U2 = importdata('u2.rpt');
%         [Jx,fx] = boundary(U2,2,Jx,fx);
        %-------------------------------------------
        if rcond(Jx) <= 1e-10
            converged = false;
            break
        end
        Jx = sparse(Jx);
        delta_x =Jx \ fx;
        x = x + delta_x;
        Ndsp = [Nxy(:,1),Ulin_to_mat(x,6)];
        DRxy = Ulin_to_mat(delta_x,6)';
        DRxy = DRxy(4:6,:);% 结点转动增量（全局坐标，列）
        Rxy = RotMat_node(Rxy,DRxy); %结点转动向量（全局坐标，列）
        Theta = DeltaOmega_to_DeltaTheta(DRxy,Theta); % 结点转动向量（局部，列）
        Ndsp(:,5:7) = Rxy';
        x = Umat_to_lin(Ndsp(:,2:7));
        
        tolerance = tol;
        
        if norm(delta_x)/norm(x) < tolerance
            waitbar(1,bar,'Solve Succeed!');
            break;
        end
        
        iterations = iterations + 1;
        waitbar(iterations/max_iter,bar,['Solving...',' Step = ',num2str(NSTP),'Interations = ',num2str(iterations)]);
    end
    
    if iterations == max_iter
        converged = false;
    elseif converged == true
        disp(['Step- ',num2str(NSTP),' Solve Succeed.iters = ',num2str(iterations),', Load Factor = ',num2str(load_factor)]);
        converged = true;
    end
    close(bar);
end
