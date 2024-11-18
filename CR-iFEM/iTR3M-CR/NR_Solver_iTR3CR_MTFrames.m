function [Ndsp, iterations,converged] = NR_Solver_iTR3CR_MTFrames(Nxy,ELEM,EELM,t,Eep,Eep2,tol,max_iter,load_factor)
    x = zeros(6*size(Nxy,1),1);
    iterations = 0;
    Nxy0 = Nxy;
    Ndsp = zeros(size(Nxy,1),7);
    Rxy = zeros(size(Nxy,1),3)';
    converged = true;
    Eep = load_factor*Eep;
    Eep2 = load_factor*Eep2;
    
    while iterations < max_iter
        
        [Jx,fx] = cr_globalKF_iTR3(Nxy,Nxy0,ELEM,Eep,Eep2,EELM,t,Ndsp);
        %-------------- 边界条件 -------------------
        FIX = importdata('fix.txt');
%         BOUND = importdata('bound.txt');
%         XSYMM = importdata('xsymm.txt');
%         YSYMM = importdata('ysymm.txt');
%         ZSYMM = importdata('zsymm.txt');
        %         U2 = importdata('u2.txt');
        [Jx,fx] = boundary(FIX,[1,2,3,4,5,6],Jx,fx);
%         [Jx,fx] = boundary(BOUND,[1,2,4,5],Jx,fx);
%         [Jx,fx] = boundary(XSYMM,[1,5,6],Jx,fx);
%         [Jx,fx] = boundary(YSYMM,[2,4,6],Jx,fx);
%         [Jx,fx] = boundary(ZSYMM,[3,4,5],Jx,fx);
        %         [Jx,fx] = boundary(U2,2,Jx,fx);
        %-------------------------------------------
        if rcond(Jx) <= 1e-16
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
        Ndsp(:,5:7) = Rxy';
        x = Umat_to_lin(Ndsp(:,2:7));
        
        tolerance = tol;
        
        if norm(delta_x)/norm(x) < tolerance
            break;
        end
        
        iterations = iterations + 1;
    end
    
    if iterations == max_iter
        converged = false;
    elseif converged == true
        converged = true;
    end
end
