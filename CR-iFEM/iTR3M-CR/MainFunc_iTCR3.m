function Ndsp = MainFunc_iTCR3(Enod,Nxy,t,EELM,Eep,Eep2)

    N_max = 1000; % 最大增量步数
    load_factor = zeros(N_max,1); % 载荷系数向量
    min_factor = 1e-4; % 最小载荷系数
    iter_max = 12; % 最大迭代次数

    % ------------------初始变量
    N = 1; % 初始增量步标记
    load_factor(1) = 0.2; % 初始载荷系数
    tol = 1e-4; % 容差
    converged = true;
    Ndsp = zeros(size(Nxy,1),7);
    Ndsp(:,1) = (1:size(Nxy))';% 初始结点位移
    % ------------------
    while N <= N_max && load_factor(N) >= min_factor
        [delta_Ndsp,iter,converged] =...
            NR_Solver_iTR3CR_MTFrames...
            (Nxy,Enod,EELM,t,Eep,Eep2,tol,iter_max,load_factor(N));

        while ~converged
            load_factor(N) = 0.5*load_factor(N);
            [delta_Ndsp,iter,converged] =...
                NR_Solver_iTR3CR_MTFrames...
                (Nxy,Enod,EELM,t,Eep,Eep2,tol,iter_max,load_factor(N));
        end
        Ndsp(:,2:4) = Ndsp(:,2:4)+delta_Ndsp(:,2:4);
        Ndsp(:,5:7) = RotMat_node(Ndsp(:,5:7)',delta_Ndsp(:,5:7)')';
        Nxy(:,2:4) = Nxy(:,2:4) + delta_Ndsp(:,2:4);

        if sum(load_factor) >= 1
            break
        end
        % ----------------------------------
        % 自动控制下一步载荷系数
        if iter <= 0.3*iter_max
            load_factor(N+1) = 1.5*load_factor(N);
        elseif iter <= 0.5*iter_max
            load_factor(N+1) = 1.2*load_factor(N);
        else
            load_factor(N+1) = load_factor(N);
        end
        if sum(load_factor) > 1
            load_factor(N+1) = 1-sum(load_factor(1:N));
        end
        % ----------------------------------
        N = N+1;
    end
    if load_factor(N) < min_factor
        disp(['Load Factor Too Small in Step - ',num2str(N),...
            '. LF = ',num2str(load_factor(N)),...
            ', Minimum = ',num2str(min_factor),'.'])
    end

end