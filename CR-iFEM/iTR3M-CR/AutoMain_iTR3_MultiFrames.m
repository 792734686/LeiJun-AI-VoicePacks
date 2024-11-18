%%
%
% ---------------输出文件目录------------------
if exist ([cd,'\file_Ndsp\'],'dir') == 7
    rmdir('file_Ndsp','s')
end
mkdir('file_Ndsp')
% ---------------应变数据目录------------------
strainTfile = [cd,'\LE\'];
TotalFrame = length(dir(strcat(strainTfile,'\','*.txt')))/2-1;

% ============== 创建循环 ================

for framenum = 1:TotalFrame

    % ---------------------------------------------
    %%
    clear Ndsp delta_Ndsp;
    %输入矩阵Nxy,Enod:
    %Enod:[单元号|节点1|节点2|节点3|节点4]
    %Nxy:[节点号|x坐标|y坐标|z坐标]
    %---------------单元信息---------------------
    Enod = importdata('Enod.txt');
    Enod = [(1:size(Enod,1))',Enod];
    Nxy = importdata('Nxy0.txt');
    Nxy = [(1:size(Nxy,1))',Nxy];

    % -----------板厚、载荷与测量单元-------------
    t = 1;      %板厚
    % EELM = importdata('EELM.rpt');
    EELM = Enod(:,1);   %输入测量应变的单元
    Eep = EptranPLN_MTFrames(strainTfile,framenum,t);    %转换为8应变分量Eep(单元数1,8)
    Eep2 = zeros(size(Eep));

    %----------------------------------------------
    %% 自动增量——NR求解器
    N_max = 1000; % 最大增量步数
    load_factor = zeros(N_max,1); % 载荷系数向量
    min_factor = 1e-4; % 最小载荷系数
    iter_max = 20; % 最大迭代次数

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
    Ndspwrite(Ndsp,framenum)
    NUwrite(Ndsp,framenum)

    disp(['Frame_',num2str(framenum),' Solve Succeed!'])

end