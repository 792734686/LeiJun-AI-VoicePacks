%%
clear;
%输入矩阵Nxy,Enod:
%Enod:[单元号|节点1|节点2|节点3|节点4]
%Nxy:[节点号|x坐标|y坐标|z坐标]
%---------------单元信息---------------------
Enod = importdata('Enod.txt');
Enod = [(1:size(Enod,1))',Enod];
Nxy = importdata('Nxy0.txt');
Nxy = [(1:size(Nxy,1))',Nxy];

% -----------载荷与测量单元-------------
% EELM = importdata('EELM.txt');
EELM = Enod(:,1);   %输入测量应变的单元
Eep = beamstrain([1,1,1,0,0,1]);    
Eep2 = zeros(size(Eep));

% ---------------输出文件目录------------------
if exist ([cd,'\file_Ndsp\'],'dir') == 7
    rmdir('file_Ndsp','s')
end
mkdir('file_Ndsp')
%----------------------------------------------
%% 自动增量——NR求解器
N_max = 100; % 最大增量步数
load_factor = zeros(N_max,1); % 载荷系数向量
min_factor = 1e-4; % 最小载荷系数
iter_max = 20; % 最大迭代次数
time = [];

% ------------------初始变量
N = 1; % 初始增量步标记
load_factor(1) = 0.1; % 初始载荷系数
tol = 1e-4; % 容差
converged = true;
Ndsp = zeros(size(Nxy,1),7);
Ndsp(:,1) = (1:size(Nxy))';
NOD_ROT = [zeros(size(Ndsp,1),3)];% 初始结点位移
tic;
% ------------------
while N <= N_max && load_factor(N) >= min_factor
    % -------------------
    [delta_Ndsp,iter,converged] =...
        NR_Solver_iBEAM3CR(Nxy,NOD_ROT,Enod,EELM,Eep,Eep2,tol,iter_max,N,load_factor(N));
    
    while ~converged
        load_factor(N) = 0.5*load_factor(N);
        [delta_Ndsp,iter,converged] =...
            NR_Solver_iBEAM3CR(Nxy,NOD_ROT,Enod,EELM,Eep,Eep2,tol,iter_max,N,load_factor(N));
    end
    Ndsp(:,2:4) = Ndsp(:,2:4)+delta_Ndsp(:,2:4);
    Ndsp(:,5:7) = RotMat_node(Ndsp(:,5:7)',delta_Ndsp(:,5:7)')';
    NOD_ROT = Ndsp(:,5:7);
    Nxy(:,2:4) = Nxy(:,2:4) + delta_Ndsp(:,2:4);
    time(N,:) = [N,sum(load_factor),Ndsp(17,2),Ndsp(17,4)];
    Ndspwrite(Ndsp,N)
    NUwrite(Ndsp,N)
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
toc;

function SE = beamstrain(switcher)
full_strain = importdata('SE.txt');
SE = zeros(size(full_strain));
for i = 1:size(switcher,2)
    SE(:,i) = switcher(i)*full_strain(:,i);
end
end
