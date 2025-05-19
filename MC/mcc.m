clc
clear
% 关闭警告
% warning off
% 定义参数
policy=0.1;
N = 2000;
stp = 30;
lamda = 0.3;
delta = 0.6;
csmd_G = 0.01;
csmd_AS = 0.01;
gama = 1.5;
mu = 0.8;
MC_rep = 20;   % 仿真次数
termi = 40;
BETA_AS_AVER = zeros(1, termi);
BETA_US_AVER = zeros(1, termi);
BETA_AG_AVER = zeros(1, termi);
% 每次仿真的结果预分配
BETA_AG = zeros(MC_rep, termi);
BETA_US = zeros(MC_rep, termi);
BETA_AS = zeros(MC_rep, termi);
% 引入邻接矩阵并转换为逻辑类型（A和B是从'ba'和'ER'加载的）
% 将邻接矩阵转成逻辑类型的作用是将邻接矩阵 A 和 B 转换为逻辑类型（布尔类型），也就是将矩阵中的所有非零元素转换为 true，零元素转换为 false。
%在后续代码中只是进行存在性检查或逻辑运算，使用逻辑类型可以避免 MATLAB 将数据看作数值进行多余的计算，比如 A(i,j) == 1 可以被简化为 A(i,j)，直接判断布尔值，而无需和数字比较。
load ba 
load ER 
A = logical(A);
B = logical(B);
% 预计算每个节点的邻居，减少循环中的查找时间
%这段代码的作用是为邻接矩阵 A 和 B 中的每一个节点找到它们的邻居，并将这些邻居存储在 A_neighbors 和 B_neighbors 两个 cell 数组中。
A_neighbors = cell(1, N);
for i = 1:N
    A_neighbors{i} = find(A(:, i));
end
B_neighbors = cell(1, N);
for i = 1:N
    B_neighbors{i} = find(B(:, i));
end
% 主循环
parfor rep = 1:MC_rep  % 仿真次数
    for l = 1:termi  % β变化
         disp(['第', num2str(rep), '次仿真，β为', num2str(l)])
        beta_R = l / 80;
        
        % 向量化初始化x
        x = rand(1, N) < csmd_G;%表示将每个随机数与阈值 csmd_G 进行比较，如果随机数小于 csmd_G，则结果为 true (1)，否则为 false (0)。
      %这样，x 是一个大小为 1 x N 的逻辑数组，其中每个元素是 true 或 false，表示是否满足小于 csmd_G 的条件。这个操作相当于以 csmd_G 为概率初始化了 x 数组中的值，x(i)=1 表示这个节点属于状态 G，x(i)=0 表示不属于状态 G。  
        % 向量化初始化m
        m = rand(1, N) < (2 * csmd_AS);
        m(x) = 1;  % 如果x(i)==1, 则m(i)=1，x 是一个逻辑数组（即 x(i)==1 的节点），这个操作对满足 x(i)==1 的节点对应的 m(i) 进行赋值，即将这些节点对应的 m(i) 强制设为 1
        x(~m) = 0; % 如果m(i)==0, 则x(i)=0，当逻辑数组（如 m）被用于索引时，m 中为 1 或 true 的元素会被选中，所以x(m)=0是不对的
        % 仿真步长
        for t = 1:stp
            % 上层传播
            % 处理m==0的节点，U节点
            m_zero_idx = find(m == 0);%查找U节点所在的列
            n = m;  % 先复制m的状态
            for i = m_zero_idx
                neighbors = A_neighbors{i};%查找U节点的邻居节点所在的列
                active_neighbors = neighbors(m(neighbors) == 1);%判断邻居节点是否为A节点，如果事就可以进行感染
                if ~isempty(active_neighbors)%检查活跃节点是否为空，如果不为空就继续下面的循环
                    p1 = rand(1, length(active_neighbors));
                    if any(p1 < lamda)
                        n(i) = 1;
                    end
                end
            end 
            % 处理m==1的节点，A节点
            m_one_idx = find(m == 1);%存储A节点所在的列
             if ~isempty(m_one_idx)%判断存储的元组是否为空，不为空则说明A节点存在
                % AS 或 AG 状态的更新
                % 向量化处理x==0和x==1的情况
                AS_idx = m_one_idx(x(m_one_idx) == 0);%AS节点
                AG_idx = m_one_idx(x(m_one_idx) == 1);%x(m_one_idx) == 1判断下层节点状态的矩阵x中，是1的节点，即G节点
                
                % AS 更新为 US 或保持 AS
                if ~isempty(AS_idx)
                    p1_AS = rand(1, length(AS_idx));%
                    n(AS_idx) = p1_AS >= delta;%如果 p1_AS(i) >= delta，则该表达式返回 true（逻辑值 1）；否则，返回 false（逻辑值 0）
                    %n(AS_idx) 表示更新 n 向量在 AS_idx 位置上的值，将这些位置的元素更新为逻辑值 1 或0
                    %，当p1_AS >= delta则赋值1，表示AS状态，当p1_AS<delta,则发生变化为US
                end
                % AG 更新为 AG 或 UG
                if ~isempty(AG_idx)
                    p1_AG = rand(1, length(AG_idx));
                    n(AG_idx) = p1_AG >= delta;%AG>= delta则赋值为1，，即AG不变，否则变成UG
                end
            end
            
            % 下层传播
            y = x;  % 先复制x的状态
            s_idx = find(x == 0);%下层节点中的S,存在AS,US
            g_idx = find(x == 1);%下层节点中的G,即AG，和上层传播留下的UG
            
            % 处理 S 状态的节点
            for i = s_idx%AS+US
                if n(i) == 1  % 这里指的是AS节点，操作AS -> AG
                    neighbors = B_neighbors{i};%AS节点的邻居节点
                    active_neighbors = neighbors(x(neighbors) == 1);%确定AS的邻居节点中的AG节点
                    if ~isempty(active_neighbors)
                        p1 = rand(1, length(active_neighbors));
                        if any(p1 < gama * beta_R)
                            y(i) = 1;
                        end
                    end
                else  %这里指的是us节点 US -> AG
                    neighbors = B_neighbors{i};
                    active_neighbors = neighbors(x(neighbors) == 1);
                    if ~isempty(active_neighbors)%检查是否不为空，不为空则进入条件分支
                        p0 = rand(1, length(active_neighbors));
                        if any(p0 < beta_R)
                            y(i) = 1;
                            n(i) = 1;
                        end
                    end
                end
            end
            
            % 处理 G 状态的节点
            for i = g_idx
                if n(i) == 1  % AG
                    p1 = rand;
                    if p1 < mu
                        y(i) = 0;  % AG 变为 AS
                    else
                        y(i) = 1;  % 保持 AG
                    end
                else  % UG
                    p0 = rand;
                    if p0 < mu
                        y(i) = 0;  % UG 变为 US
                    else
                        y(i) = 1;  % UG 变为 AG
                        n(i) = 1;
                    end
                end
            end
        %政策调控作用，US 和 AS 节点以 m 的概率直接变成 AG
        %先找到S节点
        s2_idx = find(y == 0);%下层节点中的S,存在AS,US
        %先处理AS节点
        for i=s2_idx
        if n(i)==1
            p1=rand;
            if p1<policy
                y(i)=1;
            end
        else
            p0=rand;
            if p0<policy
                y(i)=1;
                n(i)=1;
            end
        end
        end
            % 更新状态
            m = n;
            x = y;
        end  % end of stp
        
        % 统计结果
        nodeA = sum(m == 1);%A节点的总数
        BETA_A = nodeA / N;
        BETA_AG(rep, l) = sum(x) / N;%AG的密度
        BETA_AS(rep, l) = BETA_A - BETA_AG(rep, l);%AS的密度
        BETA_US(rep, l) = 1 - BETA_AS(rep, l) - BETA_AG(rep, l);%US密度
    end  % end of termi
end  % end of MC_rep

% 计算多次仿真的平均值
BETA_AS_AVER = mean(BETA_AS, 1);
BETA_AG_AVER = mean(BETA_AG, 1);
BETA_US_AVER = mean(BETA_US, 1);

% 绘图
xzhou = (1:termi) / 80;
figure;
hold on;
box on;
grid off;
set(gca, 'FontSize', 15);
plot(xzhou, BETA_AG_AVER, '-o', 'Color', [46/256 106/256 89/256]);
plot(xzhou, BETA_AS_AVER, '-^', 'Color', [176/256 76/256 72/256]);
plot(xzhou, BETA_US_AVER, '-v', 'Color', [220/256 178/256 105/256]);
set(gcf, 'DefaultTextInterpreter', 'latex');
xlabel('t', 'FontSize', 15);
ylabel('proportion', 'FontSize', 15);
legend({'$\rho^{AG}$', '$\rho^{AS}$', '$\rho^{US}$'}, 'Interpreter', 'latex', 'FontSize', 15);
hold off; 
save('BETA_AS_AVER.mat','BETA_AS_AVER')
save('BETA_AG_AVER.mat','BETA_AG_AVER')
save('BETA_US_AVER.mat','BETA_US_AVER') 
saveas(gcf, 'mc.fig')