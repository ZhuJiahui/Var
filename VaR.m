clear all;
clc;

N = 20;  % 股票数
n = 10;  % 迭代天数
RS = 200 * randn(n, N);  % 投资回报
S = 50 * rand(n, N);  % 股票价格
portfolio = zeros(n, N);  % 投资策略
A = 10000000;  % 初始投资金额

k = 1;  % 迭代次数

this_lambda = 1;  % 解参数
mu = 1000;  % 惩罚因子
c = 10;  % 期望回报  %该值的设定很重要
alpha = 0.95;  % 置信水平

rho_k = 1e-5;  % 近似度参数rho初始值

pre_weights = rand(1, N);
weights = pre_weights / sum(pre_weights);  % 初始各股权重
beta = -0.5;  % beta初始值

% 整个过程的值均保留
all_rho = zeros(n, 1);
all_weights = zeros(n, N);
all_beta = zeros(n, 1);

gamma_k = sqrt(2 / rho_k + 1 / this_lambda);
kappa_k = 1 / this_lambda / gamma_k;
while (abs(gamma_k - kappa_k) > 1e-40  && k <= n)
    
    % H_lambda要先累加之前的各个g_lambda
    H_lambda1 = 0;
    for mm = 1 : k - 1
        g_lambda_z = - all_weights(mm, :) * RS(mm, :)' - all_beta(mm);
        g_lambda_value = g_lambda_func(g_lambda_z, all_rho(mm), this_lambda);
        H_lambda1 = H_lambda1 + g_lambda_value;
    end
    
    avg_return = mean(RS(1 : k, :), 1);  % 之前若干天的平均回报
    
    % 不等式线性约束
    AA_pre = -1 * avg_return;
    A_in1 = [AA_pre, 0];
    A_in2 = [RS(k, :), 1];
    A_in3 = [-1 * RS(k, :), -1];
    
    % 四个段各自的不等式约束
    A_ineq1 = [A_in1; A_in2];
    b_ineq1 = [-c; -gamma_k];
    A_ineq2 = [A_in1; A_in2; A_in3];
    b_ineq2 = [-c; -kappa_k; gamma_k];
    %A_ineq2 = A_in1;
    %b_ineq2 = [-c];
    A_ineq3 = [A_in1; A_in2; A_in3];
    b_ineq3 = [-c; 0; kappa_k];
    A_ineq4 = [A_in1; A_in3];
    b_ineq4 = [-c; 0];
    
    
    % 等式线性约束
    A_eq = [ones(1, N), 0];
    b_eq = [1];
    
    % 初始值
    x0 = [weights, beta];
    
    % 下界
    lb = [zeros(1, N), -1000];
    % 上界
    ub = [0.1 * ones(1, N), 1000];
    options = optimset('Algorithm', 'active-set');  % 优化算法
    
    % 非线性优化
    [x1, fvalue1, exitflag1, output1] = fmincon(@(x)(x(N + 1) + mu * max((H_lambda1 + 1) / k - alpha, 0)), x0, A_ineq1, b_ineq1, A_eq, b_eq, lb, ub);
    
    [x2, fvalue2, exitflag2, output2] = fmincon(@(x)(x(N + 1) + mu * max((H_lambda1 + 1 - 0.5 * rho_k * (- x(1 : N) * RS(k, :)' - x(N + 1) - gamma_k)^2) / k - alpha, 0)), x0, A_ineq2, b_ineq2, A_eq, b_eq, lb, ub, [], options);
    
    [x3, fvalue3, exitflag3, output3] = fmincon(@(x)(x(N + 1) + mu * max((H_lambda1 + this_lambda * (- x(1 : N) * RS(k, :)' - x(N + 1))^2) / k - alpha, 0)), x0, A_ineq3, b_ineq3, A_eq, b_eq, lb, ub);
    
    [x4, fvalue4, exitflag4, output4] = fmincon(@(x)(x(N + 1) + mu * max((H_lambda1 + 0) / k - alpha, 0)), x0, A_ineq4, b_ineq4, A_eq, b_eq, lb, ub);
    
    % 取四者最小
    all_x = [x1; x2; x3; x4];
    all_fvalue = [fvalue1, fvalue2, fvalue3, fvalue4];
    [min_value, min_index] = min(all_fvalue);
    x = all_x(min_index, :);
    
    %fprintf('%f\n', min_index);
    %x = x2;
    %fprintf('%d\n', exitflag2);
    % 权值和beta的更新
    weights = x(1 : N);
    beta = x(N + 1);
    all_weights(k, :) = weights;
    all_beta(k) = beta;
    all_rho(k) = rho_k;
    
    % 计算投资策略
    if (k > 1)
        shares = (weights * 0.1 * A) ./ (S(k - 1, :));
        for nn = 1 : length(shares)
            if (shares(nn) >= 100)
                shares(nn) = ceil(shares(nn));
            else
                shares(nn) = 0;
            end
        end
        portfolio(k, :) = shares;
    end
    
     % 验证H_lambda-alpha是否大于0
    H_lambda2 = 0;
    for mm = 1 : k
        g_lambda_z2 = - all_weights(mm, :) * RS(mm, :)' - all_beta(mm);
        g_lambda_value2 = g_lambda_func(g_lambda_z2, all_rho(mm), this_lambda);
        H_lambda2 = H_lambda2 + g_lambda_value2;
    end
    %fprintf('%f', H_lambda2);
    
    k = k + 1;
    rho_k = 3 * rho_k;
    
    gamma_k = sqrt(2 / rho_k + 1 / this_lambda);
    kappa_k = 1 / this_lambda / gamma_k;
    %fprintf('%f\n', gamma_k);
    %fprintf('%f\n', kappa_k);
end


