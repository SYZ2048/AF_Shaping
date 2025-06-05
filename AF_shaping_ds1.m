% function miafis_demo()
% MIAFIS算法演示 - 复现论文6.3.4节结果

% 参数设置
clc;clear;close all;
% rng default;
N = 25;          % 波形长度
Nv = 50;         % 多普勒bin数量
gamma = 4;       % PAR参数
max_iter = 1e4; % 最大迭代次数
tol = 1e-6;      % 收敛容差
%

% 初始化随机波形
s = exp(1j * 2*pi * rand(N,1));
% 初始化记录目标函数值
obj_values = zeros(max_iter, 1);
%

% 设置干扰区域(p_k) (bin从0开始, MATLAB索引从1开始)
interference_map = zeros(N, Nv);
interference_map(3:4, 36:38) = 1; % 
interference_map(4:4, 19:20) = 1;
interference_map(2:N, 26) = 1;
%
plot_ambiguity_function(s, N, Nv, interference_map);

% 预计算A_k矩阵，N*Nv*N*N, 25x50x{25x25 double}
A = cell(N, Nv);
for r = 0:N-1
    for h = 0:Nv-1
        vh = -0.5 + h/Nv; % 归一化多普勒频率
        p_vec = exp(1j*2*pi*vh*(0:N-1)'); % 多普勒相位向量
        % circshift(s,[0,r]) 是将序列 s 循环按列移位 r 个位置。
        A{r+1,h+1} = circshift(diag(p_vec), r); % 时延+多普勒矩阵
    end
end

% 计算lambda_u(B)
lambda_u_B = 0;
for r = 0:N-1
    for h = 0:Nv-1
        if interference_map(r+1,h+1) > 0
            lambda_u_B = lambda_u_B + interference_map(r+1,h+1) * (N - r);
        end
    end
end

% MIAFIS主循环
for iter = 1:max_iter
    S = s * s';

    % 计算P矩阵
    P = zeros(N);
    for r = 0:N-1
        for h = 0:Nv-1
            if interference_map(r+1,h+1) > 0
                Ak = A{r+1,h+1};
                P = P + interference_map(r+1,h+1)/2 * (trace(Ak'*S)*Ak + trace(Ak*S)*Ak');
            end
        end
    end

    % 计算lambda_u(P)
    lambda_u_P = 0;
    for r = 0:N-1
        for h = 0:Nv-1
            if interference_map(r+1,h+1) > 0
                Ak = A{r+1,h+1};
                if r > 0
                    term = interference_map(r+1,h+1) * sqrt((N-r)*(N-1)/(2*N)) * abs(s'*Ak*s);
                else
                    term = interference_map(r+1,h+1) * N;
                end
                lambda_u_P = lambda_u_P + term;
            end
        end
    end

    % 计算z向量
    z = (P - (lambda_u_B*N + lambda_u_P)*eye(N)) * s;

    % 投影操作 (PAR和能量约束)
    s_new = project_to_constraints(z, gamma, N);

    % 计算目标函数值
    obj = 0;
    for r = 0:N-1
        for h = 0:Nv-1
            if interference_map(r+1,h+1) > 0
                Ak = A{r+1,h+1};
                obj = obj + interference_map(r+1,h+1) * abs(s'*Ak*s)^2;
            end
        end
    end
    obj_values(iter) = obj;

    % 检查收敛
    if norm(s_new - s) / norm(s) < tol
        fprintf('收敛于迭代 %d\n', iter);
        break;
    end

    s = s_new;
end

% 截断记录
obj_values = obj_values(1:iter);

% 绘制收敛曲线
figure;
semilogx(1:iter, obj_values, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('目标函数值');
title('MIAFIS收敛曲线');
grid on;

% 绘制模糊函数
[AF, delay_range, doppler_range] = plot_ambiguity_function(s, N, Nv, interference_map);

fprintf('算法完成。最终目标值: %.4f\n', obj_values(end));


function s_proj = project_to_constraints(z, gamma, N)
    % 投影操作，满足PAR和能量约束
    
    % 计算幅度
    abs_z = abs(z);
    
    % 找到满足PAR约束的投影
    beta_low = 0;
    beta_high = sqrt(gamma) / min(abs_z(abs_z > 0));
    
    % 二分搜索找到合适的beta
    for k = 1:20 % 最多20次二分搜索
        beta = (beta_low + beta_high) / 2;
        s_abs = min(beta * abs_z, sqrt(gamma));
        total_power = sum(s_abs.^2);
        
        if total_power > N
            beta_high = beta;
        else
            beta_low = beta;
        end
    end
    
    % 最终投影
    s_abs = min(beta * abs_z, sqrt(gamma));
    s_proj = sqrt(N) * (s_abs .* exp(1j * angle(z))) / norm(s_abs);
end

function [AF, delay_range, doppler_range] = plot_ambiguity_function(s, N, Nv, p)
    % 绘制模糊函数
    
    % 计算模糊函数
    delay_range = 0:N-1;
    doppler_range = -0.5:1/Nv:0.5-1/Nv;
    AF = zeros(length(delay_range), length(doppler_range));
    
    for i = 1:length(delay_range)
        r = delay_range(i);
        for j = 1:length(doppler_range)
            v = doppler_range(j);
            p_vec = exp(1j*2*pi*v*(0:N-1)');
            J = circshift(eye(N), r);
            AF(i,j) = abs(s' * J * diag(p_vec) * s)^2 / norm(s)^2;
        end
    end
    
    % 3D绘图
    figure;
    [X, Y] = meshgrid(doppler_range, delay_range);
    % surf(X, Y, 10*log10(AF), 'EdgeColor', 'none');
    surf(X, Y, AF, 'EdgeColor', 'none');
    xlabel('归一化多普勒频率');
    ylabel('时延(sample)');
    zlabel('幅度');
    title('设计的模糊函数');
    view(30, 45);
    colorbar;
    
    % 2D切片
    figure;
    subplot(2,2,1);
    imagesc(doppler_range, delay_range, 10*log10(AF));
    xlabel('归一化多普勒频率');
    ylabel('时延(sample)');
    title('模糊函数');
    clim([-30 20]);
    colorbar;
    hold on;
    % 检测p中的连通区域
    CC = bwconncomp(p);
    stats = regionprops(CC, 'BoundingBox');
    for i = 1:length(stats)
        bb = stats(i).BoundingBox; % 格式 [x, y, width, height]，索引起点为1
        % 将p矩阵坐标转换为plot的横纵坐标
        % p的行是时延(r)，列是多普勒bin
        % bb(1): 列起点，bb(2): 行起点
        x_left = doppler_range(floor(bb(1)));
        y_top = delay_range(floor(bb(2)));
        w = bb(3) * (doppler_range(2) - doppler_range(1)); % 水平宽度
        h = bb(4) * (delay_range(2) - delay_range(1));     % 垂直高度
        rectangle('Position', [x_left, y_top, w, h], ...
            'EdgeColor', 'r', 'LineWidth', 0.5);
    end
    hold off;
    
    subplot(2,2,2);
    plot(delay_range, 10*log10(AF(:, round(Nv/2))), 'LineWidth', 2);
    xlabel('时延(sample)');
    ylabel('幅度(dB)');
    title('零多普勒切面');
    grid on;
    
    subplot(2,2,3);
    plot(doppler_range, 10*log10(AF(3, :)), 'LineWidth', 2);
    xlabel('归一化多普勒频率');
    ylabel('幅度(dB)');
    title('时延=2切面');
    grid on;
    
    subplot(2,2,4);
    plot(doppler_range, 10*log10(AF(4, :)), 'LineWidth', 2);
    xlabel('归一化多普勒频率');
    ylabel('幅度(dB)');
    title('时延=3切面');
    grid on;
end