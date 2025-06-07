% function miafis_demo()
% MIAFIS算法演示 - 复现论文6.3.4节结果

% 参数设置
clc;clear;close all;
% rng default;
tic;
fs = 100e3;             % 采样频率 (Hz)
%！！！！！！ 目前的极限是250 * 200 ！！！！！！
N = 100;
Nv = 100;         % 多普勒bin数量

f_max = 500;           % 优化范围的最高频率
doppler_bins = linspace(-f_max, f_max, Nv+1);
vh_range = doppler_bins / fs; 
target_freq = 60;  % 目标频率Hz
[~, freq_bin_idx] = min(abs(doppler_bins - target_freq));  % 找最近的bin

gamma = 4;       % PAR参数
max_iter = 1e3; % 最大迭代次数
tol = 1e-6;      % 收敛容差
%

% 初始化随机波形和目标函数值
s = exp(1j * 2*pi * rand(N,1));
obj_values = zeros(max_iter, 1);
%

% 设置干扰区域(p_k) (bin从0开始, MATLAB索引从1开始)
interference_map = zeros(N, Nv);
interference_map(1:N, freq_bin_idx) = 1;
interference_map(20:N, Nv/2+1) = 1;    % 0 Hz
%
plot_ambiguity_function(s, N, Nv, interference_map, vh_range);

% data = load("100_100_1e4.mat", "s");
% s_generate = data.s;
% plot_ambiguity_function(s_generate, N, Nv, interference_map, vh_range);
% figure;plot(real(s_generate));


% 预计算A_k矩阵，N*Nv*N*N, 25x50x{25x25 double}
A = cell(N, Nv);
fprintf('A 占用内存 %.2f MB\n', whos('A').bytes / (1024^2));
for r = 0:N-1
    for h = 0:Nv-1
        vh = vh_range(h+1);   % vh = -0.5 + h/Nv; % 归一化多普勒频率
        p_vec = exp(1j*2*pi*vh*(0:N-1)'); % 多普勒相位向量
        A{r+1,h+1} = circshift(diag(p_vec), r); % 时延+多普勒矩阵 % circshift(s,[0,r]) 是将序列 s 循环按列移位 r 个位置。
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

hh = waitbar(0, '进度'); 
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

    % if mod(iter, 100) == 0
    %     elapsed_time = toc; % 计算从tic以来的耗时
    %     fprintf('第%d次迭代，已耗时 %.2f 秒。\n', iter, elapsed_time);
    %     tic; % 重新启动计时器（计下一轮的100次）
    % end

    s = s_new;

    waitbar(iter/max_iter, hh, sprintf('进度: %d%%', round(iter/max_iter*100)));
end
close(hh); % 关闭进度条


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
[AF, delay_range, doppler_range] = plot_ambiguity_function(s, N, Nv, interference_map, vh_range);

fprintf('算法完成。最终目标值: %.4f\n', obj_values(end));
toc;

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

function [AF, delay_range, doppler_range] = plot_ambiguity_function(s, N, Nv, p, vh_range)
    % 绘制模糊函数
    fs = 100e3;
    doppler_Hz = vh_range * fs;  % 横坐标：Hz
    normalize_flag = 0;

    % 计算模糊函数
    delay_range = 0:N-1;
    if normalize_flag==1
        doppler_range = -0.5:1/Nv:0.5-1/Nv;
    else
        doppler_range = vh_range;
    end
    
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
    % figure;
    % [X, Y] = meshgrid(doppler_range, delay_range);
    % % surf(X, Y, 10*log10(AF), 'EdgeColor', 'none');
    % surf(X, Y, AF, 'EdgeColor', 'none');
    % xlabel('归一化多普勒频率');
    % ylabel('时延(sample)');
    % zlabel('幅度');
    % title('设计的模糊函数');
    % view(30, 45);
    % colorbar;
    
    % 2D切片
    figure;
    subplot(2,2,1);
    if normalize_flag==1
        imagesc(doppler_range, delay_range, 10*log10(AF));
        xlabel('归一化多普勒频率');
    else
        imagesc(doppler_Hz, delay_range, 10*log10(AF));
        xlabel('多普勒频率(Hz)');
    end
    
    
    ylabel('时延(sample)');
    title('模糊函数');
    % clim([-40 20]);
    colorbar;
    hold on;
    % 检测p中的连通区域
    CC = bwconncomp(p);
    bin_width = doppler_Hz(2) - doppler_Hz(1);
    for i = 1:CC.NumObjects
        % 得到所有该区域的索引
        inds = CC.PixelIdxList{i};
        [r_idx, v_idx] = ind2sub(size(p), inds);
        r_min = min(r_idx)-1; r_max = max(r_idx)-1; % -1转为实际时延
        v_min = min(v_idx)-1; v_max = max(v_idx)-1; % -1转为实际doppler bin
        
        % 映射到坐标轴（中心居中, 所以-0.5/+0.5）
        if normalize_flag == 1
            x_left = doppler_range(v_min+1) - 0.5/Nv;
            y_bottom = delay_range(r_min+1) - 0.5;
            w = (v_max - v_min + 1) / Nv;
            h = (r_max - r_min + 1);
        else
            x_left = doppler_Hz(v_min+1) - bin_width/2;
            x_right = doppler_Hz(v_max+1) + bin_width/2;
            w = x_right - x_left;
            y_bottom = delay_range(r_min+1) - 0.5;
            h = (r_max - r_min + 1);
        end
        
        rectangle('Position', [x_left, y_bottom, w, h], ...
            'EdgeColor', 'r', 'LineWidth', 0.5, 'LineStyle', '-');
    end
    hold off;
    
    subplot(2,2,2);
    plot(delay_range, 10*log10(AF(:, round(Nv/2))), 'LineWidth', 2);
    xlabel('时延(sample)');
    ylabel('幅度(dB)');
    title('零多普勒切面');
    grid on;
    
    subplot(2,2,3);
    if normalize_flag==1
        plot(doppler_range, 10*log10(AF(3, :)), 'LineWidth', 2);
        xlabel('归一化多普勒频率');
    else
        plot(doppler_Hz, 10*log10(AF(3, :)), 'LineWidth', 2);
        xlabel('多普勒频率(Hz)');
    end
    ylabel('幅度(dB)');
    title('时延=2切面');
    grid on;
    
    subplot(2,2,4);
    if normalize_flag==1
        plot(doppler_range, 10*log10(AF(50, :)), 'LineWidth', 2);
        xlabel('归一化多普勒频率');
    else
        plot(doppler_Hz, 10*log10(AF(50, :)), 'LineWidth', 2);
        xlabel('多普勒频率(Hz)');
    end
    ylabel('幅度(dB)');
    title('时延=50切面');
    grid on;
end