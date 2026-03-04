% function miafis_demo()
% MIAFIS算法演示 - 复现论文6.3.4节结果

clc;clear;
close all;
rng default;
tic;

% 参数设置 from main.m
fs_main = 100e3;
pulse_duration = 0.1;
prf = 1/pulse_duration;
Bw = 1e3;
c0 = 1500;
signal_len = pulse_duration * fs_main;
up_sample_rate = fs_main / Bw;
N_len = signal_len  / up_sample_rate;   %！！！！！！ 目前的极限是250 * 200 ！！！！！！
distance = 20;
delta_delay = distance /c0;
delta_t = pulse_duration / N_len;

% 波形优化参数
fs = fs_main / up_sample_rate;             % 采样频率 (Hz)，需要对应main.m中的fs，fs = fs_main * N / signal_len_main
Nv = 100;         % 多普勒bin数量
f_max = 200;           % 优化范围的最高频率
doppler_bins = linspace(-f_max, f_max, Nv+1);
vh_range = doppler_bins / fs; 
target_freq = 60;  % 目标频率Hz
[~, freq_bin_idx] = min(abs(doppler_bins - target_freq));  % 找最近的bin
gamma = 1.5;       % PAR参数
max_iter = 5e3; % 最大迭代次数
tol = 3e-4;      % 收敛容差

delay_bin_idx = round(delta_delay / delta_t);
range_span = 6; % m
range_span_idx = round(range_span/c0 / delta_t );

% 设置干扰区域(p_k) (bin从0开始, MATLAB索引从1开始)
interference_map = zeros(N_len, Nv+1);
interference_map(1:N_len, freq_bin_idx:freq_bin_idx) = 1;   % Q
interference_map(20:N_len, Nv/2+1) = 1;
% interference_map(delay_bin_idx-range_span_idx:delay_bin_idx+range_span_idx+1, freq_bin_idx-5:freq_bin_idx+5) = 1;   % Q
% interference_map(delay_bin_idx-range_span_idx:delay_bin_idx+range_span_idx+1, Nv/2+1) = 1;    % 0 Hz


Algorithm = 'local';  % 'original', 'squarem', 'local'


s = exp(1j * 2*pi * rand(N_len,1)); % 初始化随机波形和目标函数值
s_random = s;
%%


% f0 = 0 - Bw/2;
% t_pulse = 0:1/fs:pulse_duration-1/fs; 
% s = exp(1i*pi*(Bw/pulse_duration)*t_pulse.^2) .* exp(1i*2*pi*f0*t_pulse);
% s = s.';
obj_values = zeros(max_iter, 1);

% data = load("100_100_5e3_local_ISLQ_60Hz.mat", "s", "interference_map");
% s_generate = data.s;
% % interference_map = data.interference_map;
% plot_ambiguity_function(s_generate, N_len, Nv, interference_map, vh_range, fs);
% figure;plot(real(s_generate));

plot_ambiguity_function(s, N_len, Nv, interference_map, vh_range, fs);

%% 
% 预计算A_k矩阵，N*Nv*N*N, 25x50x{25x25 double}
A = cell(N_len, Nv);
fprintf('A 占用内存 %.2f MB\n', whos('A').bytes / (1024^2));
for r = 0:N_len-1
    for h = 0:Nv-1
        vh = vh_range(h+1);   % vh = -0.5 + h/Nv; % 归一化多普勒频率
        p_vec = exp(1j*2*pi*vh*(0:N_len-1)'); % 多普勒相位向量
        A{r+1,h+1} = circshift(diag(p_vec), r); % 时延+多普勒矩阵 % circshift(s,[0,r]) 是将序列 s 循环按列移位 r 个位置。
    end
end

% 计算lambda_u(B)
lambda_u_B = 0;
for r = 0:N_len-1
    for h = 0:Nv-1
        if interference_map(r+1,h+1) > 0
            lambda_u_B = lambda_u_B + interference_map(r+1,h+1) * (N_len - r);
        end
    end
end

hh = waitbar(0, '进度'); 
% MIAFIS主循环
for iter = 1:max_iter
    switch Algorithm
        case 'original'
            s_new = miafis_update(s, A, interference_map, lambda_u_B, N_len, Nv, gamma);
            
        case 'squarem'
            if iter == 1
                s_prev = s;
                s_prev_obj = compute_objective(s_prev, A, interference_map, N_len, Nv);
            end
            s_new = squarem_acceleration(s, A, interference_map, lambda_u_B, N_len, Nv, gamma, @miafis_update, @compute_objective, s_prev, s_prev_obj);
            s_prev = s;
            s_prev_obj = compute_objective(s, A, interference_map, N_len, Nv);
            
        case 'local'
            s_new = local_majorization_acceleration(s, A, interference_map, lambda_u_B, N_len, Nv, gamma, @compute_objective);
    end

    % 计算目标函数值
    obj_values(iter) = compute_objective(s, A, interference_map, N_len, Nv);

    % 检查收敛
    if norm(s_new - s) / norm(s) < tol
        fprintf('收敛于迭代 %d\n', iter);
        break;
    end

    s = s_new;

    waitbar(iter/max_iter, hh, sprintf('进度: %d%%', round(iter/max_iter*100)));
end
close(hh); % 关闭进度条


% 截断记录
obj_values = obj_values(1:iter);
fprintf('算法完成。最终目标值: %.4f\n', obj_values(end));

% 绘制收敛曲线
figure;
semilogx(1:iter, obj_values, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('目标函数值');
title('MIAFIS收敛曲线');
grid on;

[AF, delay_range, doppler_range] = plot_ambiguity_function(s, N_len, Nv, interference_map, vh_range, fs);

[af,delay,doppler] = ambgfun(s, fs, prf);
[af_random,delay_random,doppler_random] = ambgfun(s_random, fs, prf);
figure();
subplot(2,1,1);
imagesc(delay*1e3,doppler,mag2db(af));
ylim([-100 100]);
xlim([-30 30]);
clim([-50 0])
cbar = colorbar;
cbar.Label.String = '(dB)';
axis xy
xlabel('Delay \tau (ms)')
ylabel('Doppler f_d (Hz)')
title('Optimized AF');
subplot(2,1,2);
imagesc(delay_random*1e3,doppler_random,mag2db(af_random));
ylim([-100 100]);
xlim([-30 30]);
clim([-50 0])
cbar = colorbar;
cbar.Label.String = '(dB)';
axis xy
xlabel('Delay \tau (ms)')
ylabel('Doppler f_d (Hz)')
title('Random AF');

figure();
[afmag, delay] = ambgfun(s, fs, prf,Cut='Doppler', CutValue=0);
[afmag_random, ~] = ambgfun(s_random, fs, prf,Cut='Doppler', CutValue=0);
afmag_db = mag2db(afmag);
afmag_random_db = mag2db(afmag_random);
plot(delay*1500, afmag_db, 'LineWidth', 1);hold on;
plot(delay*1500, afmag_random_db, 'LineWidth', 1);
% xlabel('Delay \tau (ms)');
xlabel('Range (m)')
ylabel('Response');
title('匹配滤波结果（多普勒切面）');
legend('优化波形','随机波形','Location','best');

grid on;


%% 绘制模糊函数
filename = 'range7_13.mat';
% data = load(filename, "s");
% s = data.s;
% ----------------------- Test code start    ----------------------- 
% rho = 1.9;               % Variable exponent parameter
% alpha = 1e3;             % Frequency modulation term
% prf = 1 / pulse_duration;
% fc = 0;
% gsfm_finst = @(t)(fc + (Bw/2)*sin(2*pi*alpha*t.^rho));
% gsfmwav = phased.CustomFMWaveform('SampleRate',fs,...
%     'PulseWidth',pulse_duration,'FrequencyModulation',gsfm_finst,'PRF',prf);
% gsfm_samples = gsfmwav();
% s = gsfm_samples;
% ----------------------- Test code end      ----------------------- 

[AF, delay_range, doppler_range] = plot_ambiguity_function(s, N_len, Nv, interference_map, vh_range, fs);

toc;

save(filename);   % 只保存变量 varName

%% -------------------------------- Test
% 读取特定文件中的波形，和随机波形进行模糊函数性能比对
clc;
% clear;
close all;
rng default;

filename = 'WISL1.mat';
data = load(filename);
s = data.xWeCAN;
[AF, delay_range, doppler_range] = plot_ambiguity_function(s_random, N_len, Nv, interference_map, vh_range, fs);
[AF, delay_range, doppler_range] = plot_ambiguity_function(s, N_len, Nv, interference_map, vh_range, fs);


[af,delay,doppler] = ambgfun(s, fs, prf);
[af_random,delay_random,doppler_random] = ambgfun(s_random, fs, prf);
figure();
subplot(2,1,1);
imagesc(delay*1e3,doppler,mag2db(af));
ylim([-100 100]);
xlim([-30 30]);
clim([-50 0])
cbar = colorbar;
cbar.Label.String = '(dB)';
axis xy
xlabel('Delay \tau (ms)')
ylabel('Doppler f_d (Hz)')
title('Optimized AF');
subplot(2,1,2);
imagesc(delay_random*1e3,doppler_random,mag2db(af_random));
ylim([-100 100]);
xlim([-30 30]);
clim([-50 0])
cbar = colorbar;
cbar.Label.String = '(dB)';
axis xy
xlabel('Delay \tau (ms)')
ylabel('Doppler f_d (Hz)')
title('Random AF');


figure();
[afmag, delay] = ambgfun(s, fs, prf,'Cut','Doppler');
[afmag_random, ~] = ambgfun(s_random, fs, prf,'Cut','Doppler');
afmag_db = mag2db(afmag);
afmag_random_db = mag2db(afmag_random);
plot(delay*1500, afmag_db, 'LineWidth', 1);hold on;
plot(delay*1500, afmag_random_db, 'LineWidth', 1);
% xlabel('Delay \tau (ms)');
xlabel('Range (m)')
ylabel('Response');
xlim([-40 40])
title('匹配滤波结果（0多普勒切面）');
legend('优化波形','随机波形','Location','best');
grid on;

% -------------------------------- Test Code End --------------------------------


function s_new = miafis_update(s, A, interference_map, lambda_u_B, N, Nv, gamma)
    % MIAFIS核心更新步骤
    S = s * s';
    P = zeros(N);
    lambda_u_P = 0;
    
    for r = 0:N-1
        for h = 0:Nv-1
            if interference_map(r+1,h+1) > 0
                Ak = A{r+1,h+1};
                P = P + interference_map(r+1,h+1)/2 * (trace(Ak'*S)*Ak + trace(Ak*S)*Ak');
                
                if r > 0
                    term = interference_map(r+1,h+1) * sqrt((N-r)*(N-1)/(2*N)) * abs(s'*Ak*s);
                else
                    term = interference_map(r+1,h+1) * N;
                end
                lambda_u_P = lambda_u_P + term;
            end
        end
    end

    z = (P - (lambda_u_B*N + lambda_u_P)*eye(N)) * s;
    s_new = project_to_constraints(z, gamma, N);
end

function s_new = squarem_acceleration(s, A, interference_map, lambda_u_B, N, Nv, gamma, update_func, obj_func, s_prev, s_prev_obj)
    % SQUAREM加速方法
    s1 = update_func(s, A, interference_map, lambda_u_B, N, Nv, gamma);
    s2 = update_func(s1, A, interference_map, lambda_u_B, N, Nv, gamma);
    
    q = s1 - s;
    v = s2 - s1 - q;
    alpha = -norm(q)/norm(v);
    
    % 初始更新
    s_new = project_to_constraints(s - 2*alpha*q + alpha^2*v, gamma, N);
    new_obj = obj_func(s_new, A, interference_map, N, Nv);
    
    % 回溯直到单调性满足
    while new_obj > s_prev_obj
        alpha = (alpha - 1)/2;
        s_new = project_to_constraints(s - 2*alpha*q + alpha^2*v, gamma, N);
        new_obj = obj_func(s_new, A, interference_map, N, Nv);
        
        if abs(alpha + 1) < 1e-6
            break;  % 防止无限循环
        end
    end
end

function s_new = local_majorization_acceleration(s, A, interference_map, lambda_u_B, N, Nv, gamma, obj_func)
    % 局部优化加速方法
    % 初始t值较小
    t = 0.1 * (lambda_u_B * N);  % 从小的t值开始
    
    while true
        [s_new, P, lambda_u_P] = miafis_update_with_P(s, A, interference_map, lambda_u_B, N, Nv, gamma);
        z = (P - t*eye(N)) * s;
        s_candidate = project_to_constraints(z, gamma, N);
        
        % 检查单调性
        current_obj = obj_func(s, A, interference_map, N, Nv);
        candidate_obj = obj_func(s_candidate, A, interference_map, N, Nv);
        
        if candidate_obj <= current_obj || t >= (lambda_u_B * N + lambda_u_P)
            s_new = s_candidate;
            break;
        else
            t = min(2 * t, (lambda_u_B * N + lambda_u_P));  % 增加t值
        end
    end
end

function [s_new, P, lambda_u_P] = miafis_update_with_P(s, A, interference_map, lambda_u_B, N, Nv, gamma)
    % 返回P和lambda_u_P的MIAFIS更新
    S = s * s';
    P = zeros(N);
    lambda_u_P = 0;
    
    for r = 0:N-1
        for h = 0:Nv-1
            if interference_map(r+1,h+1) > 0
                Ak = A{r+1,h+1};
                P = P + interference_map(r+1,h+1)/2 * (trace(Ak'*S)*Ak + trace(Ak*S)*Ak');
                
                if r > 0
                    term = interference_map(r+1,h+1) * sqrt((N-r)*(N-1)/(2*N)) * abs(s'*Ak*s);
                else
                    term = interference_map(r+1,h+1) * N;
                end
                lambda_u_P = lambda_u_P + term;
            end
        end
    end
    
    z = (P - (lambda_u_B*N + lambda_u_P)*eye(N)) * s;
    s_new = project_to_constraints(z, gamma, N);
end

function obj = compute_objective(s, A, interference_map, N, Nv)
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
end

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

function [AF, delay_range, doppler_range] = plot_ambiguity_function(s, N, Nv, p, vh_range, fs)
    % 绘制模糊函数
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

    %%%%%%% [新增功能] 计算interference_map区域内的旁瓣均值 %%%%%
    if nargin >= 4 && ~isempty(p)
        % 确保 p 是逻辑掩码
        mask = logical(p);
        
        % 检查维度一致性
        if isequal(size(AF), size(mask))
            % 提取掩码区域内的功率值
            vals_in_region = AF(mask);
            
            if ~isempty(vals_in_region)
                % 1. 计算线性域的平均功率
                mean_power = mean(vals_in_region);
                % 2. 转换为 dB
                mean_sidelobe_dB = 10 * log10(mean_power);
                
                % 3. 命令行输出
                fprintf('\n------------------------------------------------\n');
                fprintf('【统计结果】\n');
                fprintf('  Interference Map 区域内的旁瓣均值: %.4f dB\n', mean_sidelobe_dB);
                fprintf('------------------------------------------------\n\n');
            else
                disp('提示: Interference Map 区域为空，无法计算均值。');
            end
        else
            warning('Interference Map (p) 的尺寸与生成的 AF 不匹配，跳过均值计算。');
        end
    end
    %%%%%%% [新增结束] %%%%%
    
    % 2D切片
    figure;
    % subplot(2,2,1);
    if normalize_flag==1
        imagesc(doppler_range, delay_range, 10*log10(AF));
        xlabel('归一化多普勒频率');
    else
        imagesc(doppler_Hz, delay_range, 10*log10(AF));
        xlabel('多普勒频率(Hz)');
    end
    ylabel('时延(sample)');
    title('模糊函数');
    clim([-40 20]);
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
    
    figure;
    subplot(2,2,2);
    plot(delay_range, 10*log10(AF(:, round(Nv/2))), 'LineWidth', 1);
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
        plot(doppler_range, 10*log10(AF(23, :)), 'LineWidth', 2);
        xlabel('归一化多普勒频率');
    else
        plot(doppler_Hz, 10*log10(AF(23, :)), 'LineWidth', 2);
        xlabel('多普勒频率(Hz)');
    end
    ylabel('幅度(dB)');
    title('时延=50切面');
    grid on;
end