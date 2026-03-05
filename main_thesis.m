clc;clear;
close all;
rng default;

% -------------------------------- 仿真参数设置 --------------------------------
SNR = 10;        %   SNR - 信噪比(dB)
TSR = -15;       %   TSR - 目标-混响比(dB)

fs = 100e3;
fc = 15e3;
c0 = 1500;
Tp = 0.1;       % 脉冲宽度
Bw = 1e3;
analysis_duration = 0.5;    % 发射间隔
Waveform_type = 'AFShaping';    % AFShaping or Random


% 探测全过程参数
num_pings = 10;               % 总发射次数 (Ping数)
PRI = analysis_duration;      % 脉冲重复周期 (0.5s)

% 混响散射点生成
N_r = 1e3;  % 混响散射点数量
[ranges, delays, reverb_dopplers] = generate_scatters(Tp, analysis_duration, N_r, fc, c0);

% % ==================== 目标参数和运动真值生成 ====================
target_range = 100;          % 初始距离 100 m
target_velocity = 2;        % 初始速度 2 m/s (对应多普勒频移 40 Hz)
target_accel = 1.0;          % 设定的加速度 1.0 m/s^2
ping_accel_start = 4;        % 设定从第几个 Ping 开始加速
v_max = 4;                  % 最大速度 4 m/s (对应 80 Hz)

true_target_range    = zeros(1, num_pings);
true_target_velocity = zeros(1, num_pings);
true_target_doppler  = zeros(1, num_pings);
target_delay = 2*target_range/c0;                 % 双程延迟 (s)
target_doppler = 2*target_velocity*fc/c0;         % 多普勒频移 (Hz)
for ping_idx = 1:num_pings
    %  设定当前 Ping 的目标速度 (模拟: 匀速 -> 加速 -> 匀速)
    if ping_idx == 1
        current_v = target_velocity;
    else
        % 获取上一时刻的速度
        last_v = true_target_velocity(ping_idx - 1);
        % 判断是否处于加速阶段
        if ping_idx >= ping_accel_start && last_v < v_max
            % 运动学公式: v = v0 + a * t (这里的 t 是脉冲重复周期 PRI)
            current_v = last_v + target_accel * PRI;
            
            % 防止加速超过设定的最大速度
            current_v = min(current_v, v_max);
        else
            % 匀速阶段（还没开始加速，或者已经达到最大速度）
            current_v = last_v;
        end
    end

    %  目标真实距离 (径向靠近声呐)
    if ping_idx == 1
        current_r = target_range;
    else % 当前距离 = 上一 Ping 距离 - 速度 * 脉冲重复周期
        current_r = true_target_range(ping_idx-1) - true_target_velocity(ping_idx-1) * PRI;
    end
    
    % 计算真实多普勒频移
    current_fd = 2 * current_v * fc / c0;
    
    true_target_range(ping_idx)    = current_r;
    true_target_velocity(ping_idx) = current_v;
    true_target_doppler(ping_idx)  = current_fd;
end
% 可视化目标的动态运动过程
figure('Name', 'Target Dynamic Scenario');
subplot(2,1,1);
plot(1:num_pings, true_target_velocity, '-bo', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Ping 序号'); ylabel('速度 (m/s)');
title('目标径向速度变化曲线');
grid on; ylim([1.5 4.5]);
subplot(2,1,2);
plot(1:num_pings, true_target_doppler, '-ro', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
xlabel('Ping 序号'); ylabel('多普勒频移 (Hz)');
title('目标多普勒频移变化曲线 (覆盖波形库范围)');
grid on; ylim([35 85]);
% 3. 距离变化曲线
figure;
plot(1:num_pings, true_target_range, '-go', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
xlabel('Ping 序号'); 
ylabel('距离 (m)');
title('目标径向距离变化曲线 (向声呐靠近)');
grid on;
% ==================== 目标参数和运动真值生成END ====================


%%
close all;
fprintf('---------------- 发射波形: %s ----------------\n', Waveform_type);

load('Waveform_Lib_40to80_step10.mat');
search_doppler_window = 15;     % 多普勒匹配滤波后的局部搜索窗口大小

% 用于记录估计值和波形调用记录
estimated_doppler = zeros(1, num_pings);    % 速度估计值
estimated_range = zeros(1, num_pings);      % 距离估计值
selected_waveform_idx = zeros(1, num_pings);    % 各个Ping选择的波形
target_peak_responses = zeros(1, num_pings);    % 各个Ping的峰值响应

t_pulse = 0:1/fs:Tp-1/fs; % 发射信号时间向量
t_analysis = 0:1/fs:analysis_duration-1/fs; % 分析周期时间向量

signal_len = Tp * fs;
up_sample_rate = fs / Bw;
N_len = signal_len  / up_sample_rate;
fprintf('N_len: %d\n', N_len);

for ping_idx = 1:num_pings
    fprintf('正在处理第 %d/%d 个 Ping...\n', ping_idx, num_pings);
    

    % ---------------- 波形选择 ----------------
    if ping_idx == 1
        % 初始时刻：假设已知目标的初始状态，计算初始预期多普勒
        expected_doppler = true_target_doppler(1); 
    else
        % 后续时刻：直接将上一帧在匹配滤波中找到的多普勒峰值作为当前预期
        expected_doppler = estimated_doppler(ping_idx - 1);
    end
    [~, best_match_idx] = min(abs(library_dopplers - expected_doppler));
    selected_waveform_idx(ping_idx) = best_match_idx;
    current_tx_doppler_label = library_dopplers(best_match_idx);
    if strcmp(Waveform_type, 'AFShaping')
        s_baseband = Waveform_Library{best_match_idx};
    elseif strcmp(Waveform_type, 'Random')
        s_baseband = exp(1j * 2*pi * rand(N_len,1)); 
    end

    s_up = resample(s_baseband, up_sample_rate, 1);
    tx_signal = s_up' .* exp(1i*2*pi*fc*t_pulse);

    % ---------------- 生成环境回波 ----------------
    current_true_range = true_target_range(ping_idx);
    current_true_doppler = true_target_doppler(ping_idx);
    current_true_delay = 2 * current_true_range / c0;

    % 生成目标回波
    target_echo = generate_target_echo(tx_signal, t_pulse, t_analysis, current_true_delay, current_true_doppler, fs);
    % 生成混响
    reverberation = generate_reverberation(tx_signal, t_pulse, t_analysis, N_r, ranges, delays, reverb_dopplers, fs);
    % 调整信噪比/混响比
    [target_echo, reverberation, noise] = adjust_signal_levels(target_echo, reverberation, SNR, TSR);
    % 合成最终接收信号
    rx_signal = target_echo + reverberation + noise;

    % ---------------- 单Ping多普勒匹配滤波 ----------------
    [doppler_time_result, doppler_axis, time_axis] = matched_filter_processing(tx_signal, rx_signal, fs, c0);

    % ---------------- 在局部多普勒波门内提取目标 ----------------
    % 找到预期多普勒对应的中心索引
    [~, center_dop_idx] = min(abs(doppler_axis - expected_doppler));
    
    % 计算多普勒分辨率和搜索范围(bins)
    doppler_resolution = doppler_axis(2) - doppler_axis(1);
    dop_idx_span = round(search_doppler_window / doppler_resolution);
    
    % 划定搜索索引边界 (防止越界)
    search_idx_start = max(1, center_dop_idx - dop_idx_span);
    search_idx_end   = min(length(doppler_axis), center_dop_idx + dop_idx_span);
    
% 截取局部搜索波门数据
    search_area = abs(doppler_time_result(search_idx_start:search_idx_end, :));
    
    % 在波门内寻找能量最大值作为当前 Ping 的目标测量值
    [max_val, linear_idx] = max(search_area(:));
    [local_dop_idx, local_range_idx] = ind2sub(size(search_area), linear_idx);
    
    % 映射回全局索引
    global_dop_idx = search_idx_start + local_dop_idx - 1;
    
    % 记录当前帧估计的 多普勒 和 距离
    estimated_doppler(ping_idx) = doppler_axis(global_dop_idx);
    estimated_range(ping_idx) = time_axis(local_range_idx); 
    target_peak_responses(ping_idx) = 20*log10(max_val);
    
    if strcmp(Waveform_type, 'AFShaping')
        % 打印当前 Ping 的状态信息
        fprintf('    -> 真值: [%5.1f m, %5.1f Hz] | 估计: [%5.1f m, %5.1f Hz] | 发射波形: %d Hz\n', ...
            current_true_range, current_true_doppler, ...
            estimated_range(ping_idx), estimated_doppler(ping_idx), ...
            current_tx_doppler_label);
    elseif strcmp(Waveform_type, 'Random')
        fprintf('    -> 真值: [%5.1f m, %5.1f Hz] | 估计: [%5.1f m, %5.1f Hz] | 发射波形: Fix Random\n', ...
            current_true_range, current_true_doppler, ...
            estimated_range(ping_idx), estimated_doppler(ping_idx));
    end

        
    % ---------------- 基础可视化当前 Ping 的结果 ----------------
    % 将当前 Ping 的所有真值、估计值、波门信息传入画图函数
    plot_results(doppler_time_result, doppler_axis, time_axis, ...
        current_true_range, current_true_doppler, ...
        estimated_range(ping_idx), estimated_doppler(ping_idx), ...
        expected_doppler, search_doppler_window, ...
        ping_idx, num2str(current_tx_doppler_label));
end
% sonar_signal_simulation('Shaping_Q', SNR, TSR, fc, fs, c0, Tp, analysis_duration, N_r, ranges, delays, dopplers);
%%

function sonar_signal_simulation(waveform_type, SNR, TSR, ...
                                    fc, fs, c, pulse_duration, analysis_duration, ...
                                    N_r, ranges, delays, dopplers)
    % 声呐信号仿真与处理（支持波形切换）
    % 输入参数：
    %   waveform_type - 波形类型
    disp(['---------------- 使用 ', waveform_type, ' 信号 ----------------']);


    t_pulse = 0:1/fs:pulse_duration-1/fs; % 发射信号时间向量
    t_analysis = 0:1/fs:analysis_duration-1/fs; % 分析周期时间向量
    prf = 1 / pulse_duration;       % 脉冲重复频率（Hz）
    Bw = 1e3;
    fprintf('Bandwidth  %.2f kHz \t\t Range Resolution (m): %.3f \n', Bw/1e3, c/Bw/2);
    signal_len = pulse_duration * fs;
    up_sample_rate = fs / Bw;
    N_len = signal_len  / up_sample_rate;
    fprintf('N_len: %d\n', N_len);

    % 发射信号 (默认为线性调频信号)
    switch (waveform_type)
        case 'Shaping_Q'
            data = load("100_100_5e3_local_ISLQ_60Hz.mat", "s");
            s_generate = data.s;    % 100*1
            s = resample(s_generate, up_sample_rate, 1); % 10000*1
            tx_signal = s' .* exp(1i*2*pi*fc*t_pulse);  % t_pulse是(0:N-1)/fs 1*10000
        case 'Random'
            s_generate = exp(1j * 2*pi * rand(N_len,1));
            s = resample(s_generate, up_sample_rate, 1);
            tx_signal = s' .* exp(1i*2*pi*fc*t_pulse);  % t_pulse是(0:N-1)/fs 1*10000
            fprintf('Random Signal Bandwidth (kHz): %.2f\n', fs/up_sample_rate/1e3);
        otherwise
            error('不支持的波形类型');
    end
    plot_AF(tx_signal, fs, prf, waveform_type);
    % PAR
    instantaneous_power = abs(tx_signal).^2;   % 瞬时功率
    peak_power = max(instantaneous_power);      % 峰值功率
    average_power = mean(instantaneous_power);  % 平均功率
    fprintf('峰均功率比 (PAPR) = %.2f \n', peak_power / average_power);
    
    % 生成目标回波
    target_echo = generate_target_echo(tx_signal, t_pulse, t_analysis, target_delay, target_doppler, fs);

    % 生成混响
    reverberation = generate_reverberation(tx_signal, t_pulse, t_analysis, N_r, ranges, delays, dopplers, fs);

    % 调整信号强度
    [target_echo, reverberation, noise] = adjust_signal_levels(target_echo, reverberation, SNR, TSR);

    % 合成接收信号
    received_signal = target_echo + reverberation + noise;

    % 可视化时域波形
    plot_waveforms(t_pulse, t_analysis, tx_signal, target_echo, reverberation, received_signal, fs, fc, waveform_type);

    % 匹配滤波处理
    [doppler_time_result, doppler_axis, time_axis] = matched_filter_processing(tx_signal, received_signal, fs, c);

    % 绘制结果
    plot_results(doppler_time_result, doppler_axis, time_axis, target_velocity, target_range, fc, waveform_type);

    % ----------------------- Test code start    ----------------------- 
    % ----------------------- Test code end      -----------------------
end


%%
function [ranges, delays, dopplers] = generate_scatters(pulse_duration, analysis_duration, N_r, fc, c)
    max_range = c * (analysis_duration - pulse_duration)/2;
    ranges = (0.03 + (1 - 0.03) * rand(1, N_r)) * max_range;   % 随机距离
    delays = 2 * ranges / c;             % 双程延迟
    
    % 双指数分布的多普勒频移
    s = 20;   % s \in [6, 20], s 越大，频移越轻微
    mu = log(10^(s/20));
    reverberation_knots = sign(rand(1, N_r)-0.5) .* exprnd(1/mu, 1, N_r);
    dopplers = 2 * reverberation_knots * 0.514444 * fc / c;
    % ---------------- 测试代码 ----------------
    % Doppler 分布可视化
    figure;
    histogram(dopplers, 100, 'Normalization', 'pdf'); % 使用归一化的直方图（表示概率密度）
    title('Doppler Shift Distribution (Double-sided Exponential)');
    xlabel('Doppler Shift (Hz)');
    ylabel('Probability Density');
    legend('Simulated Doppler Shifts');
    grid on;
    % 可视化：混响点 Doppler-Range 散点图
    figure;
    scatter(ranges, dopplers, 20, 'filled');  % 每个点代表一个散射体
    xlabel('Range (m)');
    ylabel('Doppler Shift (Hz)');
    title('Reverberation Scatterers: Range vs. Doppler');
    grid on;
    % 画出 attenuation vs. range
    attenuations = 1 ./ log(ranges + exp(1));
    figure;
    plot(ranges, attenuations, '.');
    xlabel('Range (m)');
    ylabel('Attenuation Factor');
    title('Attenuation vs. Range');
    grid on;

    % ---------------- 测试代码 ----------------end
end

%% 目标回波生成函数
function target_echo = generate_target_echo(tx_signal, t_pulse, t_analysis, delay, doppler, fs)
    % 初始化全零信号
    target_echo = zeros(size(t_analysis));
    
    % 计算延迟样本数
    delay_samples = round(delay * fs);
    pulse_samples = length(t_pulse);
    
    % 确保回波完全落在分析周期内
    if delay_samples + pulse_samples > length(t_analysis)
        warning('目标回波超出分析周期，结果将被截断');
        pulse_samples = length(t_analysis) - delay_samples;
        if pulse_samples <= 0
            error('目标延迟过长，无法在分析周期内观察到回波');
        end
    end
    
    % 生成带多普勒频移的回波
    t_echo = t_pulse(1:pulse_samples);
    echo_signal = tx_signal(1:pulse_samples) .* exp(1i*2*pi*doppler*t_echo);
    
    % 将回波放入正确位置
    target_echo(delay_samples+1:delay_samples+pulse_samples) = echo_signal;
end

%% 混响生成函数
function reverberation = generate_reverberation(tx_signal, t_pulse, t_analysis, N_r, ranges, delays, dopplers, fs)
    % 初始化混响信号
    reverberation = zeros(size(t_analysis));
    reverberation_no_decay = zeros(size(t_analysis));

    % 对每个散射体生成回波并叠加
    for i = 1:N_r
        echo = generate_target_echo(tx_signal, t_pulse, t_analysis, delays(i), dopplers(i), fs);
        attenuation = 1 / log(ranges(i) + exp(1));      % attenuations = 1 ./ ((ranges + 1e-6) .^ attenuation_n);
        reverberation = reverberation + attenuation * echo;
        reverberation_no_decay = reverberation_no_decay + echo;
    end

    % ---------------- 测试代码 ----------------
    % t_ms = t_analysis * 1e3;
    % figure;
    % plot(t_ms, real(reverberation_no_decay), 'b', 'DisplayName', 'No Distance Attenuation'); hold on;
    % plot(t_ms, real(reverberation), 'r', 'DisplayName', 'With Distance Attenuation');
    % xlabel('Time (ms)');
    % ylabel('Amplitude');
    % title('Reverberation Comparison (With vs Without Distance Attenuation)');
    % legend('Location', 'best');
    % grid on;
    % ---------------- 测试代码 ----------------end
end


%% 信号强度调整
function [target_echo, reverberation, noise] = adjust_signal_levels(...
    target_echo, reverberation, SNR, TSR)
    
    % 计算基础功率
    P_target = mean(abs(target_echo).^2);
    P_reverb = mean(abs(reverberation).^2);
    
    % 根据TSR调整目标回波强度
    target_echo = target_echo * 10^(TSR/20) * sqrt(P_reverb/P_target);
    P_target = mean(abs(target_echo).^2);
    
    % 根据RSR调整混响强度
    % reverberation = reverberation * 10^(RSR/20);
    % P_reverb = mean(abs(reverberation).^2);
    
    % 根据SNR生成噪声
    noise_power = (P_target) / (10^(SNR/10));
    noise = sqrt(noise_power/2) * (randn(size(target_echo)) + 1i*randn(size(target_echo)));
end

%% 波形可视化函数
function plot_waveforms(t_pulse, t_analysis, tx_signal, target_echo, reverberation, received_signal, fs, fc, waveform_type)
    % 时域绘图
    figure('Name', ['Waveform Analysis - ', waveform_type]);
    
    % 发射信号时域
    subplot(4,1,1);
    plot(t_pulse*1e3, real(tx_signal));
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('Transmitted Signal (Time Domain)');
    grid on;
    
    % 目标回波时域
    subplot(4,1,2);
    plot(t_analysis*1e3, real(target_echo));
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('Target Echo (Time Domain)');
    % ylim([-400 400]);
    grid on;
    
    % 混响时域
    subplot(4,1,3);
    plot(t_analysis*1e3, real(reverberation));
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('Reverberation (Time Domain)');
    % ylim([-400 400]);
    grid on;
    
    % 接收信号时域
    subplot(4,1,4);
    plot(t_analysis*1e3, real(received_signal));
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('Received Signal (Time Domain)');
    % ylim([-400 400]);
    grid on;
    
    % 添加总标题
    sgtitle(['Signal Components Analysis - ', waveform_type, ' Waveform'], 'FontSize', 14);
end
%% 匹配滤波处理函数
function [doppler_time_result, doppler_axis, time_axis] = matched_filter_processing(tx_signal, rx_signal, fs, c)
    % 参数设置
    doppler_bins = 256;                  % 多普勒维度的分辨率
    doppler_range = 256;                 % 多普勒搜索范围 (Hz)
    
    % 时间延迟匹配滤波
    matched_filter = conj(fliplr(tx_signal));
    % time_compressed = conv(rx_signal, matched_filter, 'same');
    time_compressed = conv(rx_signal, matched_filter, 'full');
    
    % 创建多普勒频移矩阵
    doppler_axis = linspace(-doppler_range, doppler_range, doppler_bins+1);
    doppler_time_result = zeros(doppler_bins+1, length(time_compressed));
    
    % 对每个多普勒频移进行补偿和匹配滤波
    for k = 1:doppler_bins+1
        doppler_compensation = exp(-1i*2*pi*doppler_axis(k)*(0:length(rx_signal)-1)/fs);
        compensated_signal = rx_signal .* doppler_compensation(1:length(rx_signal));
        % doppler_time_result(k, :) = abs(conv(compensated_signal, matched_filter, 'same'));
        doppler_time_result(k, :) = abs(conv(compensated_signal, matched_filter, 'full'));
    end
    
    M = length(matched_filter);                     % 匹配滤波器长度
    N = length(rx_signal);                          % 接收信号长度
    L_full = N + M - 1;  % full 卷积输出长度
    valid_range = M : (L_full - M + 1);  % 有效匹配结果索引
    doppler_time_result = doppler_time_result(:, valid_range);

    time_axis = (0:N - M) / fs * c / 2;
end

%% 结果绘制函数
function plot_AF(sig, fs, prf, waveform_type)
    % [afmag_NAF, delay_NAF, doppler] = ambgfun(tx_signal, fs, prf);
    [afmag, delay, doppler] = ambgfun(sig, fs, prf);
    
    % 3D AF
    % figure();
    % mesh(delay,doppler,afmag);
    % colorbar;
    % colormap('jet');
    % set(gcf,'color','w');
    
    % 2D AF
    % figure();
    % imagesc(delay*1e3,doppler,mag2db(afmag))
    % xlabel('Delay \tau (ms)')
    % ylabel('Doppler f_d (Hz)')
    % colorbar;
    % colormap('jet');
    % clim([-50 0])
    % ylim([-500 500])
    

    % Q function
    % figure();
    % Q_func = sum(abs(afmag).^2, 2);   % 对每个多普勒下所有的时延积分，得到 Q 函数
    % Q_func = Q_func / max(Q_func);    % 归一化
    % Q_dB = 10 * log10(Q_func);        % 转为 dB
    % plot(doppler, Q_dB);
    % title('Q function dB');
    % xlim([-500 500])
    
    % 匹配滤波结果图（0多普勒切面）
    figure();
    [afmag, delay] = ambgfun(sig, fs, prf,'Cut','Doppler');
    afmag_db = mag2db(afmag);
    plot(delay*1500, afmag_db, 'LineWidth', 1);
    % xlabel('Delay \tau (ms)');
    xlabel('Range (m)')
    ylabel('Response (dB)');
    title(['Matched Filter Zero-Doppler Cut - ', waveform_type], 'Interpreter', 'none');
    ylim([-60 0])
    % xlim([-30 30])
    grid on;
end

%% 基础结果绘制函数 (单 Ping 可视化)
function plot_results(doppler_time_result, doppler_axis, time_axis, ...
                      true_r, true_fd, est_r, est_fd, exp_fd, dop_window, ...
                      ping_idx, waveform_type)
                  
    % 为每一 Ping 创建一个独立的新窗口
    figure('Name', ['Ping ', num2str(ping_idx), ' Result']);
    
    % 绘制多普勒-时间(距离)热力图
    imagesc(time_axis, doppler_axis, 20*log10(abs(doppler_time_result)));
    axis xy;
    xlabel('Range (m)');
    ylabel('Doppler Shift (Hz)');
    title(['Ping ', num2str(ping_idx), ' | Tx Waveform: ', waveform_type, 'Hz'], 'Interpreter', 'none');
    colorbar;
    colormap('jet');
    
    % 调整对比度，只显示最高 15dB 的动态范围，过滤掉大部分底噪
    c = clim;
    clim([c(2)-10 c(2)]);
    
    hold on;
    
    % 1. 标出目标的真实位置（红色空心大圆圈，Size=14）
    plot(true_r, true_fd, 'ro', 'MarkerSize', 18, 'LineWidth', 2, 'DisplayName', 'Target True Pos');
    
    % 2. 标出波门内搜索到的峰值估计位置（黑色空心小正方形，Size=8）
    plot(est_r, est_fd, 's', 'Color', [0.9, 0.9, 0.9], 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Estimated Peak');
    
    % 3. 画出多普勒搜索波门的上下界限（白色虚线）
    gate_upper = exp_fd + dop_window;
    gate_lower = exp_fd - dop_window;
    yline(gate_upper, 'w--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(gate_lower, 'w--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    legend('Location', 'northeast');
    
    % 限制坐标轴显示范围，放大目标附近的区域，便于观察波形凹陷和波门
    % xlim([max(0, true_r - 40), true_r + 40]);
    xlim([0 150]);
    ylim([0, 200]);
    grid on;
    hold off;
end
% ----------------------- Test code start    ----------------------- 
function beamwidth = calculate_3dB_beamwidth(range_axis, response)
    max_response = max(response);
    half_power = max_response - 3;
    
    % 找到响应高于半功率点的所有点
    above_half_power = response >= half_power;
    
    if any(above_half_power)
        % 找到第一个和最后一个高于半功率点的索引
        first_idx = find(above_half_power, 1, 'first');
        last_idx = find(above_half_power, 1, 'last');
        
        % 计算波束宽度
        beamwidth = range_axis(last_idx) - range_axis(first_idx);
    else
        beamwidth = NaN;
    end
end
% ----------------------- Test code end      -----------------------