clc;clear;
% close all;
rng default;

sonar_signal_simulation('CW', 20, 0, 10);
sonar_signal_simulation('LFM', 20, 0, 10);
sonar_signal_simulation('Shaping', 20, 0, 10);
sonar_signal_simulation('Random', 20, 0, 10);

function sonar_signal_simulation(waveform_type, SNR, TSR, RSR)
    % 声呐信号仿真与处理（支持波形切换）
    % 输入参数：
    %   waveform_type - 波形类型：'CW'（连续波）或 'LFM'（线性调频）
    %   SNR - 信噪比(dB)
    %   TSR - 目标-混响比(dB)
    %   RSR - 混响-噪声比(dB)
    if nargin < 4
        RSR = 15; if nargin < 3
            TSR = 10; if nargin < 2
                SNR = 20; if nargin < 1
                    waveform_type = 'LFM';
                end
            end
        end
    end
    
    % 声呐信号仿真与处理
    % 仿真参数设置
    fs = 100e3;             % 采样频率 (Hz)
    pulse_duration = 0.1;   % 发s射脉冲持续时间 (s)
    analysis_duration = 0.5; % 总分析时长 (s)
    signal_len = pulse_duration * fs;
    t_pulse = 0:1/fs:pulse_duration-1/fs; % 发射信号时间向量
    t_analysis = 0:1/fs:analysis_duration-1/fs; % 分析周期时间向量
    
    c = 1500;               % 声速 (m/s)
    
    % 发射信号 (默认为线性调频信号)
    switch (waveform_type)
        case 'CW'
            % 连续波信号
            disp('使用CW(连续波)信号');
            fc = 15e3;              % 中心频率 (Hz)
            tx_signal = exp(1i*2*pi*fc*t_pulse);    
        case 'LFM'
            % 线性调频信号
            disp('使用LFM(线性调频)信号');
            fc = 15e3;              % 中心频率 (Hz)
            Bw = 1e3;        % 带宽 (Hz)
            tx_signal = exp(1i*pi*(Bw/pulse_duration)*t_pulse.^2) .* exp(1i*2*pi*fc*t_pulse); % 1xN complex double
        case 'Shaping'
            disp('使用AF Shaping信号');
            fc = 15e3; % 统一的中心频率
            data = load("100_100_1e4.mat", "s");
            s_generate = data.s;
            s = resample(s_generate, signal_len, length(s_generate));
            tx_signal = (s .* exp(1i*2*pi*fc*t_pulse.'))';  % t_pulse是(0:N-1)/fs
        case 'Random'
            disp('使用Random信号');
            fc = 15e3; % 统一的中心频率
            s = exp(1j * 2*pi * rand(signal_len,1));
            tx_signal = (s .* exp(1i*2*pi*fc*t_pulse.'))';
        % -------------------------------- TODO --------------------------------
        % case 'AF_Shaping'
        % 以CW、LFM为基础进行波形优化后进行比对
        %
        % 
        % ------------------------------------------------------------------------
        otherwise
            error('不支持的波形类型，请输入CW或LFM');
    end
    
    % 目标参数
    target_range = 50;      % 目标距离 (m)
    target_velocity = 3;     % 目标速度 (m/s), 正值为远离
    
    % 生成目标回波
    target_delay = 2*target_range/c;                 % 双程延迟 (s)
    target_doppler = 2*target_velocity*fc/c;         % 多普勒频移 (Hz)
    fprintf('Target Delay (s): %.1f\n', target_delay);
    fprintf('Target Doppler (Hz): %.2f\n', target_doppler);
    target_echo = generate_target_echo(tx_signal, t_pulse, t_analysis, target_delay, target_doppler, fs);
    
    % 生成混响
    N_r = 1000;             % 散射体数量
    reverberation = generate_reverberation(tx_signal, t_pulse, t_analysis, N_r, fc, c, fs);
    
    % 调整信号强度
    [target_echo, reverberation, noise] = adjust_signal_levels(...
        target_echo, reverberation, SNR, TSR, RSR);
    
    % 合成接收信号
    received_signal = target_echo + reverberation + noise;
    
    plot_waveforms(t_pulse, t_analysis, tx_signal, target_echo, reverberation, received_signal, fs, fc, waveform_type);

    % 匹配滤波处理
    [doppler_time_result, doppler_axis, time_axis] = matched_filter_processing(tx_signal, received_signal, fs, c);
    
    % 绘制结果
    plot_results(doppler_time_result, doppler_axis, time_axis, target_velocity, target_range, fc);

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
function reverberation = generate_reverberation(tx_signal, t_pulse, t_analysis, N_r, fc, c, fs)
    % 初始化混响信号
    reverberation = zeros(size(t_analysis));
    reverberation_no_decay = zeros(size(t_analysis));
    
    % 生成散射体参数
    max_range = c * (t_analysis(end) - t_pulse(end))/2;
    ranges = rand(1, N_r) * max_range;   % 随机距离
    delays = 2 * ranges / c;             % 双程延迟
    
    % 双指数分布的多普勒频移
    beta = 0.5; % Lower is Thinner
    dopplers = sign(rand(1, N_r)-0.5) .* exprnd(beta, 1, N_r) * fc/c;

    % ---------------- 测试代码 ----------------
    % % Doppler 分布可视化
    % figure;
    % histogram(dopplers, 100, 'Normalization', 'pdf'); % 使用归一化的直方图（表示概率密度）
    % title('Doppler Shift Distribution (Double-sided Exponential)');
    % xlabel('Doppler Shift (Hz)');
    % ylabel('Probability Density');
    % legend('Simulated Doppler Shifts');
    % grid on;
    % 可视化：混响点 Doppler-Range 散点图
    % figure;
    % scatter(ranges, dopplers, 20, 'filled');  % 每个点代表一个散射体
    % xlabel('Range (m)');
    % ylabel('Doppler Shift (Hz)');
    % title('Reverberation Scatterers: Range vs. Doppler');
    % grid on;
    % % 画出 attenuation vs. range
    % attenuations = 1 ./ log(ranges + exp(1));
    % figure;
    % plot(ranges, attenuations, '.');
    % xlabel('Range (m)');
    % ylabel('Attenuation Factor');
    % title('Attenuation vs. Range');
    % grid on;

    % ---------------- 测试代码 ----------------end
     

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

%% 噪声生成函数
function noise = generate_noise(signal, SNR)
    signal_power = mean(abs(signal).^2);
    noise_power = signal_power / (10^(SNR/10));
    noise = sqrt(noise_power/2) * (randn(size(signal)) + 1i*randn(size(signal)));
end

%% 信号强度调整
function [target_echo, reverberation, noise] = adjust_signal_levels(...
    target_echo, reverberation, SNR, TSR, RSR)
    
    % 计算基础功率
    P_target = mean(abs(target_echo).^2);
    P_reverb = mean(abs(reverberation).^2);
    
    % 根据TSR调整目标回波强度
    target_echo = target_echo * 10^(TSR/20) * sqrt(P_reverb/P_target);
    P_target = mean(abs(target_echo).^2);
    
    % 根据RSR调整混响强度
    reverberation = reverberation * 10^(RSR/20);
    P_reverb = mean(abs(reverberation).^2);
    
    % 根据SNR生成噪声
    noise_power = (P_target + P_reverb) / (10^(SNR/10));
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
    doppler_bins = 128;                  % 多普勒维度的分辨率
    doppler_range = 128;                 % 多普勒搜索范围 (Hz)
    
    % 时间延迟匹配滤波
    matched_filter = conj(fliplr(tx_signal));
    % time_compressed = conv(rx_signal, matched_filter, 'same');
    time_compressed = conv(rx_signal, matched_filter, 'full');
    
    % 创建多普勒频移矩阵
    doppler_axis = linspace(-doppler_range, doppler_range, doppler_bins);
    doppler_time_result = zeros(doppler_bins, length(time_compressed));
    
    % 对每个多普勒频移进行补偿和匹配滤波
    for k = 1:doppler_bins
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
function plot_results(doppler_time_result, doppler_axis, time_axis, target_velocity, target_range, fc)
    figure;
    
    % 绘制多普勒-时间(距离)图
    imagesc(time_axis, doppler_axis, 20*log10(abs(doppler_time_result)));
    axis xy;
    xlabel('Range (m)');
    ylabel('Doppler Shift (Hz)');
    title('Matched Filter Output: Doppler vs Range');
    colorbar;
    colormap('jet');
    c = clim;          % 获取当前 color axis 的范围
    new_min = 90;      % 想要设定的下限
    clim([new_min c(2)]);
    
    % 标记目标位置
    hold on;
    target_doppler = 2*target_velocity*fc/1500; % 计算理论多普勒频移
    plot(target_range, target_doppler, 'wo', 'MarkerSize', 10, 'LineWidth', 2);
    legend('Target Location');
    

    % % 寻找所有局部最大值位置
    % data = 20*log10(abs(doppler_time_result));
    % % 设置峰值检测阈值（高于最大值的-6dB）
    % threshold = max(data(:)) - 0.1;
    % % 寻找所有超过阈值的峰值
    % [rows, cols] = find(data >= threshold);
    % % 标记所有检测到的峰值位置
    % if ~isempty(rows)
    %     for i = 1:length(rows)
    %         plot(time_axis(cols(i)), doppler_axis(rows(i)), 'wp', ...
    %             'MarkerSize', 12, 'LineWidth', 1.5, ...
    %             'MarkerFaceColor', 'w');
    %     end
    % end
    % % 标记理论目标位置
    % target_doppler = 2*target_velocity*20e3/1500; % 计算理论多普勒频移
    % h_target = plot(target_range, target_doppler, 'wx', ...
    %     'MarkerSize', 12, 'LineWidth', 2);
    % % 创建图例
    % if ~isempty(rows)
    %     legend([h_target, plot(nan, nan, 'wp', 'MarkerFaceColor', 'w')], ...
    %         {'Theoretical Target', 'Detected Peaks'}, ...
    %         'Location', 'best');
    % else
    %     legend(h_target, 'Theoretical Target', 'Location', 'best');
    % end
    % hold off;

end