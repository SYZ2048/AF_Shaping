clc;clear;close all;

sonar_signal_simulation('CW');
% sonar_signal_simulation('LFM');


function sonar_signal_simulation(waveform_type)
    % 声呐信号仿真与处理（支持波形切换）
    % 输入参数：
    %   waveform_type - 波形类型：'CW'（连续波）或 'LFM'（线性调频）
    if nargin < 1
        waveform_type = 'LFM'; % 默认使用LFM信号
    end

    SNR = 10; %   SNR - 信噪比(dB)，控制噪声强度
    TSR = 20; %   TSR - 目标-混响比(Target-to-Reverberation Ratio, dB)
    RSR = 15; %   RSR - 混响-噪声比(Reverberation-to-Noise Ratio, dB)
    
    % 声呐信号仿真与处理
    % 仿真参数设置
    fs = 100e3;             % 采样频率 (Hz)
    T = 0.1;                % 信号持续时间 (s)
    t = 0:1/fs:T-1/fs;      % 时间向量
    fc = 20e3;              % 中心频率 (Hz)
    c = 1500;               % 声速 (m/s)
    
    % 发射信号 (默认为线性调频信号)
    switch upper(waveform_type)
        case 'CW'
            % 连续波信号
            disp('使用CW(连续波)信号');
            tx_signal = exp(1i*2*pi*fc*t);
            B = 1/T;        % CW的有效带宽近似为1/T
            
        case 'LFM'
            % 线性调频信号
            disp('使用LFM(线性调频)信号');
            B = 5e3;        % 带宽 (Hz)
            tx_signal = exp(1i*pi*(B/T)*t.^2) .* exp(1i*2*pi*fc*t);
            
        otherwise
            error('不支持的波形类型，请输入CW或LFM');
    end
    
    % 目标参数
    target_range = 50;      % 目标距离 (m)
    target_velocity = 3;     % 目标速度 (m/s), 正值为远离
    target_rcs = 1;          % 目标雷达截面积
    
    % 混响参数
    N_r = 1000;             % 散射体数量
    
    % 生成目标回波
    target_delay = 2*target_range/c;                 % 双程延迟 (s)
    target_doppler = 2*target_velocity*fc/c;         % 多普勒频移 (Hz)
    target_echo = generate_target_echo(tx_signal, t, target_delay, target_doppler, target_rcs, fs);
    
    % 生成混响
    reverberation = generate_reverberation(tx_signal, t, N_r, 1, fc, c, fs);
    
    % 计算信号功率
    P_tx = mean(abs(tx_signal).^2);
    P_target_base = mean(abs(target_echo).^2);
    P_reverb_base = mean(abs(reverberation).^2);
    target_echo = target_echo * 10^(TSR/20) * sqrt(P_reverb_base/P_target_base);
    reverberation = reverberation * 10^(RSR/20);

    % 生成高斯白噪声
    noise = generate_noise(tx_signal, SNR);
    
    % 合成接收信号
    received_signal = target_echo + reverberation + noise;
    
    plot_waveforms(t, tx_signal, target_echo, reverberation, received_signal, fs, fc, waveform_type);

    % 匹配滤波处理
    [doppler_time_result, doppler_axis, time_axis] = matched_filter_processing(tx_signal, received_signal, fs, fc, c);
    
    % 绘制结果
    plot_results(doppler_time_result, doppler_axis, time_axis, target_velocity, target_range);

end
%% 目标回波生成函数
function target_echo = generate_target_echo(tx_signal, t, delay, doppler, rcs, fs)
    % 计算延迟样本数
    delay_samples = round(delay * fs);
    
    % 确保延迟不超过信号长度
    max_delay = length(t) - 1;
    if delay_samples > max_delay
        delay_samples = max_delay;
    end
    
    % 创建延迟信号
    delayed_signal = [zeros(1, delay_samples), tx_signal(1:end-delay_samples)];
    
    % 应用多普勒频移
    t_shifted = t(1:end-delay_samples);
    doppler_shift = delayed_signal .* exp(1i*2*pi*doppler*t);
    
    % 应用RCS
    target_echo = sqrt(rcs) * doppler_shift;
end

%% 混响生成函数
function reverberation = generate_reverberation(tx_signal, t, N_r, level, fc, c, fs)
    % 初始化混响信号
    reverberation = zeros(size(tx_signal));
    
    % 生成散射体参数
    ranges = rand(1, N_r) * c * t(end)/2;           % 随机距离
    delays = 2 * ranges / c;                         % 双程延迟
    rcs_values = abs(randn(1, N_r)) * level;         % 随机RCS
    
    % 双指数分布的多普勒频移 (以0为中心)
    beta = 10; % 控制分布形状
    dopplers = sign(rand(1, N_r)-0.5) .* exprnd(beta, 1, N_r) * fc/c;
    
    % 对每个散射体生成回波并叠加
    for i = 1:N_r
        echo = generate_target_echo(tx_signal, t, delays(i), dopplers(i), rcs_values(i), fs);
        reverberation = reverberation + echo;
    end
end

%% 噪声生成函数
function noise = generate_noise(signal, SNR)
    signal_power = mean(abs(signal).^2);
    noise_power = signal_power / (10^(SNR/10));
    noise = sqrt(noise_power/2) * (randn(size(signal)) + 1i*randn(size(signal)));
end

%% 波形可视化函数
function plot_waveforms(t, tx_signal, target_echo, reverberation, received_signal, fs, fc, waveform_type)
    % 时域绘图
    figure('Name', ['Waveform Analysis - ', waveform_type], 'Position', [100, 100, 1200, 800]);
    
    % 发射信号时域
    subplot(3,2,1);
    plot(t*1e3, real(tx_signal));
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('Transmitted Signal (Time Domain)');
    xlim([0, max(t)*1e3]);
    grid on;
    
    % 目标回波时域
    subplot(3,2,3);
    plot(t*1e3, real(target_echo));
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('Target Echo (Time Domain)');
    % xlim([0, max(t)*1e3]);
    grid on;
    
    % 混响时域
    subplot(3,2,5);
    plot(t*1e3, real(reverberation));
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('Reverberation (Time Domain)');
    % xlim([0, max(t)*1e3]);
    grid on;
    
    % 接收信号时域
    subplot(3,2,2);
    plot(t*1e3, real(received_signal));
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('Received Signal (Time Domain)');
    % xlim([0, max(t)*1e3]);
    grid on;
    
    % 频域分析
    N = length(tx_signal);
    f = (-N/2:N/2-1)*(fs/N)/1e3; % 频率轴(kHz)
    
    % 发射信号频域
    subplot(3,2,4);
    tx_spectrum = abs(fftshift(fft(tx_signal)));
    plot(f, 20*log10(tx_spectrum/max(tx_spectrum)));
    xlabel('Frequency (kHz)');
    ylabel('Normalized Magnitude (dB)');
    title('Transmitted Signal (Frequency Domain)');
    xlim([fc/1e3-10, fc/1e3+10]);
    grid on;
    
    % 接收信号频域
    subplot(3,2,6);
    rx_spectrum = abs(fftshift(fft(received_signal)));
    plot(f, 20*log10(rx_spectrum/max(rx_spectrum)));
    xlabel('Frequency (kHz)');
    ylabel('Normalized Magnitude (dB)');
    title('Received Signal (Frequency Domain)');
    xlim([fc/1e3-10, fc/1e3+10]);
    grid on;
    
    % 添加总标题
    sgtitle(['Signal Components Analysis - ', waveform_type, ' Waveform'], 'FontSize', 14);
end
%% 匹配滤波处理函数
function [doppler_time_result, doppler_axis, time_axis] = matched_filter_processing(tx_signal, rx_signal, fs, fc, c)
    % 参数设置
    fft_length = 1024;
    doppler_bins = 128;                  % 多普勒维度的分辨率
    doppler_range = 100;                 % 多普勒搜索范围 (Hz)
    
    % 时间延迟匹配滤波
    matched_filter = conj(fliplr(tx_signal));
    time_compressed = conv(rx_signal, matched_filter, 'same');
    
    % 创建多普勒频移矩阵
    doppler_axis = linspace(-doppler_range, doppler_range, doppler_bins);
    doppler_time_result = zeros(doppler_bins, length(time_compressed));
    
    % 对每个多普勒频移进行补偿和匹配滤波
    for k = 1:doppler_bins
        doppler_compensation = exp(-1i*2*pi*doppler_axis(k)*(0:length(rx_signal)-1)/fs);
        compensated_signal = rx_signal .* doppler_compensation(1:length(rx_signal));
        doppler_time_result(k, :) = abs(conv(compensated_signal, matched_filter, 'same'));
    end
    
    % 时间轴 (转换为距离)
    time_axis = (0:length(time_compressed)-1)/fs * c/2;
end

%% 结果绘制函数
function plot_results(doppler_time_result, doppler_axis, time_axis, target_velocity, target_range)
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
    new_min = 115;      % 想要设定的下限
    clim([new_min c(2)]);
    
    % 标记目标位置
    hold on;
    target_doppler = 2*target_velocity*20e3/1500; % 计算理论多普勒频移
    plot(target_range, target_doppler, 'wx', 'MarkerSize', 10, 'LineWidth', 2);
    legend('Target Location');
    hold off;
    
    % 绘制3D视图
    % figure;
    % [X, Y] = meshgrid(time_axis, doppler_axis);
    % surf(X, Y, 20*log10(abs(doppler_time_result)), 'EdgeColor', 'none');
    % xlabel('Range (m)');
    % ylabel('Doppler Shift (Hz)');
    % zlabel('Amplitude (dB)');
    % title('3D View of Matched Filter Output');
    % colormap('jet');
    % c = caxis;          % 获取当前 color axis 的范围
    % new_min = 50;      % 想要设定的下限
    % caxis([new_min c(2)]);

end