clc;clear;
close all;
rng default;

SNR = 10;
TSR = -14;
RSR = 0;

tic
sonar_signal_simulation('Test', SNR, TSR, RSR);
%%
sonar_signal_simulation('Shaping_Q', SNR, TSR, RSR);
%%
sonar_signal_simulation('Shaping', SNR, TSR, RSR);
%%
sonar_signal_simulation('Shaping_ISL', SNR, TSR, RSR);
%%
sonar_signal_simulation('SFM', SNR, TSR, RSR);
%%
% sonar_signal_simulation('CW', SNR, TSR, RSR);
%%
% sonar_signal_simulation('LFM', SNR, TSR, RSR);
%%
% sonar_signal_simulation('GC', SNR, TSR, RSR);
%%
sonar_signal_simulation('GSFM', SNR, TSR, RSR);
%%
sonar_signal_simulation('Random', SNR, TSR, RSR);
%%
sonar_signal_simulation('BPSK', SNR, TSR, RSR);
%%
toc

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
    prf = 1 / pulse_duration;       % 脉冲重复频率（Hz）
    analysis_duration = 0.5; % 总分析时长 (s)
    signal_len = pulse_duration * fs;
    t_pulse = 0:1/fs:pulse_duration-1/fs; % 发射信号时间向量
    t_analysis = 0:1/fs:analysis_duration-1/fs; % 分析周期时间向量
    N_len = 100;
    p = signal_len / N_len;
    c = 1500;               % 声速 (m/s)
    fc = 15e3; % 统一的中心频率
    Bw = 1e3;        % 带宽 (Hz)

    % 发射信号 (默认为线性调频信号)
    switch (waveform_type)
        case 'CW'
            % 连续波信号
            disp('使用CW(连续波)信号');
            tx_signal = exp(1i*2*pi*fc*t_pulse);    
        case 'LFM'
            % 线性调频信号
            disp('使用LFM(线性调频)信号');
            f0 = fc - Bw/2;
            tx_signal = exp(1i*pi*(Bw/pulse_duration)*t_pulse.^2) .* exp(1i*2*pi*f0*t_pulse); % 1xN complex double
        case 'BPSK'
            % BPSK信号
            disp('使用BPSK信号');
            code_length = 13;       % 码元长度(可根据需要调整)
            chip_duration = pulse_duration/code_length;  % 每个码元的持续时间
            binary_seq = '1111100110101';   % randi([0 1], 1, code_length);
            phase_seq = binary_seq * pi;    % 将二进制序列转换为相位(0->0°, 1->180°)
            % 创建时间向量
            t_chip = 0:1/fs:chip_duration-1/fs;
            tx_signal = [];
            % 生成BPSK信号
            for i = 1:code_length
                chip_signal = exp(1i*(2*pi*fc*t_chip + phase_seq(i)));
                tx_signal = [tx_signal chip_signal];
            end
            tx_signal = tx_signal(1:min(length(tx_signal), length(t_pulse)));
            if length(tx_signal) < length(t_pulse)
                tx_signal = [tx_signal zeros(1, length(t_pulse)-length(tx_signal))];
            end
        case 'SFM'
            disp('使用SFM信号');
            fm = 100;
            beta = Bw / fm / 2 - 1; % B = 2fm*(1+beta)
            sfm_finst = @(t) fc+beta*fm*cos(2*pi*fm*t);
            sfmwav = phased.CustomFMWaveform(...
                'SampleRate',fs, ...
                'PulseWidth',pulse_duration, ...
                'FrequencyModulation', sfm_finst, ...
                'PRF', prf);
            sfm_samples = sfmwav();
            tx_signal = sfm_samples.';
        % -------------------------------- TODO --------------------------------
        % case 'AF_Shaping'
        % 以CW、LFM为基础进行波形优化后进行比对
        case 'Shaping'
            disp('使用AF Shaping信号');
            data = load("100_100_5e3_local_ISLQ_60Hz.mat",'s', 'interference_map');
            s_generate = data.s;
            inference = data.interference_map;
            % plot_AF(s_generate, fs/p, prf);
            s = resample(s_generate, signal_len, length(s_generate));
            tx_signal = s' .* exp(1i*2*pi*fc*t_pulse);  % t_pulse是(0:N-1)/fs 1*10000
        case 'Shaping_Q'
            disp('使用Shaping_Q信号');
            data = load("100_100_5e3_local_Q_60Hz.mat", "s");
            s_generate = data.s;    % 100*1
            s = resample(s_generate, p, 1); % 10000*1
            tx_signal = s' .* exp(1i*2*pi*fc*t_pulse);  % t_pulse是(0:N-1)/fs 1*10000

            % 调制解调、升降采样的AF比对测试代码
            % plot_AF(s_generate, fs/p, prf);
            % plot_AF(s, fs, prf);
            % rx_signal = tx_signal .* exp(-1i*2*pi*fc*t_pulse); 
            % tx_signal_downsample = resample(rx_signal, 1, p);
            % figure;plot(real(s_generate));
            % figure;plot(real(s));
            % figure;plot(real(tx_signal));
            % figure;plot(real(rx_signal));
            % figure;plot(real(tx_signal_downsample));
            % plot_AF(tx_signal_downsample, fs/p, prf);
        case 'Shaping_ISL'
            disp('使用Shaping_ISL信号');
            data = load("100_100_5e3_local_ISL_60Hz.mat", "s");
            s_generate = data.s;
            s = resample(s_generate, signal_len, length(s_generate));
            tx_signal = s' .* exp(1i*2*pi*fc*t_pulse);  % t_pulse是(0:N-1)/fs 1*10000

        case 'GC'
            disp('使用GC信号');
            GC_N_components = 20;            % 谐波数量（频点数量）
            GC_f0 = fc - Bw/2;
            GC_s = 250;             % s 为起始的频率间隔(Hz)
            GC_r = 1.05;              % r 为一个稍大于1的比例系数
            % 计算GC频点
            GC_frq = zeros(1, GC_N_components);            % f_i = f_{i-1}+sr^{i-2}
            GC_frq(1,1) = GC_f0;
            max_freq = fc + Bw/2;
            for i = 2:GC_N_components
                next_f = GC_frq(i-1) + GC_s * GC_r^(i-2);
                if next_f >= max_freq
                    GC_N_components = i-1; % 更新实际分量数量
                    break;
                end
                GC_frq(1,i) = next_f;
            end
            GC_frq = GC_frq(1:GC_N_components);
            % Phase
            GC_phase = zeros(1, GC_N_components);
            method = 'Newman'; % 可选 'Narahashi' 或 'Newman'
            for i = 1:GC_N_components
                switch method
                    case 'Narahashi'
                        % Narahashi 公式 (φ_i = π(i-1)(i-2)/(N-1))
                        GC_phase(i) = pi * (i-1) .* (i-2) / (GC_N_components-1);

                    case 'Newman'
                        % Newman 公式 (φ_i = π(i-1)^2/N)
                        GC_phase(i) = pi * (i-1).^2 / GC_N_components;

                    otherwise
                        % 默认随机相位（保留原有逻辑）
                        GC_phase = 2*pi * rand(1, GC_N_components);
                end
            end

            gc_samples = zeros(1, signal_len);
            for n = 1:GC_N_components
                gc_samples = gc_samples + exp(1j*(2*pi*GC_frq(n)*t_pulse + GC_phase(n)));
            end
            tx_signal = gc_samples/max(abs(gc_samples));
        case 'GSFM'
            disp('使用GSFM信号');
            rho = 1.9;               % Variable exponent parameter
            alpha = 1e3;             % Frequency modulation term
            gsfm_finst = @(t)(fc + (Bw/2)*sin(2*pi*alpha*t.^rho));
            gsfmwav = phased.CustomFMWaveform('SampleRate',fs,...
                'PulseWidth',pulse_duration,'FrequencyModulation',gsfm_finst,'PRF',prf);
            gsfm_samples = gsfmwav();
            tx_signal = gsfm_samples.';
        % ------------------------------------------------------------------------
        case 'Random'
            disp('使用Random信号');
            s_generate = exp(1j * 2*pi * rand(100,1));
            s = resample(s_generate, p, 1);
            tx_signal = s' .* exp(1i*2*pi*fc*t_pulse);  % t_pulse是(0:N-1)/fs 1*10000
        case 'Test'
            if ~exist('matlab.mat', 'file')
                disp('File Test.mat does not exist. Exiting function.');
                return; % 如果文件不存在，直接返回
            end
            disp('使用Test信号');
            data = load("Test.mat", "s");
            s_generate = data.s;
            s = resample(s_generate, signal_len, length(s_generate));
            tx_signal = s' .* exp(1i*2*pi*fc*t_pulse);  % t_pulse是(0:N-1)/fs 1*10000
        otherwise
            error('不支持的波形类型，请输入CW或LFM');
    end
    % plot_AF(tx_signal, fs, prf);
    % PAR
    instantaneous_power = abs(tx_signal).^2;   % 瞬时功率
    peak_power = max(instantaneous_power);      % 峰值功率
    average_power = mean(instantaneous_power);  % 平均功率
    fprintf('峰均功率比 (PAPR) = %.2f \n', peak_power / average_power);
    % % 频谱图
    % figure;
    % f = (-signal_len/2:signal_len/2-1)*(fs/signal_len);
    % X_shifted = fftshift(fft(tx_signal));
    % X = abs(X_shifted) / max(abs(X_shifted));
    % plot(f, X);xlim([fc-Bw fc+Bw]);
    % title([waveform_type, ' Spectrum']);


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

    % 可视化时域波形
    % plot_waveforms(t_pulse, t_analysis, tx_signal, target_echo, reverberation, received_signal, fs, fc, waveform_type);

    % 匹配滤波处理
    [doppler_time_result, doppler_axis, time_axis] = matched_filter_processing(tx_signal, received_signal, fs, c);

    % 绘制结果
    plot_results(doppler_time_result, doppler_axis, time_axis, target_velocity, target_range, fc, waveform_type);

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
    ranges = (0.03 + (1 - 0.03) * rand(1, N_r)) * max_range;   % 随机距离
    delays = 2 * ranges / c;             % 双程延迟
    
    % 双指数分布的多普勒频移
    s = 20;   % s \in [6, 20], s 越大，频移越轻微
    mu = log(10^(s/20));
    reverberation_knots = sign(rand(1, N_r)-0.5) .* exprnd(1/mu, 1, N_r);
    dopplers = 2 * reverberation_knots * 0.514444 * fc / c;
    
    % ---------------- 测试代码 ----------------
    % Doppler 分布可视化
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
    target_echo, reverberation, SNR, TSR, ~)
    
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
function plot_results(doppler_time_result, doppler_axis, time_axis, target_velocity, target_range, fc, waveform_type)
    figure;
    
    % 绘制多普勒-时间(距离)图
    imagesc(time_axis, doppler_axis, 20*log10(abs(doppler_time_result)));
    axis xy;
    xlabel('Range (m)');
    ylabel('Doppler Shift (Hz)');
    title(['Matched Filter Output: Doppler vs Range - ',waveform_type], 'Interpreter', 'none');
    colorbar;
    colormap('jet');
    c = clim;          % 获取当前 color axis 的范围
    new_min = 90;      % 想要设定的下限
    clim([c(2)-8 c(2)]);
    xlim([0 150])
    
    % % 标记目标位置
    % hold on;
    % target_doppler = 2*target_velocity*fc/1500; % 计算理论多普勒频移
    % plot(target_range, target_doppler, 'wo', 'MarkerSize', 10, 'LineWidth', 2);
    % legend('Target Location');
    

    % 寻找所有局部最大值位置
    data = 20*log10(abs(doppler_time_result));
    % 设置峰值检测阈值（高于最大值的-6dB）
    threshold = max(data(:)) - 0.1;
    % 寻找所有超过阈值的峰值
    [rows, cols] = find(data >= threshold);
    % 标记所有检测到的峰值位置
    hold on;
    if ~isempty(rows)
        for i = 1:length(rows)
            plot(time_axis(cols(i)), doppler_axis(rows(i)), 'ro', ...
                'MarkerSize', 12, 'LineWidth', 1.5);
        end
    end
    % 标记理论目标位置
    target_doppler = 2*target_velocity*fc/1500; % 计算理论多普勒频移
    h_target = plot(target_range, target_doppler, 'wo', ...
        'MarkerSize', 10, 'LineWidth', 2);
    % 创建图例
    if ~isempty(rows)
        legend([h_target, plot(nan, nan, 'ro', 'MarkerFaceColor', 'w')], ...
            {'Theoretical Target', 'Detected Peaks'}, ...
            'Location', 'best');
    else
        legend(h_target, 'Theoretical Target', 'Location', 'best');
    end
    hold off;

end

function plot_AF(sig, fs, prf)
    % [afmag_NAF, delay_NAF, doppler] = ambgfun(tx_signal, fs, prf);
    c = 1500;
    fc = 15e3;
    [afmag, delay, doppler] = ambgfun(sig, fs, prf);
    
    % 3D AF
    % figure();
    % mesh(delay,doppler,afmag);
    % colorbar;
    % colormap('jet');
    % set(gcf,'color','w');
    
    % 2D AF
    figure();
    imagesc(delay*1e3,doppler,mag2db(afmag))
    xlabel('Delay \tau (ms)')
    ylabel('Doppler f_d (Hz)')
    colorbar;
    colormap('jet');
    clim([-50 0])
    ylim([-500 500])
    

    % Q function
    figure();
    Q_func = sum(abs(afmag).^2, 2);   % 对每个多普勒下所有的时延积分，得到 Q 函数
    Q_func = Q_func / max(Q_func);    % 归一化
    Q_dB = 10 * log10(Q_func);        % 转为 dB
    plot(doppler, Q_dB);
    title('Q function dB');
    xlim([-500 500])
end