% Sampling rate according to Shannon et al's study followed
fs = 10000; % Hz

% Filter bank parameters
nb = [1, 2, 4, 8, 16]; %Frequency Bands
% Spread across same total frequency range 0 to 4 kHz
f_low = 0; % Hz
f_high = 4000; % Hz

filter_order = 3;   %Mentioned in original study that signal was split into frequency bands using third-order elliptical filters
filter_ripple = 15; %Adjacent filters overlapping at the point where the output from each filter was 15 dB down from the level in the pass-band

env_filter_order = 2;   
env_filter_ripple = 6;  %Filter ripple of -6 dB per octave

% Pre-emphasis filter constructed using butter (applied according to original study)
pre_emphasis_fc = 1200; % Hz of low pass filter
[b_pre,a_pre] = butter(1, pre_emphasis_fc/(fs/2), 'high');

% Bandpass filter bank
for n = 1:length(nb)
    % Compute center frequencies for each band
    fc = linspace(f_low, f_high, nb(n)+2);
    fc = fc(2:end-1);
    
    wn = [fc(1)/(fs/2) fc(end)/(fs/2)];
    % Compute filter parameters for each band using third-order elliptical filters
    [b, a] = ellip(filter_order, filter_ripple, 40, wn, 'bandpass');
    
    % Save filter coefficients for each band
    filter_bank(n).b = b;
    filter_bank(n).a = a;
end

% Envelope extraction filters as specified
env_fc = [16, 160];
% For each cutoff frequency
for n = 1:length(env_fc)
    [env_b, env_a] = ellip(env_filter_order, env_filter_ripple, 40, env_fc(n)/(fs/2), 'low');

    % Save filter coefficients for each filter
    env_filter(n).b = env_b;
    env_filter(n).a = env_a;
end

% Generate noise carrier through using white noise (according to original study)
t = 0:1/fs:1;
noise_carrier = randn(length(t), 1);

% Apply filter bank and envelope extraction to noise carrier for each
% frequency band
for n = 1:length(nb)
    % Make sure that index in position 1 does not exceed bounds if nb <=1
    if nb(n) > 1
        % Apply filter bank
        y = zeros(length(noise_carrier), nb(n));
        for k = 1:nb(n)
            y(:, k) = filter(filter_bank(n).b, filter_bank(n).a, noise_carrier);
        end
    % Skip the for loop if not greater than 1
    else
        y = filter(filter_bank(n).b, filter_bank(n).a, noise_carrier);
    end
    % Apply envelope extraction
    for k = 1:length(env_fc)
        env = abs(hilbert(filter(env_filter(k).b, env_filter(k).a, y)));
    end
    
    % Combine envelopes and filter with lowpass filter to remove high
    % frequency components
    y_env = sum(env, 2);
    [b_lp, a_lp] = butter(6, 4000/(fs/2), 'low');
    y_env_lp = filter(b_lp, a_lp, y_env);
    
    % Amplify to achieve desired output level (75 dBA for continuous speech
    % according to original study)
    y_env_lp = y_env_lp * 10^(75/20);
    
    % Plot spectrogram for each band
    figure;
    spectrogram(y_env_lp, 256, [], [], fs, 'yaxis');
    title(['Filter Bank with ', num2str(nb(n)), ' bands']);
end