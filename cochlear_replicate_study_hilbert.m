% Sampling rate according to Shannon et al's study followed
fs = 10000; % Hz

% Filter bank parameters
nb = [1, 2, 4, 8, 16]; %Frequency Bands
% Spread across same total frequency range 0 to 4 kHz
f_low = 0; % Hz
f_high = 4000; % Hz

filter_order = 3;   %Mentioned in original study that signal was split into frequency bands using third-order elliptical filters
filter_ripple = 15; %Adjacent filters overlapping at the point where the output from each filter was 15 dB down from the level in the pass-band

% Pre-emphasis filter constructed using butter (applied according to original study)
pre_emphasis_fc = 1200; % Hz of low pass filter
[b_pre,a_pre] = butter(1, pre_emphasis_fc/(fs/2), 'high');

% Bandpass filter bank and Hilbert transform envelope extraction
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
    
    % Apply filter bank and Hilbert transform envelope extraction
    y = filter(filter_bank(n).b, filter_bank(n).a, noise_carrier);
    env = abs(hilbert(y));
    
    % Combine envelopes and filter with lowpass filter to remove high
    % frequency components
    y_env_lp = filter(b_lp, a_lp, env);
    
    % Amplify to achieve desired output level (75 dBA for continuous speech
    % according to original study)
    y_env_lp = y_env_lp * 10^(75/20);
    
    % Plot spectrogram for each band
    figure;
    spectrogram(y_env_lp, 256, [], [], fs, 'yaxis');
    title(['Filter Bank with Hilbert Transform with ', num2str(nb(n)), ' bands']);
end
