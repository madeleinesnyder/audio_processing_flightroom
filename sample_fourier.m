%% Fourier

load("sample_audio.mat")

mic = 1;

% Sample signal
Fs = fs; % Sampling frequency in Hz
t = 0:1/Fs:1; % Time vector (1 second)
signal = recbuf(:,mic); % Sine wave

% Perform Fourier Transform using FFT
n = length(recbuf); % Length of the signal
fTransform = fft(signal); % FFT of the signal
fAxis = (0:n-1)*(Fs/n); % Frequency range
magnitude = abs(fTransform/n); % Magnitude of the FFT

% Plot Fourier Transform
figure();
plot(fAxis, magnitude);
title('Magnitude of FFT of the Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
axis([0 100 0 max(magnitude)]); % Limiting x-axis for better visibility

% Generate Spectrogram
figure();
spectrogram(signal, 256, 250, 256, Fs, 'yaxis');



% Create a custom colormap (transition from black to red)
customCMap = [linspace(0, 1, 256)', zeros(256, 2)];  % Red channel increases, others stay at 0

% Apply the custom colormap
colormap(customCMap);

% Add color bar and title
colorbar;
title('Spectrogram with Custom Colormap');

