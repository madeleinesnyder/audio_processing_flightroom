
aze = cell(1,4)
for i=1:4
    aze{i} = (event_locs{i}-event_vec_start);
    aze{i} = aze{i}(aze{i}>=0);
    aze{i}(1) = 1;
    aze{i} = aze{i}(aze{i}<length(sum_whitenedSignal_allmics));
end

% Below from fb_pretty_sonogram code in FinchScope repo
overlap=2000;
tscale=2;
N=2048;
nfft=2^nextpow2(N);
low=2.9;
high=10;

t_=-N/2+1:N/2;
sigma=(tscale/1e3)*motu_Fs;
w = exp(-(t_/sigma).^2);
dw = -2*w.*(t_/(sigma^2));
    
% Calculate spectrogram
[S,F,T,P]=spectrogram(whitenedSignal_mm{1},w,overlap,nfft,motu_Fs);
[S2]=spectrogram(whitenedSignal_mm{1},dw,overlap,nfft,motu_Fs);
%[S,F,T,P]=spectrogram(whitenedSignal_mm{1},w,overlap,nfft,motu_Fs);
%[S2]=spectrogram(whitenedSignal_mm{1},dw,overlap,nfft,motu_Fs);
IMAGE=100*((abs(S)+abs(S2))/2);
IMAGE=log(IMAGE+2);
IMAGE(IMAGE>high)=high;
IMAGE(IMAGE<low)=low;
IMAGE=IMAGE-low;
IMAGE=IMAGE/(high-low);
IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
IMAGE_MODS = log(abs(IMAGE)+1e+2);
% Above from fb_pretty_sonogram code in FinchScope repo
IMAGE_MODS_SUM = sum(IMAGE_MODS);

spectrogram_ee_idxs = floor(aze{1}/(length(sum_whitenedSignal_allmics)/length(T)));
spectrogram_ee = NaN(1,length(T));
if spectrogram_ee_idxs(1) == 0; spectrogram_ee_idxs(1) = 1; end;
spectrogram_ee(spectrogram_ee_idxs) = 60000;

% Plot the spectrogram with the peaks highlighted
figure(); 
subplot(4,1,1); hold on; 
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gcf, 'Color', 'k'); 
colormap(hot)
imagesc(T,F,IMAGE_MODS); %set(gca,'YDir','normal');
scatter(T, spectrogram_ee',5,'y'); % Mark the peaks
ylabel('Frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('Spectrogram of whitened data','FontWeight','bold','Color','w');
axis tight;
hold off;

% Calculate spectrogram
[S,F,T,P]=spectrogram(whitenedSignal_mm{2},w,overlap,nfft,motu_Fs);
[S2]=spectrogram(whitenedSignal_mm{2},dw,overlap,nfft,motu_Fs);
%[S,F,T,P]=spectrogram(whitenedSignal_mm{1},w,overlap,nfft,motu_Fs);
%[S2]=spectrogram(whitenedSignal_mm{1},dw,overlap,nfft,motu_Fs);
IMAGE=100*((abs(S)+abs(S2))/2);
IMAGE=log(IMAGE+2);
IMAGE(IMAGE>high)=high;
IMAGE(IMAGE<low)=low;
IMAGE=IMAGE-low;
IMAGE=IMAGE/(high-low);
IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
IMAGE_MODS = log(abs(IMAGE)+1e+2);
% Above from fb_pretty_sonogram code in FinchScope repo
IMAGE_MODS_SUM = sum(IMAGE_MODS);

spectrogram_ee_idxs = floor(aze{2}/(length(sum_whitenedSignal_allmics)/length(T)));
spectrogram_ee = NaN(1,length(T));
if spectrogram_ee_idxs(1) == 0; spectrogram_ee_idxs(1) = 1; end;
spectrogram_ee(spectrogram_ee_idxs) = 60000;

% Plot the spectrogram with the peaks highlighted
subplot(4,1,2); hold on; 
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gcf, 'Color', 'k'); 
colormap(hot)
imagesc(T,F,IMAGE_MODS); %set(gca,'YDir','normal');
scatter(T, spectrogram_ee',5,'y'); % Mark the peaks
ylabel('Frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('Spectrogram of whitened data','FontWeight','bold','Color','w');
axis tight;
hold off;

% Calculate spectrogram
[S,F,T,P]=spectrogram(whitenedSignal_mm{3},w,overlap,nfft,motu_Fs);
[S2]=spectrogram(whitenedSignal_mm{3},dw,overlap,nfft,motu_Fs);
%[S,F,T,P]=spectrogram(whitenedSignal_mm{1},w,overlap,nfft,motu_Fs);
%[S2]=spectrogram(whitenedSignal_mm{1},dw,overlap,nfft,motu_Fs);
IMAGE=100*((abs(S)+abs(S2))/2);
IMAGE=log(IMAGE+2);
IMAGE(IMAGE>high)=high;
IMAGE(IMAGE<low)=low;
IMAGE=IMAGE-low;
IMAGE=IMAGE/(high-low);
IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
IMAGE_MODS = log(abs(IMAGE)+1e+2);
% Above from fb_pretty_sonogram code in FinchScope repo
IMAGE_MODS_SUM = sum(IMAGE_MODS);

spectrogram_ee_idxs = floor(aze{3}/(length(sum_whitenedSignal_allmics)/length(T)));
spectrogram_ee = NaN(1,length(T));
if spectrogram_ee_idxs(1) == 0; spectrogram_ee_idxs(1) = 1; end;
spectrogram_ee(spectrogram_ee_idxs) = 60000;

% Plot the spectrogram with the peaks highlighted
subplot(4,1,3); hold on; 
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gcf, 'Color', 'k'); 
colormap(hot)
imagesc(T,F,IMAGE_MODS); %set(gca,'YDir','normal');
scatter(T, spectrogram_ee',5,'y'); % Mark the peaks
ylabel('Frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('Spectrogram of whitened data','FontWeight','bold','Color','w');
axis tight;
hold off;

% Calculate spectrogram
[S,F,T,P]=spectrogram(whitenedSignal_mm{4},w,overlap,nfft,motu_Fs);
[S2]=spectrogram(whitenedSignal_mm{4},dw,overlap,nfft,motu_Fs);
%[S,F,T,P]=spectrogram(whitenedSignal_mm{1},w,overlap,nfft,motu_Fs);
%[S2]=spectrogram(whitenedSignal_mm{1},dw,overlap,nfft,motu_Fs);
IMAGE=100*((abs(S)+abs(S2))/2);
IMAGE=log(IMAGE+2);
IMAGE(IMAGE>high)=high;
IMAGE(IMAGE<low)=low;
IMAGE=IMAGE-low;
IMAGE=IMAGE/(high-low);
IMAGE=63*(IMAGE); % 63 is the cutoff for GIF
IMAGE_MODS = log(abs(IMAGE)+1e+2);
% Above from fb_pretty_sonogram code in FinchScope repo
IMAGE_MODS_SUM = sum(IMAGE_MODS);

spectrogram_ee_idxs = floor(aze{4}/(length(sum_whitenedSignal_allmics)/length(T)));
spectrogram_ee = NaN(1,length(T));
if spectrogram_ee_idxs(1) == 0; spectrogram_ee_idxs(1) = 1; end;
spectrogram_ee(spectrogram_ee_idxs) = 60000;

% Plot the spectrogram with the peaks highlighted
subplot(4,1,4); hold on; 
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gcf, 'Color', 'k'); 
colormap(hot)
imagesc(T,F,IMAGE_MODS); %set(gca,'YDir','normal');
scatter(T, spectrogram_ee',5,'y'); % Mark the peaks
ylabel('Frequency kHz','FontWeight','bold');
xlabel('time (s)','FontWeight','bold');
title('Spectrogram of whitened data','FontWeight','bold','Color','w');
axis tight;
hold off;