clc
clear all
close all
info = edfinfo('r01.edf');
data = edfread('r01.edf');
direct = edfread('r01.edf', 'SelectedSignals', 'Direct_1');
directarray = table2array(direct);
directecg = vertcat(directarray{:});
figure(1)
subplot(2,1,1)
plot(directecg(1:1920))% reference
title('Thoraicic ECG');
abdo = edfread('r01.edf', 'SelectedSignals', 'Abdomen_1');
abdoarray = table2array(abdo);%input to filter
abdoecg = vertcat(abdoarray{:});
subplot(2,1,2);
plot(abdoecg(1:1920));
mecg=abdoecg(1:1920);
title('Abdomen ECG')
filt = dsp.FIRFilter('Numerator',fir1(10,.25));
[a ,d]= swt(abdoecg(1:1920),5,'db4');
a=a(2,:);
d=d(2,:);
figure(2)
subplot(1,2,1);
plot(a);
title('approximate coeff');
subplot(1,2,2);
plot(d);
title('detailed coeff');
rls = dsp.RLSFilter('ForgettingFactor',0.99);
ta=transpose(a);
td=transpose(d);
[y1,e] = rls(directecg(1:1920),ta);
[y2,e2]=rls(directecg(1:1920),td);
figure(3)
subplot(2,1,1)
plot(1:1920, [mecg(1:1920),y1]);
title('System Identification of an FIR filter - a');
legend('Desired', 'Output');
xlabel('time index');
ylabel('signal value');
subplot(2,1,2)
plot(1:1920,[abdoecg(1:1920),y2])
title('System Identification of an FIR filter- d');
legend('Desired', 'Output');
xlabel('time index');
ylabel('signal value');
l=y1;
thr=3;
ythard = wthresh(l,'h',thr);
y2hard= wthresh(y2,'h',thr);
ytsoft = wthresh(l,'s',thr);
y2soft=wthresh(y2,'s',thr);
recon=iswt(transpose(ythard),transpose(y2hard),'db1');
plot(recon);
title('extracted fetal ecg');
P_original_signal = var(abdoecg(1:1920));
P_filtered_signal = var(recon(1:1920));
P_noise_before = var(abdoecg(1:1920)-transpose(recon(1:1920))); 
P_noise_after = var(transpose(recon(1:1920))-abdoecg(1:1920));
SNR_before_dB = 10 * log10(P_original_signal / P_noise_before);
SNR_after_dB = 10 * log10(P_filtered_signal / P_noise_after);
fprintf('SNR before filtering: %.2f dB\n', abs(SNR_before_dB));
fprintf('SNR after filtering: %.2f dB\n', abs(SNR_after_dB));
correlation_matrix = corrcoef(directecg(1:1920), recon);
correlation_coefficient = abs(correlation_matrix(1, 2));
disp(['Correlation Coefficient: ', num2str(correlation_coefficient)]);
ECG_inverted = -recon;
[~,locs_Rwave] = findpeaks(recon,'MinPeakHeight',0.5,'MinPeakDistance',200);
[~,locs_Swave] = findpeaks(ECG_inverted,'MinPeakHeight',0.5,'MinPeakDistance',200);
figure
hold on
plot(1:1920,recon)
plot(locs_Rwave,recon(locs_Rwave),'rv','MarkerFaceColor','r')
plot(locs_Swave,recon(locs_Swave),'rs','MarkerFaceColor','b')
grid on
legend('ECG Signal','R waves','S waves')
xlabel('Samples')
ylabel('Voltage(mV)')
title('R wave and S wave in Noisy ECG Signal')
rls = dsp.RLSFilter('ForgettingFactor',0.99);
ta1=transpose(a1);
td1=transpose(d1);
[y1,e] = rls(abdoecg(1:1920),ta1);
[y2,e2]=rls(directecg(1:1920),td1);
trecon=transpose(recon);
[y3,e] = rls(trecon,ta1);
[y4,e2]=rls(trecon,td1);
l=y3;
thr=2;
ythard1 = wthresh(l,'h',thr);
y4hard= wthresh(y4,'h',thr);
y3soft1 = wthresh(l,'s',thr);
y4soft=wthresh(y4,'s',thr);
recon1=iswt(transpose(ythard1),transpose(y4hard),'db1');
figure
plot(recon1);
title('extracted maternal ecg');
correlation_matrix1 = corrcoef(directecg(1:1920),recon1);
correlation_coefficient1 = correlation_matrix(1, 2);
disp(['Correlation Coefficient1: ', num2str(correlation_coefficient1)]);

%Extracting maternal ECG
filt = dsp.FIRFilter('Numerator',fir1(10,.25));
[a1,d1]= swt(abdoecg(1:1920),5,'db1');
a1=a1(2,:);
d1=d1(2,:);
rls = dsp.RLSFilter('ForgettingFactor',0.99);
ta1=transpose(a1);
td1=transpose(d1);
[y1,e] = rls(abdoecg(1:1920),ta1);
[y2,e2]=rls(directecg(1:1920),td1);
trecon=transpose(recon);
[y3,e] = rls(trecon,ta1);
[y4,e2]=rls(trecon,td1);
l=y3;
thr=2;
ythard1 = wthresh(l,'h',thr);
y4hard= wthresh(y4,'h',thr);
y3soft1 = wthresh(l,'s',thr);
y4soft=wthresh(y4,'s',thr);
recon1=iswt(transpose(ythard1),transpose(y4hard),'db1');
figure
plot(recon1);
title('extracted maternal ecg');
correlation_matrix1 = corrcoef(directecg(1:1920),recon1);
correlation_coefficient1 = correlation_matrix(1, 2);
disp(['Correlation Coefficient1: ', num2str(correlation_coefficient1)]);

