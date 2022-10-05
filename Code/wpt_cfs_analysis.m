%% Analysis of individual spectra %%

cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Data'

data1 = readtable('cellobiose_conv.txt');
data1 = data1(((data1.ppm > 1.0) & (data1.ppm <= 5.8)), :);
data1.intensity = data1.intensity / sum(data1.intensity);

data2 = readtable('fructose_conv.txt');
data2 = data2(((data2.ppm > 1.0) & (data2.ppm <= 5.8)), :);
data2.intensity = data2.intensity / sum(data2.intensity);

data3 = readtable('sucrose_conv.txt');
data3 = data3(((data3.ppm > 1.0) & (data3.ppm <= 5.8)), :);
data3.intensity = data3.intensity / sum(data3.intensity);

% undecimated wavelet packet transform
n1 = 7; n2 = 7; n3 = 8;
wpdata1 = wpdec(data1.intensity, n1, 'db9');
wpdata2 = wpdec(data2.intensity, n2, 'db9');
wpdata3 = wpdec(data3.intensity, n3, 'db9');

cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Analysis'

rwpc1 = wprcoef(wpdata1, 2^n1-1);
figure;
plot(data1.ppm, rwpc1, 'b','linewidth',1);
axis tight;
writetable(table(data1.ppm, rwpc1), 'cellobiose_wpt_db9L7.txt')

rwpc2 = wprcoef(wpdata2, 2^n2-1);
figure;
plot(data2.ppm, rwpc2, 'b','linewidth',1);
axis tight;
writetable(table(data2.ppm, rwpc2), 'fructose_wpt_db9L7.txt')

rwpc3 = wprcoef(wpdata3, 2^n3-1);
figure;
plot(data3.ppm, rwpc3, 'b','linewidth',1);
axis tight;
writetable(table(data3.ppm, rwpc3), 'sucrose_wpt_db9L8.txt')

%% Analysis of noise-free mixture data %%
cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Data'

mix1 = readtable('cfs_mix0_vf.txt');
mix1 = mix1(((mix1.ppm > 1.0) & (mix1.ppm <= 5.8)), :);

% undecimated wavelet packet transform
m1 = 8;
wpmix1 = wpdec(mix1.intensity, m1, 'db9');

cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Analysis'

rwpcm1 = wprcoef(wpmix1, 2^m1-1);
figure;
plot(mix1.ppm, rwpcm1, 'b','linewidth',1);
axis tight;
writetable(table(mix1.ppm, rwpcm1), 'cfs_wpt_db9L8.txt')

% Analysis of the data slice between 3.78 ppm and 3.94 ppm %
mix1A = mix1(((mix1.ppm > 3.78) & (mix1.ppm <= 3.94)), :);

m1A = 7;
wpmix1A = wpdec(mix1A.intensity, m1A, 'db9');

rwpcm1A = wprcoef(wpmix1A, 2^m1A-1);
figure;
plot(mix1A.ppm, rwpcm1A, 'b','linewidth',1);
axis tight;
writetable(table(mix1A.ppm, rwpcm1A), 'cfs_slice_wpt_db9L7.txt')

%% Analysis of noisy mixture data %%
cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Data'

mix2 = readtable('cfs_mix0_noisy_vf.txt');
mix2 = mix2(((mix2.ppm > 1.0) & (mix2.ppm <= 5.8)), :);

% undecimated wavelet packet transform
m2 = 8;
wpmix2 = wpdec(mix2.intensity, m2, 'db9');

cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Analysis'

rwpcm2 = wprcoef(wpmix2, 2^m2-1);
figure;
plot(mix2.ppm, rwpcm2, 'b','linewidth',1);
axis tight;
writetable(table(mix2.ppm, rwpcm2), 'cfs_noisy_wpt_db9L8.txt')

% Analysis of the data slice between 3.78 ppm and 3.94 ppm %
mix2A = mix2(((mix2.ppm > 3.78) & (mix2.ppm <= 3.94)), :);

m2A = 7;
wpmix2A = wpdec(mix2A.intensity, m2A, 'db9');

rwpcm2A = wprcoef(wpmix2A, 2^m2A-1);
figure;
plot(mix2A.ppm, rwpcm2A, 'b','linewidth',1);
axis tight;
writetable(table(mix2A.ppm, rwpcm2A), 'cfs_noisy_slice_wpt_db9L7.txt')