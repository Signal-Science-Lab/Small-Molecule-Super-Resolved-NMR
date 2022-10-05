%% Directories for reading and writing data %%
rdir = 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Data';
wdir = 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Analysis';

%% Analysis of individual spectra %%
files = ["glutamine","glycine","isoleucine","leucine","threonine","valine"];

for f = files
    cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Data'
    data = readtable(strjoin([f,'_conv.txt'],''));
    data.intensity = rescale(data.intensity);
    % undecimated wavelet packet transform
    n = 7;
    wpdata = wpdec(data.intensity, n, 'db9');
    cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Analysis'
    rwpc = wprcoef(wpdata, 2^n-1);
    figure;
    plot(data.ppm, rwpc, 'b','linewidth',1);
    writetable(table(data.ppm, rwpc), strjoin([f,"_wpt_db9L",int2str(n),".txt"],''));
end

%% Analysis of the mixture-1 %%
cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Data'
data = readtable('metabol_mix_noisy_vf.txt');
data.intensity = rescale(data.intensity);
% undecimated wavelet packet transform
n = 7;
wpdata = wpdec(data.intensity, n, 'db9');
cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Analysis'
rwpc = wprcoef(wpdata, 2^n-1);
figure;
plot(data.ppm, rwpc, 'b','linewidth',1);
writetable(table(data.ppm, rwpc), 'metabol_mix_wpt_db9L7.txt');

%% Analysis of the mixture-2 %%
cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Data'
data = readtable('metabol_mix2_noisy_vf.txt');
data.intensity = rescale(data.intensity);
% undecimated wavelet packet transform
n = 7;
wpdata = wpdec(data.intensity, n, 'db9');
cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift_vf\Analysis'
rwpc = wprcoef(wpdata, 2^n-1);
figure;
plot(data.ppm, rwpc, 'b','linewidth',1);
writetable(table(data.ppm, rwpc), 'metabol_mix2_wpt_db9L7.txt');