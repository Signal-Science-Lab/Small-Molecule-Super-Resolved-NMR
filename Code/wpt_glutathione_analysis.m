cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift\Data'

data = readtable('glutathione_conv.txt');
data = data(((data.ppm > 0.5) & (data.ppm <= 5.0)), :);
data.intensity = data.intensity / max(data.intensity);

% undecimated wavelet packet transform
n = 7;
wpdata = wpdec(data.intensity, n, 'db9');

rwpc = wprcoef(wpdata, 2^n-1);
figure;
plot(data.ppm, rwpc, 'b','linewidth',1);

cd 'C:\Users\as836\Documents\SignalScienceLab\SciAdv_NMR_Shift\Analysis'
writetable(table(data.ppm, rwpc), 'gluta_wpt_db9L7.txt');

for i = 4:7
    wpdata = wpdec(data.intensity, i, 'db9');    
    rwpc = wprcoef(wpdata, 2^i-1);
    fname = join(['gluta_wpt_db9L',int2str(i),'.txt']);
    writetable(table(data.ppm, rwpc), fname);
end