function [simag] = HilbertTransfer(sreal)
S = fft(sreal);
for i = 1:length(S)
    if i < length(S)/2
        S(i) = S(i) * -1i;
    elseif i > length(S)/2
        S(i) = S(i) * 1i;
    end
end
y = ifft(S,'symmetric'); % 由于实数序列的FFT满足共轭对称性，所以反FFT时可以指定'symmetric'参数，告诉Matlab我这个S是共轭对称的，从而使Matlab知道ifft的结果应该是实序列
simag = sreal + 1i*y;
end
