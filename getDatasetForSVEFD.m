% 参考文献：
%   [1] 曲志昱,毛校洁,侯长波.基于奇异值熵和分形维数的雷达信号识别[J].系统工程与电子技术,2018,40(02):303-307.
%% generate trainset for SVEFD
clc;clear all;
disp("generating trainset for SVEFD ...");
SignalType = 2;
SNRmin = -5;  
SNRmax = 18;
SNRstep = 1;
Samples = 50;
getDataset(1.1,[1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],SignalType,SNRmin,SNRmax,SNRstep,Samples,nan);
copyfile("训练集","训练集(CWD)");
movefile("训练集","训练集(FFT)");
FourierAnalysis("训练集(FFT)");
TimeFrequencyAnalysis("训练集(CWD)",2,nan,1);
getFractalDimension([1 1],["训练集(FFT)"]);    % 计算频域序列的分形维数
Preproccess([0 0 0 1],["训练集(CWD)"],4,0.15); % 将训练集(CWD)文件夹内的灰度图像的小于pixel_th的像素置为0
getSingularValueEntropy(["训练集(CWD)"]);
Preproccess([0 1 0 0],["训练集(CWD)","训练集(FFT)"],4,nan);
ConcatenateFeatureVector(["训练集(FFT)" "训练集(CWD)"],"训练集");
%% generate testset for SVEFD
