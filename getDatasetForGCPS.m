% 参考文献：
%   [1] Z. Yu, J. Tang and Z. Wang, "GCPS: A CNN Performance Evaluation Criterion for Radar Signal Intrapulse Modulation Recognition," in IEEE Communications Letters, vol. 25, no. 7, pp. 2290-2294, July 2021, doi: 10.1109/LCOMM.2021.3070151.
%% generate trainset for GCPS
clc;clear all;
disp("generating trainset for GCPS ...");
SignalType = 2;
SNRmin = -4;  
SNRmax = 14;
SNRstep = 2;
Samples = 250;
getDataset(1.1,[1 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],SignalType,SNRmin,SNRmax,SNRstep,Samples,nan);
TimeFrequencyAnalysis("训练集",1,[224 224],1);
Preproccess([0 1 1],["训练集"],5);
%% generate testset for GCPS
clc;clear all;
disp("generating testset for GCPS ...");
SignalType = 2;
SNRmin = -10;  
SNRmax = 10;
SNRstep = 1;
Samples = 200;
getDataset(2.1,[1 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],SignalType,SNRmin,SNRmax,SNRstep,Samples,nan);
TimeFrequencyAnalysis("测试集",1,[224 224],1);
Preproccess([0 0 1],["测试集"],5);
%% 为观察无噪声情况下的Grad-CAM制作样本
clc;clear all;
SignalType = 2;
getDataset(2.1, [1 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0], 2, inf, 1, nan)
TimeFrequencyAnalysis("测试集",1,[224 224],1);
Preproccess([0 0 1],["测试集"],5);