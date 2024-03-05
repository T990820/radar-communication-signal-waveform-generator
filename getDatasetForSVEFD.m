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
Preproccess([0 1 0 0],["训练集(CWD)","训练集(FFT)"],4,nan); % 按照相同的映射打乱训练集(CWD)和训练集(FFT)两文件夹中的文件
ConcatenateFeatureVector(["训练集(FFT)" "训练集(CWD)"],"训练集");
Preproccess([1 0 0 0],["训练集"],4,nan);
%% generate testset for SVEFD
clc;clear all;
disp("generating trainset for SVEFD ...");
SignalType = 2;
SNRmin = -15;  
SNRmax = 10;
SNRstep = 1;
Samples = 100;
file_path_length = 
getDataset(2.1,[1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],SignalType,SNRmin,SNRmax,SNRstep,Samples,nan);
copyfile("测试集","测试集(CWD)");
movefile("测试集","测试集(FFT)");
FourierAnalysis("测试集(FFT)");
TimeFrequencyAnalysis("测试集(CWD)",2,nan,1);
getFractalDimension([1 1],["测试集(FFT)"]);    % 计算频域序列的分形维数
Preproccess([0 0 0 1],["测试集(CWD)"],nan,0.15); % 将测试集(CWD)文件夹内的灰度图像的小于pixel_th的像素置为0
getSingularValueEntropy(["测试集(CWD)"]);
for snr = SNRmin:SNRstep:SNRmax
    ConcatenateFeatureVector([strcat("测试集(FFT)\Image",num2str(snr),"dB\") strcat("测试集(CWD)\Image",num2str(snr),"dB\")],strcat("测试集\",num2str(snr),"dB\"));
end
Preproccess([1 0 0 0],["测试集"],nan,nan);
%% generate validation set for SVEFD
clc;clear all;
disp("generating trainset for SVEFD ...");
SignalType = 2;
SNRmin = -5;  
SNRmax = 18;
SNRstep = 1;
Samples = 50;
getDataset(1.1,[1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0],SignalType,SNRmin,SNRmax,SNRstep,Samples,nan);
copyfile("训练集","验证集(CWD)");
movefile("训练集","验证集(FFT)");
FourierAnalysis("验证集(FFT)");
TimeFrequencyAnalysis("验证集(CWD)",2,nan,1);
getFractalDimension([1 1],["验证集(FFT)"]);    % 计算频域序列的分形维数
Preproccess([0 0 0 1],["验证集(CWD)"],4,0.15); % 将验证集(CWD)文件夹内的灰度图像的小于pixel_th的像素置为0
getSingularValueEntropy(["验证集(CWD)"]);
Preproccess([0 1 0 0],["验证集(CWD)","验证集(FFT)"],4,nan); % 按照相同的映射打乱验证集(CWD)和验证集(FFT)两文件夹中的文件
ConcatenateFeatureVector(["验证集(FFT)" "验证集(CWD)"],"验证集");
Preproccess([1 0 0 0],["验证集"],4,nan);