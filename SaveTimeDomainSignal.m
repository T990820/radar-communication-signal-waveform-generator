function [] = SaveTimeDomainSignal(mode,signal,noise,snr,snrth,modulation,file_index,global_index,file_nums)
% Instruments:
%   mode
%       mode = 1:   为图像识别网络的训练集生成时域信号，并将每个类的样本存到单独的文件夹里
%       mode = 1.1: 为图像识别网络的训练集生成时域信号，并将所有样本存到一个文件夹里，样本文件命名为'BPSK_0005.mat'
%       mode = 2:   为图像识别网络的测试集生成时域信号
%       mode = 3:   为图像去噪网络的训练集生成时域信号
%       mode = 4:   为图像去噪网络的测试集生成时域信号
%       mode = 5:   为噪声估计网络的训练集生成时域信号
%       mode = 7:   为复合调制信号识别任务的训练集生成时域信号，并将每个类的样本存到单独的文件夹中
%       Mode = 9:   为多分量调制信号去噪任务的训练集生成时域信号
%   modulation
%       字符串类型变量，表示signal的调制方式
%   file_nums：文件的总数，用来给文件名前补0
s = signal;
file_path_length = length(num2str(file_nums)); % 文件名的长度
if mode == 1 || mode == 7
    save(['训练集\\' modulation '\\' num2str(file_index) '.mat'],"s")
elseif mode == 1.1
    file_name = sprintf(['%0' num2str(file_path_length) 's'],num2str(global_index));
    save(['训练集\' modulation '_' file_name '.mat'],"s");
elseif mode == 2 || mode == 8
    save(['测试集\\' num2str(snr) 'dB\\' modulation '\\' num2str(file_index) '.mat'],"s")
elseif mode == 2.1 || mode == 6
    file_name = sprintf(['%0' num2str(file_path_length) 's'],num2str(global_index));
    save(['测试集\\' num2str(snr) 'dB\\' modulation '_' file_name '.mat'],"s");
elseif mode == 3
    file_name = sprintf(['%0' num2str(file_path_length) 's'],num2str(global_index));
    s = signal - noise;
    save(['Pure\\' modulation '_' file_name '.mat'], "s")
    s = signal;
    save(['Dirty\\' modulation '_' file_name '.mat'],"s")
    s = noise;
    save(['Noise\\' modulation '_' file_name '.mat'],"s")
elseif mode == 4
    file_name = sprintf(['%0' num2str(file_path_length) 's'],num2str(global_index));
    s = signal - noise;
    save(['Pure\\'  num2str(snr) 'dB\\' modulation '_' file_name '.mat'],"s");
    s = signal;
    save(['Dirty\\' num2str(snr) 'dB\\' modulation '_' file_name '.mat'],"s");
    s = noise;
    save(['Noise\\' num2str(snr) 'dB\\' modulation '_' file_name '.mat'],"s");
elseif mode == 5
    if snr > snrth
        save(['训练集\\低噪声水平\\' num2str(global_index) '.mat'],"s")
    else
        save(['训练集\\高噪声水平\\' num2str(global_index) '.mat'],"s")
    end
elseif mode == 9
    file_name = sprintf(['%0' num2str(file_path_length) 's'],num2str(global_index));
    s = signal - noise;
    save(['Pure\\' modulation '\\' file_name '_' modulation '.mat'], "s")
    s = signal;
    save(['Dirty\\' modulation '\\' file_name '_' modulation '.mat'],"s")
    s = noise;
    save(['Noise\\' modulation '\\' file_name '_' modulation '.mat'],"s")
end
