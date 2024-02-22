% Function: getDataset
% Author:   Nongkai Tian
% Date:     2023/12/28 16:00
% Instruments:
%   为生成的数据集创建文件夹，有如下两种调用模式
%       CreateFolders(Mode,Modulations,SNRmin,SNRmax)
%       CreateFolders(Mode,Modulations,SNR)
%   Mode
%       Mode = 1:   为图像识别网络的训练集生成时域信号，并将每个类的样本存到单独的文件夹里
%       Mode = 1.1: 为图像识别网络的训练集生成时域信号，并将所有样本存到一个文件夹里，样本文件命名为'BPSK_0005.mat'
%       Mode = 2:   为图像识别网络的测试集生成时域信号
%       Mode = 3:   为图像去噪网络的训练集生成时域信号
%       Mode = 4:   为图像去噪网络的测试集生成时域信号
%       Mode = 5:   为噪声估计网络的训练集生成时域信号
function [] = CreateFolders(varargin)
Mode = varargin{1};
Modulations = varargin{2};
if nargin == 4
    SNRmin = varargin{3}; SNRmax = varargin{4};
    SNR = SNRmin:SNRmax;
elseif nargin == 3
    SNR = varargin{3};
end
ModulationsStruct = struct( ...
    'a','NS',...
    'b','EQFM', ...
    'c','LFM', ...
    'd','BPSK', ...
    'e','QPSK', ...
    'f','FSK', ...
    'g','SFM', ...
    'h','COSTAS', ...
    'i','TANFM', ...
    'j','FRANK', ...
    'k','P1', ...
    'l','P2', ...
    'm','P3', ...
    'n','P4', ...
    'o','T1', ...
    'p','T2', ...
    'q','T3', ...
    'r','T4', ...
    's','LFM-SFM', ...
    't','EQFM-SFM');
if Mode == 1
    mkdir('训练集');
    for i = 1:length(Modulations)
        if Modulations(i) == 1
            mkdir(['训练集\' ModulationsStruct.(char('a'+i-1))]);
        end
    end
elseif Mode == 1.1
    mkdir('训练集\');
elseif Mode == 2
    mkdir('测试集');
    for i = SNR
        for j = 1:length(Modulations)
            if Modulations(j) == 1
                mkdir(['测试集\\Image' num2str(i) 'dB\\' ModulationsStruct.(char('a'+j-1))]);
            end
        end
    end
elseif Mode == 2.1
    mkdir('测试集');
    for i = SNR
        mkdir(['测试集\\Image' num2str(i) 'dB\\']);
    end
elseif Mode == 3
    mkdir("Pure\\");mkdir("Dirty\\");mkdir("Noise\\");
elseif Mode == 4
    for snr = SNR
        mkdir(['Pure\\'  num2str(snr) 'dB']);
        mkdir(['Dirty\\' num2str(snr) 'dB']);
        mkdir(['Noise\\' num2str(snr) 'dB']);
    end
elseif Mode == 5
    mkdir('训练集');mkdir('训练集\\低噪声水平');mkdir('训练集\\高噪声水平');
end