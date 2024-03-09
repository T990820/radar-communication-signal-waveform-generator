% Function: getDataset
% Author:   Nongkai Tian
% Date:     2024/1/17 12:13
% Instruments:
%   有如下的2种调用形式：
%       getDataset(Mode, Modulations, SignalType, SNRmin, SNRmax, SNRstep, Samples, snrth)
%       getDataset(Mode, Modulations, SignalType, SNR, Samples, snrth)
%   SignalType
%       SignalType = 1: real signal
%       SignalType = 2: imag signal
%   Mode
%       Mode = 1:   为图像识别网络的训练集生成时域信号，并将每个类的样本存到单独的文件夹里
%       Mode = 1.1: 为图像识别网络的训练集生成时域信号，并将所有样本存到一个文件夹里，样本文件命名为'BPSK_0005.mat'
%       Mode = 2:   为图像识别网络的测试集生成时域信号，将每dB不同调制方式的样本存到不同的文件夹里
%       Mode = 2.1：为图像识别网络的测试集生成时域信号，将每dB不同调制方式的样本存到一个文件夹里，样本文件命名为'BPSK_0005.mat'
%       Mode = 3:   为图像去噪网络的训练集生成时域信号
%       Mode = 4:   为图像去噪网络的测试集生成时域信号
%       Mode = 5:   为噪声估计网络的训练集生成时域信号
%   SNRmin：信噪比的最小值
%   SNRmax：信噪比的最大值
%   Modulations: 23 elements represents NS, EQFM, LFM, BPSK, QPSK, BFSK, SFM, COSTAS, TANFM, FRANK, P1, P2, P3, P4, T1, T2, T3, T4, LFM-SFM, EQFM-SFM, DDC-MASK, DDC-MPSK, DDC-MFSK respectively. 1 means generate, 0 means do not generate
% Reference:
%   [1] https://www.radartutorial.eu/08.transmitters/Frank%20Code.en.html
%   [2] J. E. Fielding, "Polytime coding as a means of pulse compression," in IEEE Transactions on Aerospace and Electronic Systems, vol. 35, no. 2, pp. 716-721, April 1999, doi: 10.1109/7.766951.
%   [3] https://www.radartutorial.eu/08.transmitters/Barker%20Code.en.html
function [] = getDataset(varargin)
disp("generating time domain signal ...");
if nargin == 8
    Mode = varargin{1};
    Modulations = varargin{2};
    SignalType  = varargin{3};
    SNRmin      = varargin{4};
    SNRmax      = varargin{5};
    SNRstep     = varargin{6};
    Samples     = varargin{7};
    snrth       = varargin{8};
    SNR         = SNRmin:SNRstep:SNRmax;
elseif nargin == 6
    Mode        = varargin{1};
    Modulations = varargin{2};
    SignalType  = varargin{3};
    SNR         = varargin{4};
    Samples     = varargin{5};
    snrth       = varargin{6};
end
warning ('off', 'all');
if size(Samples,1)==1 && size(Samples,2)==1 % 如果Samples是个数而不是矩阵，则将其拓展为向量
    Samples = ones(1,length(Modulations))*Samples;
end
AllFiles = dir(pwd);
disp("deleting existing folders ...");
for i = 3:length(AllFiles) % delete all generated folders
    if length(AllFiles(i).name) <= 2
        disp('DANGER!')
        return
    end
    if AllFiles(i).isdir == 1
        rmdir(AllFiles(i).name,'s')
    end
end
CreateFolders(Mode,Modulations,SNR);
param_dict = getParameterDict("信号参数.xlsx");      % 参数字典
fs         = getParameter(param_dict, '总体', 'fs'); % sampling frequency
n          = getParameter(param_dict, '总体', 'n');  % 默认采样点数(可覆盖)
TimeWidth  = n / fs;                                 % 时域宽度
v          = getParameter(param_dict, '总体', 'v');  % 信号电平
barkers = struct( ...
    'barker_5',  [1 1 1 0 1], ...              % 5位Barker码 [3]
    'barker_7',  [1 1 1 0 0 1 0], ...          % 7位Barker码 [3]
    'barker_11', [1 1 1 0 0 0 1 0 0 1 0], ...  % 11位Barker码[3]
    'barker_13', [1 1 1 1 1 0 0 1 1 0 1 0 1]); % 13位Barker码[3]
M_for_DDC_MASK = [4];
M_for_DDC_MPSK = [2 4];
M_for_DDC_MFSK = [2 4 8];
if Modulations(21) == 0
    M_for_DDC_MASK = [];
end
if Modulations(22) == 0
    M_for_DDC_MPSK = [];
end
if Modulations(23) == 0
    M_for_DDC_MFSK = [];
end
samples_num = getSampleNums(Modulations, Samples, SNR, M_for_DDC_MASK, M_for_DDC_MPSK, M_for_DDC_MFSK); % 将要生成的样本的数量
h          = waitbar(0,'Initializing','name','Time-Frequency Analysis');
sg = SignalGenerator(fs,v,n,TimeWidth);
fprintf('params for all modulated signals:\n\tsampling frequency(fs):\t\t%.2fHz\n\tsignal Vcc(v):\t\t\t\t%.2fV\n\tdefault sampling points(n):\t%d\n',fs,v,n);
fprintf('%d samples are expected to be generated\n',samples_num);
GlobalIndex    = 0; % 全局索引
%% NS
if Modulations(1) == 1
    FileIndex = 1; % 局部索引
    fc_min = getParameter(param_dict, 'NS', 'fc_min');
    fc_max = getParameter(param_dict, 'NS', 'fc_max');
    n_min = getParameter(param_dict, 'NS', 'n_min');
    n_max = getParameter(param_dict, 'NS', 'n_max');
    fprintf('params for NS:\n\tcarrier frequency(fc):\t[%.2fHz, %.2fHz]\n\tsamples number(N₀):\t\t[%d, %d]\n',fc_min*fs,fc_max*fs,n_min,n_max);
    for snr = SNR
        for i = 1:Samples(1)
            waitbar(GlobalIndex/samples_num,h,['NS-' num2str(snr) 'dB-' num2str(i)]);
            fc = (fc_min + (fc_max - fc_min) * rand()) * fs; % 载波频率
            n = floor(n_min + (n_max - n_min) * rand());     % 采样点数
            time_width = n / fs;                             % 时域宽度
            [s,noise] = sg.generateNS(time_width, n, fc,snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'NS',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% EQFM
if Modulations(2) == 1
    FileIndex = 1; % 局部索引
    f_valley_min = getParameter(param_dict, 'EQFM', 'f_valley_min');
    f_valley_max = getParameter(param_dict, 'EQFM', 'f_valley_max');
    band_width_min = getParameter(param_dict, 'EQFM', 'band_width_min');
    band_width_max = getParameter(param_dict, 'EQFM', 'band_width_max');
    n_min = getParameter(param_dict, 'EQFM', 'n_min');
    n_max = getParameter(param_dict, 'EQFM', 'n_max');
    fprintf('params for EQFM:\n\tvalley frequency(fmin):\t[%.2fHz, %.2fHz]\n\tsamples number(N₀):\t\t[%d, %d]\n\tband width(B):\t\t\t[%.2fHz, %.2fHz]\n',f_valley_min*fs,f_valley_max*fs,n_min,n_max,band_width_min*fs,band_width_max*fs);
    for snr = SNR
        for i = 1:Samples(2)
            waitbar(GlobalIndex/samples_num,h,['EQFM-' num2str(snr) 'dB-' num2str(i)]);
            f_valley = (f_valley_min + (f_valley_max - f_valley_min) * rand())*fs;         % 谷值频率
            n = floor(n_min + (n_max - n_min) * rand());                                   % 采样点数
            band_width = (band_width_min + (band_width_max - band_width_min) * rand())*fs; % 频域宽度
            while(f_valley+band_width>0.5*fs)
                f_valley = (f_valley_min + (f_valley_max - f_valley_min) * rand())*fs;
                band_width = (band_width_min + (band_width_max - band_width_min) * rand())*fs;
            end
            time_width = n / fs;                                                           % 时域宽度
            [s,noise] = sg.generateEQFM(f_valley,n,band_width,time_width,snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'EQFM',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% LFM
if Modulations(3) == 1
    FileIndex = 1;
    f0_min = getParameter(param_dict, 'LFM', 'f0_min');
    f0_max = getParameter(param_dict, 'LFM', 'f0_max');
    band_width_min = getParameter(param_dict, 'LFM', 'band_width_min');
    band_width_max = getParameter(param_dict, 'LFM', 'band_width_max');
    n_min = getParameter(param_dict, 'LFM', 'n_min');
    n_max = getParameter(param_dict, 'LFM', 'n_max');
    fprintf('params for LFM:\n\tinitial frequency(f0):\t[%.2fHz, %.2fHz]\n\tsamples number(N₀):\t\t[%d, %d]\n\tband width(B):\t\t\t[%.2fHz, %.2fHz]\n',f0_min*fs,f0_max*fs,n_min,n_max,band_width_min*fs,band_width_max*fs);
    for snr = SNR
        for i = 1:Samples(3)
            waitbar(GlobalIndex/samples_num,h,['LFM-' num2str(snr) 'dB-' num2str(i)]);
            f0 = (f0_min+(f0_max-f0_min)*rand())*fs;                                        % 起始频率
            band_width = (band_width_min + (band_width_max - band_width_min) * rand()) *fs; % 频域宽度
            while(f0+band_width>0.5*fs)
                f0 = (f0_min+(f0_max-f0_min)*rand())*fs;
                band_width = (band_width_min + (band_width_max - band_width_min) * rand()) *fs;
            end
            n = floor(n_min + (n_max - n_min) * rand());                                    % 采样点数
            [s,noise] = sg.generateLFM(snr,f0,band_width,n,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'LFM',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% BPSK
if Modulations(4) == 1
    FileIndex = 1;
    fc_min = getParameter(param_dict, 'BPSK', 'fc_min');
    fc_max = getParameter(param_dict, 'BPSK', 'fc_max');
    n_min = getParameter(param_dict, 'BPSK', 'n_min');
    n_max = getParameter(param_dict, 'BPSK', 'n_max');
    for snr = SNR
        for i = 1:Samples(4)
            waitbar(GlobalIndex/samples_num,h,['BPSK-' num2str(snr) 'dB-' num2str(i)]);
            % 随机生成5/7/11/13位Barker码
            barker_length = [5 7 11 13];
            random_length = barker_length(randperm(length(barker_length)));
            barker = barkers.(['barker_' num2str(random_length(1))]);
            fc = (fc_min+rand()*(fc_max-fc_min))*fs;     % 载波频率
            n = floor(n_min + (n_max - n_min) * rand()); % 采样点数
            time_width = n / fs;                         % 时域宽度
            t = linspace(0,time_width/length(barker),n/length(barker));
            [s,noise] = sg.generateBPSK(barker,fc,t,snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'BPSK',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% QPSK
if Modulations(5) == 1
    FileIndex = 1;
    code = [0 0 0 0 0 0 0 0 ...
        0 0 0 1 1 0 1 1 ...
        0 0 1 0 0 0 1 0 ...
        0 0 1 1 1 0 0 1]; % 16-bit Frank Code[1]
    for snr = SNR
        for i = 1:Samples(5)
            waitbar(GlobalIndex/samples_num,h,['QPSK-' num2str(snr) 'dB-' num2str(i)]);
            [s,noise] = sg.generateQPSK(snr,code,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'QPSK',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% BFSK
if Modulations(6) == 1
    FileIndex = 1;
    fc_min = getParameter(param_dict, 'BFSK', 'fc_min');
    fc_max = getParameter(param_dict, 'BFSK', 'fc_max');
    use_barker = getParameter(param_dict, 'BFSK', 'use_barker');
    use_random_code = getParameter(param_dict, 'BFSK', 'use_random_code');
    assert(GlobalIndex~=samples_num, "码元选择不正确！")
    for snr = SNR
        for i=1:Samples(6)
            waitbar(GlobalIndex/samples_num,h,['FSK-' num2str(snr) 'dB-' num2str(i)]);
            if use_barker == 1
                % 随机生成5/7/11/13位Barker码
                barker_length = [5 7 11 13];
                random_length = barker_length(randperm(length(barker_length)));
                code = barkers.(['barker_' num2str(random_length(1))]);
            elseif use_random_code ~= 0
                code = rand(1,use_random_code) > 0.5;
                while(sum(code)==0 || sum(code)==length(code))
                    code = rand(1,use_random_code) > 0.5;
                end
            end
            fc = (fc_min+rand()*(fc_max-fc_min))*fs; % 载波频率
            [s,noise] = sg.generateBFSK(fc,code,snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'FSK',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% SFM
if Modulations(7) == 1
    FileIndex = 1;
    f_peak_min = getParameter(param_dict, 'SFM', 'f_peak_min');
    f_peak_max = getParameter(param_dict, 'SFM', 'f_peak_max');
    f_valley_min = getParameter(param_dict, 'SFM', 'f_valley_min');
    f_valley_max = getParameter(param_dict, 'SFM', 'f_valley_max');
    for snr = SNR
        for i = 1:Samples(7)
            waitbar(GlobalIndex/samples_num,h,['SFM-' num2str(snr) 'dB-' num2str(i)]);
            fmin = (f_valley_min+(f_valley_max-f_valley_min)*rand())*fs;
            fmax = (f_peak_min+(f_peak_max-f_peak_min)*rand())*fs;
            ff   = (1+rand(1)*2)/TimeWidth; % frequency's frequency, frequency changes 1 to 3 times per period
            [s,noise] = sg.generateSFM(ff,fmin,fmax,snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'SFM',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% COSTAS
if Modulations(8) == 1
    FileIndex = 1;
    COSTAS_sequence_length_min = getParameter(param_dict, 'COSTAS', 'COSTAS_sequence_length_min');
    COSTAS_sequence_length_max = getParameter(param_dict, 'COSTAS', 'COSTAS_sequence_length_max');
    f_min = getParameter(param_dict, 'COSTAS', 'f_min')*fs;
    f_max = getParameter(param_dict, 'COSTAS', 'f_max')*fs;
    fprintf('params for COSTAS:\n\tcarrier frequency(fc):\t[%.2fHz, %.2fHz]\n\tcostas array length(N₀):\t[%d, %d]\n\t\n',f_min,f_max,COSTAS_sequence_length_min,COSTAS_sequence_length_max);
    for snr = SNR
        for i = 1:Samples(8)
            waitbar(GlobalIndex/samples_num,h,['COSTAS-' num2str(snr) 'dB-' num2str(i)]);
            COSTAS_sequence_length = floor(COSTAS_sequence_length_min + rand() * (COSTAS_sequence_length_max - COSTAS_sequence_length_min));
            [s,noise] = sg.generateCOSTAS(COSTAS_sequence_length,f_min,f_max,snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'COSTAS',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% TANFM
if Modulations(9) == 1
    FileIndex = 1;
    for snr = SNR
        for i = 1:Samples(9)
            waitbar(GlobalIndex/samples_num,h,['TANFM-' num2str(snr) 'dB-' num2str(i)]);
            [s,noise] = sg.generateTANFM(snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'TANFM',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% Frank[1]
if Modulations(10) == 1
    FileIndex = 1;
    fc_min = getParameter(param_dict, 'FRANK', 'fc_min');
    fc_max = getParameter(param_dict, 'FRANK', 'fc_max');
    M_min = getParameter(param_dict, 'FRANK', 'M_min');
    M_max = getParameter(param_dict, 'FRANK', 'M_max');
    n_min = getParameter(param_dict, 'FRANK', 'n_min');
    n_max = getParameter(param_dict, 'FRANK', 'n_max');
    for snr = SNR
        for i = 1:Samples(10)
            waitbar(GlobalIndex/samples_num,h,['Frank-' num2str(snr) 'dB-' num2str(i)]);
            fc = (fc_min+rand()*(fc_max-fc_min))*fs;     % 载波频率
            n = floor(n_min + (n_max - n_min) * rand()); % 采样点数
            time_width = n / fs;                         % 时域宽度
            M = ceil(M_min + (M_max - M_min) * rand());  % 相位点数
            [s,noise] = sg.generateFrank(M,fc,n,time_width,snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'FRANK',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% P1
if Modulations(11) == 1
    FileIndex = 1;
    N_min = getParameter(param_dict, 'P1', 'N_min');
    N_max = getParameter(param_dict, 'P1', 'N_max');
    for snr = SNR
        for i = 1:Samples(11)
            waitbar(GlobalIndex/samples_num,h,['P1-' num2str(snr) 'dB-' num2str(i)]);
            N = floor(N_min + (N_max - N_min) * rand()); % 采样点数
            [s,noise] = sg.generateP1(snr,N,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'P1',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% P2
if Modulations(12) == 1
    FileIndex = 1;
    for snr = SNR
        for i = 1:Samples(12)
            waitbar(GlobalIndex/samples_num,h,['P2-' num2str(snr) 'dB-' num2str(i)]);
            [s,noise] = sg.generateP2(snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'P2',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% P3
if Modulations(13) == 1
    FileIndex = 1;
    for snr = SNR
        for i = 1:Samples(13)
            waitbar(GlobalIndex/samples_num,h,['P3-' num2str(snr) 'dB-' num2str(i)]);
            [s,noise] = sg.generateP3(snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'P3',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% P4
if Modulations(14) == 1
    FileIndex = 1;
    for snr = SNR
        for i = 1:Samples(14)
            waitbar(GlobalIndex/samples_num,h,['P4-' num2str(snr) 'dB-' num2str(i)]);
            [s,noise] = sg.generateP4(snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'P4',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% T1[2]
if Modulations(15) == 1
    FileIndex = 1;
    fc_min = getParameter(param_dict, 'T1', 'fc_min');
    fc_max = getParameter(param_dict, 'T1', 'fc_max');
    k_min = getParameter(param_dict, 'T1', 'k_min');
    k_max = getParameter(param_dict, 'T1', 'k_max');
    n_min = getParameter(param_dict, 'T1', 'n_min');
    n_max = getParameter(param_dict, 'T1', 'n_max');
    for snr = SNR
        for i = 1:Samples(15)
            waitbar(GlobalIndex/samples_num,h,['T1-' num2str(snr) 'dB-' num2str(i)]);
            fc = (fc_min+rand()*(fc_max-fc_min))*fs;
            k = floor(k_min+rand()*(k_max-k_min));
            n = floor(n_min+rand()*(n_max-n_min));
            [s,noise] = sg.generateT1(snr,fc,k,n,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'T1',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% T2[2]
if Modulations(16) == 1
    FileIndex = 1;
    fc_min = getParameter(param_dict, 'T2', 'fc_min');
    fc_max = getParameter(param_dict, 'T2', 'fc_max');
    k_min = getParameter(param_dict, 'T2', 'k_min');
    k_max = getParameter(param_dict, 'T2', 'k_max');
    n_min = getParameter(param_dict, 'T2', 'n_min');
    n_max = getParameter(param_dict, 'T2', 'n_max');
    for snr = SNR
        for i = 1:Samples(16)
            waitbar(GlobalIndex/samples_num,h,['T2-' num2str(snr) 'dB-' num2str(i)]);
            fc = (fc_min+rand()*(fc_max-fc_min))*fs;
            k = floor(k_min+rand()*(k_max-k_min));
            n = floor(n_min+rand()*(n_max-n_min));
            [s,noise] = sg.generateT2(snr,fc,k,n,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'T2',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% T3[2]
if Modulations(17) == 1
    FileIndex = 1;
    fc_min = getParameter(param_dict, 'T3', 'fc_min');
    fc_max = getParameter(param_dict, 'T3', 'fc_max');
    n_min = getParameter(param_dict, 'T3', 'n_min');
    n_max = getParameter(param_dict, 'T3', 'n_max');
    for snr = SNR
        for i = 1:Samples(17)
            waitbar(GlobalIndex/samples_num,h,['T3-' num2str(snr) 'dB-' num2str(i)]);
            fc = (fc_min+rand()*(fc_max-fc_min))*fs;
            n = floor(n_min+rand()*(n_max-n_min));
            [s,noise] = sg.generateT3(snr,fc,n,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'T3',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% T4[2]
if Modulations(18) == 1
    FileIndex = 1;
    fc_min = getParameter(param_dict, 'T4', 'fc_min');
    fc_max = getParameter(param_dict, 'T4', 'fc_max');
    n_min = getParameter(param_dict, 'T4', 'n_min');
    n_max = getParameter(param_dict, 'T4', 'n_max');
    for snr = SNR
        for i = 1:Samples(18)
            waitbar(GlobalIndex/samples_num,h,['T4-' num2str(snr) 'dB-' num2str(i)]);
            fc = (fc_min+rand()*(fc_max-fc_min))*fs;
            n = floor(n_min+rand()*(n_max-n_min));
            [s,noise] = sg.generateT4(snr,fc,n,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'T4',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% LFM-SFM
if Modulations(19) == 1
    FileIndex = 1;
    for snr = SNR
        for i = 1:Samples(19)
            waitbar(GlobalIndex/samples_num,h,['LFM-SFM-' num2str(snr) 'dB-' num2str(i)]);
            [s,noise] = sg.generateLFM_SFM(snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'LFM-SFM',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% EQFM-SFM
if Modulations(20) == 1
    FileIndex = 1;
    for snr = SNR
        for i = 1:Samples(20)
            waitbar(GlobalIndex/samples_num,h,['EQFM-SFM-' num2str(snr) 'dB-' num2str(i)]);
            [s,noise] = sg.generateEQFM_SFM(snr,SignalType);
            SaveTimeDomainSignal(Mode,s,noise,snr,snrth,'EQFM-SFM',FileIndex,GlobalIndex,samples_num);
            FileIndex = FileIndex + 1;
            GlobalIndex = GlobalIndex + 1;
        end
    end
end
%% DDC_MASK
if Modulations(21) == 1
    for M = M_for_DDC_MASK
        FileIndex = 1;
        for snr = SNR
            for i = 1:Samples(21)
                waitbar(GlobalIndex/samples_num,h,['DDC-' num2str(M) 'ASK-' num2str(snr) 'dB-' num2str(i)]);
                symbol_rate = 12500;
                [s,noise] = sg.generateDDC_MASK(snr,1,M,symbol_rate,2);
                SaveTimeDomainSignal(Mode,s,noise,snr,snrth,['DDC-' num2str(M) 'ASK'],FileIndex,GlobalIndex,samples_num);
                FileIndex = FileIndex + 1;
                GlobalIndex = GlobalIndex + 1;
            end
        end
    end
end
%% DDC_MPSK
if Modulations(22) == 1
    for M = M_for_DDC_MPSK
        FileIndex = 1;
        for snr = SNR
            for i = 1:Samples(22)
                waitbar(GlobalIndex/samples_num,h,['DDC-' num2str(M) 'PSK-' num2str(snr) 'dB-' num2str(i)]);
                symbol_rate = 12500;
                [s,noise] = sg.generateDDC_MPSK(snr,1,M,symbol_rate,2);
                SaveTimeDomainSignal(Mode,s,noise,snr,snrth,['DDC-' num2str(M) 'PSK'],FileIndex,GlobalIndex,samples_num);
                FileIndex = FileIndex + 1;
                GlobalIndex = GlobalIndex + 1;
            end
        end
    end
end
%% DDC_MFSK
if Modulations(23) == 1
    for M = M_for_DDC_MFSK
        FileIndex = 1;
        for snr = SNR
            for i = 1:Samples(21)
                waitbar(GlobalIndex/samples_num,h,['DDC-' num2str(M) 'FSK-' num2str(snr) 'dB-' num2str(i)]);
                symbol_rate = 12500;
                [s,noise] = sg.generateDDC_MFSK(snr,1,M,100000,symbol_rate,2);
                SaveTimeDomainSignal(Mode,s,noise,snr,snrth,['DDC-' num2str(M) 'FSK'],FileIndex,GlobalIndex,samples_num);
                FileIndex = FileIndex + 1;
                GlobalIndex = GlobalIndex + 1;
            end
        end
    end
end
%% MQAM
if Modulations(24) == 1
    for M = M_for_DDC_MFSK
        FileIndex = 1;
        for snr = SNR
            for i = 1:Samples(21)
                waitbar(GlobalIndex/samples_num,h,[num2str(M) 'QAM-' num2str(snr) 'dB-' num2str(i)]);
                [s,noise] = sg.generateMQAM();
                SaveTimeDomainSignal(Mode,s,noise,snr,snrth,[num2str(M) 'QAM'],FileIndex,GlobalIndex,samples_num);
                FileIndex = FileIndex + 1;
                GlobalIndex = GlobalIndex + 1;
            end
        end
    end
end
delete(h);close all;
assert(GlobalIndex==samples_num, "样本个数计算不正确！")
end
