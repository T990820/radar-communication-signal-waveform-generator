function getMultiComponentSignal(Mode,Modulations,SNR)
%   Mode = 7: 为复合调制信号识别任务的训练集生成时域信号，并将每个类的样本存到单独的文件夹中
all_modulations = {"NS", "EQFM", "LFM", "VTFM", "BPSK", "QPSK", "BFSK", "SFM", "COSTAS", "TANFM", ...
    "FRANK", "P1", "P2", "P3", "P4", "T1", "T2", "T3", "T4", "LFM-SFM", ...
    "EQFM-SFM", "DDC-MASK", "DDC-MPSK", "DDC-MFSK", "MQAM"};
switch Mode
    case 7
        for i = 1 : length(Modulations)
            for j = i+1 : length(Modulations)
                if Modulations(i) == 1 && Modulations(j) == 1
                    addSingleComponentSignals(['训练集\' convertStringsToChars(all_modulations{i})], ...
                        ['训练集\' convertStringsToChars(all_modulations{j})], ...
                        ['训练集\' convertStringsToChars(all_modulations{i}) '+' convertStringsToChars(all_modulations{j})]);
                end
            end
        end
        deleteSingleComponentSignals("训练集");
    case 8
        for snr = SNR
            for i = 1 : length(Modulations)
                for j = i+1 : length(Modulations)
                    if Modulations(i) == 1 && Modulations(j) == 1
                        addSingleComponentSignals(['测试集\' num2str(snr) 'dB\' convertStringsToChars(all_modulations{i})], ...
                            ['测试集\' num2str(snr) 'dB\' convertStringsToChars(all_modulations{j})], ...
                            ['测试集\' num2str(snr) 'dB\' convertStringsToChars(all_modulations{i}) '+' convertStringsToChars(all_modulations{j})]);
                    end
                end
            end
            deleteSingleComponentSignals(['测试集\' num2str(snr) 'dB\']);
        end
    case 9
        for i = 1 : length(Modulations)
            for j = i+1 : length(Modulations)
                if Modulations(i) == 1 && Modulations(j) == 1
                    addSingleComponentSignals(['Dirty\' convertStringsToChars(all_modulations{i})], ...
                        ['Dirty\' convertStringsToChars(all_modulations{j})], ...
                        ['Dirty\' convertStringsToChars(all_modulations{i}) '+' convertStringsToChars(all_modulations{j})]);
                    addSingleComponentSignals(['Pure\' convertStringsToChars(all_modulations{i})], ...
                        ['Pure\' convertStringsToChars(all_modulations{j})], ...
                        ['Pure\' convertStringsToChars(all_modulations{i}) '+' convertStringsToChars(all_modulations{j})]);
                    addSingleComponentSignals(['Noise\' convertStringsToChars(all_modulations{i})], ...
                        ['Noise\' convertStringsToChars(all_modulations{j})], ...
                        ['Noise\' convertStringsToChars(all_modulations{i}) '+' convertStringsToChars(all_modulations{j})]);
                end
            end
        end
        deleteSingleComponentSignals('Noise');
        deleteSingleComponentSignals('Pure');
        deleteSingleComponentSignals('Dirty');
    otherwise
        return
end
end

function addSingleComponentSignals(i_source_folder, j_source_folder, destination_folder)
file_index = 1;
mkdir(destination_folder);
all_i_signals = dir(i_source_folder);
all_j_signals = dir(j_source_folder);
assert(length(all_i_signals) == length(all_j_signals), "两单分量信号样本个数不一致！")
for k = 3:length(all_i_signals)
    i_signal = load([all_i_signals(k).folder '\' all_i_signals(k).name]).('s');
    j_signal = load([all_j_signals(k).folder '\' all_j_signals(k).name]).('s');
    s = i_signal + j_signal;
    save([destination_folder '\' num2str(file_index) '.mat'],"s");
    file_index = file_index + 1;
end
end

function deleteSingleComponentSignals(folder)
all_folders = dir(folder);
for i = 3:length(all_folders)
    if ~contains(all_folders(i).name,"+")
        rmdir([all_folders(i).folder '\' all_folders(i).name], 's')
    end
end
end