function getMultiComponentSignal(Mode,Modulations)
%   Mode = 7: 为复合调制信号识别任务的训练集生成时域信号，并将每个类的样本存到单独的文件夹中
all_modulations = {"NS", "EQFM", "LFM", "VTFM", "BPSK", "QPSK", "BFSK", "SFM", "COSTAS", "TANFM", ...
                   "FRANK", "P1", "P2", "P3", "P4", "T1", "T2", "T3", "T4", "LFM-SFM", ...
                   "EQFM-SFM", "DDC-MASK", "DDC-MPSK", "DDC-MFSK", "MQAM"};
switch Mode
    case 7
        for i = 1 : length(Modulations)
            for j = i+1 : length(Modulations)
                if Modulations(i) == 1 && Modulations(j) == 1
                    file_index = 1;
                    mkdir(['训练集\' convertStringsToChars(all_modulations{i}) '+' convertStringsToChars(all_modulations{j})]);
                    all_i_signals = dir(['训练集\' convertStringsToChars(all_modulations{i})]);
                    all_j_signals = dir(['训练集\' convertStringsToChars(all_modulations{j})]);
                    assert(length(all_i_signals) == length(all_j_signals), "两单分量信号样本个数不一致！")
                    for k = 3:length(all_i_signals)
                        i_signal = load([all_i_signals(k).folder '\' all_i_signals(k).name]).('s');
                        j_signal = load([all_j_signals(k).folder '\' all_j_signals(k).name]).('s');
                        s = i_signal + j_signal;
                        save(['训练集\' convertStringsToChars(all_modulations{i}) '+' convertStringsToChars(all_modulations{j}) '\' num2str(file_index) '.mat'],"s");
                        file_index = file_index + 1;
                    end
                end
            end
        end
    otherwise
        return
end
end