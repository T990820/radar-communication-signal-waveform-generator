% Function: getDataset
% Author:   Nongkai Tian
% Date:     2024/1/17 12:13
% Instruments: 在getDataset.m中调用，用来获取要生成的样本的个数
function [res] = getSampleNums(Modulations, Samples, SNR, M_for_DDC_MASK, M_for_DDC_MPSK, M_for_DDC_MFSK)
    assert(length(Modulations)==length(Samples), "表示包含的调制方式的向量长度必须要等于表示各调制方式每dB信噪比下生成的样本数的向量长度！")
    res = abs(Modulations)*Samples'*length(SNR);
    if ~isempty(M_for_DDC_MASK)
        res = res + Samples(21)*(length(M_for_DDC_MASK)-1)*length(SNR);
    end
    if ~isempty(M_for_DDC_MPSK)
        res = res + Samples(22)*(length(M_for_DDC_MPSK)-1)*length(SNR);
    end
    if ~isempty(M_for_DDC_MFSK)
        res = res + Samples(23)*(length(M_for_DDC_MFSK)-1)*length(SNR);
    end
end