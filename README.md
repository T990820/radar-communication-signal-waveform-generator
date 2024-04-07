# 雷达/通信信号波形发生与特征提取代码
&emsp;&emsp;本仓库的代码用于雷达/通信信号波形的仿真与特征提取，各个函数的用法如下所示。
## getDataset
&emsp;&emsp;getDataset函数用于生成特定信噪比范围内的时域采样信号，有如下的2种调用形式。
```matlab
getDataset(Mode, Modulations, SignalType, SNRmin, SNRmax, SNRstep, Samples, snrth);
getDataset(Mode, Modulations, SignalType, SNR, Samples, snrth);
```
各个参数的实际含义如下所示。
### Mode
| __Mode__ | __意义__ |
|:---------:|:---------|
| 1 | 为时频图像识别网络的训练集生成时域信号，并将每个类的样本存到单独的文件夹里 |
| 1.1 | 为时频图像识别网络的训练集生成时域信号，并将所有样本存到一个文件夹里，样本文件命名为'BPSK_0005.mat' |
| 2 | 为时频图像识别网络的测试集生成时域信号，将每dB不同调制方式的样本存到不同的文件夹里 |
| 2.1 | 为时频图像识别网络的测试集生成时域信号，将每dB不同调制方式的样本存到一个文件夹里，样本文件命名为'BPSK_0005.mat' |
| 3 | 为时频图像去噪网络的训练集生成时域信号 |
| 4 | 为时频图像去噪网络的测试集生成时域信号 |
| 5 | 为噪声估计网络的训练集生成时域信号 |
### Modulations
&emsp;&emsp;Modulations是一个长度为23的向量，其中的每个元素代表一种调制方式，1表示生成该种调制方式的信号，0表示不生成该种调制方式的信号。
| __Modulations分量__ | __意义__ |
|:---------:|:---------:|
| Modulations(1) | 单音信号(NS) |
| Modulations(2) | 偶二次调频信号(EQFM) |
| Modulations(3) | 线性调频信号(LFM) |
| Modulations(4) | 双线性调频信号(VTFM) |
| Modulations(5) | 二进制相移键控信号(BPSK) |
| Modulations(6) | 四进制相移键控信号(QPSK) |
| Modulations(7) | 二进制频移键控信号(QFSK) |
| Modulations(8) | 正弦调频信号(SFM) |
| Modulations(9) | COSTAS |
| Modulations(10) | 正切调频信号(TANFM) |
| Modulations(11) | FRANK多相编码信号(FRANK) |
| Modulations(12) | P1多相编码信号(P1) |
| Modulations(13) | P1多相编码信号(P2) |
| Modulations(14) | P1多相编码信号(P3) |
| Modulations(15) | P1多相编码信号(P4) |
| Modulations(16) | T1多时编码信号(T1) |
| Modulations(17) | T2多时编码信号(T2) |
| Modulations(18) | T3多时编码信号(T3) |
| Modulations(19) | T4多时编码信号(T4) |
| Modulations(20) | LFM与SFM混合调制信号(LFM-SFM) |
| Modulations(21) | EQFM与SFM混合调制信号(EQFM-SFM) |
| Modulations(22) | 下变频后的多进制幅移键控信号(DDC-MASK) |
| Modulations(23) | 下变频后的多进制相移键控信号(DDC-MPSK) |
| Modulations(24) | 下变频后的多进制频移键控信号(DDC-MFSK) |
| Modulations(25) | 多进制正交幅度调制(MQAM) |
### SignalType
&emsp;&emsp;SignalType=1表示生成实信号，SignalType=2表示生成复指数信号
### SNRmin
&emsp;&emsp;信噪比的最小值
### SNRmax
&emsp;&emsp;信噪比的最大值
### SNRstep
&emsp;&emsp;信噪比的步进
### Samples
&emsp;&emsp;Samples是一个长度为23的向量，其中的每个元素代表一种调制方式每dB信噪比下生成的样本个数。Samples中23个分量与实际调制信号的对应关系与Modulations的一致。
### SNR
&emsp;&emsp;数据集中包含的信噪比数组，一般用于信噪比非均匀采样的情况下，比如SNR=[-14 -12 -10 -6 -4 -2 0 2]
### snrth
&emsp;&emsp;信噪比门限，仅用于为噪声估计网络生成训练集/测试集。
