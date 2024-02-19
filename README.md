# 雷达/通信信号波形发生与时频分析代码
&emsp;&emsp;本仓库的代码用于雷达/通信信号波形的仿真与时频分析，各个函数的用法如下所示。
## getDataset
&emsp;&emsp;getDataset函数有如下的2种调用形式。
```matlab
getDataset(Mode, Modulations, SignalType, SNRmin, SNRmax, SNRstep, Samples, snrth);
getDataset(Mode, Modulations, SignalType, SNR, Samples, snrth);
```
各个参数的实际含义如下所示。
### Mode
| __Mode__ | __意义__ |
|:---------:|:---------:|
| 1 | 为时频图像识别网络的训练集生成时域信号，并将每个类的样本存到单独的文件夹里 |
| 1.1 | 为时频图像识别网络的训练集生成时域信号，并将所有样本存到一个文件夹里，样本文件命名为'BPSK_0005.mat' |
| 2 | 为时频图像识别网络的测试集生成时域信号，将每dB不同调制方式的样本存到不同的文件夹里 |
| 2.1 | 为时频图像识别网络的测试集生成时域信号，将每dB不同调制方式的样本存到一个文件夹里，样本文件命名为'BPSK_0005.mat' |
| 3 | 为时频图像去噪网络的训练集生成时域信号 |
| 4 | 为时频图像去噪网络的测试集生成时域信号 |
| 5 | 为噪声估计网络的训练集生成时域信号 |
### Modulations
&emsp;&emsp;Modulations是一个长度为23的向量，其中的每个元素代表一种调制方式，1表示生成，0表示不生成。
| Modulations(1) | 单音信号(NS) |
|:---------:|:---------:|
| Modulations(2) | 偶二次调频信号(EQFM) |
### SNRmin
&emsp;&emsp;信噪比的最小值
### SNRmax
&emsp;&emsp;信噪比的最大值
