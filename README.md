# 雷达/通信信号波形发生与时频分析代码
&emsp;&emsp;本仓库的代码用于雷达/通信信号波形的仿真与时频分析，各个函数的用法如下所示。
## getDataset
&emsp;&emsp;getDataset函数有如下的2种调用形式。
```
getDataset(Mode, Modulations, SignalType, SNRmin, SNRmax, SNRstep, Samples, snrth);
getDataset(Mode, Modulations, SignalType, SNR, Samples, snrth);
```
各个参数的实际含义如下所示。
### Mode
Mode     | 意义
-------- | -----
1  | 为时频图像识别网络的训练集生成时域信号，并将每个类的样本存到单独的文件夹里
手机  | $12
导管  | $1
### SNRmin
&emsp;&emsp;信噪比的最小值
### SNRmax
&emsp;&emsp;信噪比的最大值
