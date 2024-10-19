%
% @brief        对folder文件内的时域信号进行时频分析
% @param[in]    folder 字符串，表示.mat文件格式的时域信号的存储路径
% @param[in]    tfa 函数句柄，表示时频分析所要用到的Matlab函数
% @param[in]    flag_delete_mat 标志位，表示时频分析结束之后要不要删除.mat文件
%
function [] = TimeFrequencyAnalysis(folder,tfa,shape,flag_delete_mat)
addpath("tftb-0.2\mfiles\");
h     = waitbar(0,'Initializing','name','Time-Frequency Analysis');
files = getAllFiles(folder);
param_dict = getParameterDict("信号参数.xlsx");
fs = getParameter(param_dict, '总体', 'fs');
type = getParameter(param_dict, '总体', 'type');
for i = 1:length(files)
    s = load(files{i}).('s');
    [p,~,~] = tfa(s');
    img = abs(p);
    img = imresize(img(end/2:end,:),shape);
    img = img / max(max(img));
    split_files = strsplit(files{i},'.');
    imwrite(img,[split_files{1} '.jpg']);
    waitbar(i/length(files),h,['正在生成第' num2str(i) '/' num2str(length(files)) '张STFT时频图像']);
end
if flag_delete_mat
    waitbar(1,h,'正在删除.mat文件');
    for i = 1:length(files)
        delete(files{i});
    end
end
waitbar(1,h,'完成');
close(h)
