function [] = TimeFrequencyAnalysis(folder,mode,Shape,delete_mat)
% Instruments:
%   folder
%       待处理的文件夹，例如'训练集'
%   mode
%       mode = 1: 对folder内的雷达时域信号做STFT
%       mode = 2: 求folder内的雷达时域信号的CWD
%       mode = 3: 求folder内的雷达时域信号的CTFD
h     = waitbar(0,'Initializing','name','Time-Frequency Analysis');
files = getAllFiles(folder);
param_dict = getParameterDict("信号参数.xlsx");
fs = getParameter(param_dict, '总体', 'fs');
type = getParameter(param_dict, '总体', 'type');
if mode == 1
    time_resolution = 128/fs;
    overlap = 90;
    for i = 1:length(files)
        s = load(files{i}).('s');
        % [p,~,~]=pspectrum(s,fs,'spectrogram','TimeResolution',time_resolution,'Overlap',overlap,'Leakage',0.1);
        % if type == 1
        %     img=flipud(imresize(p,Shape));
        % else
        %     img=flipud(imresize(p(end/2:end,:),Shape));
        % end
        [p,~,~]=tfrstft(s');
        if type == 2
            img = abs(p);
            img = imresize(img(end/2:end,:),Shape);
            img = img / max(max(img));
        end
        split_files = strsplit(files{i},'.');
        imwrite(img,[split_files{1} '.jpg']);
        waitbar(i/length(files),h,['正在生成第' num2str(i) '/' num2str(length(files)) '张STFT时频图像']);
    end
elseif mode == 2 % 基于tfrcw()的CWD计算过程很慢
    for i = 1:length(files)
        s = load(files{i}).('s');
        [img,~,~] = tfrcw(s'); % tfrcw只能分析列向量形式的信号
        img = imresize(img,Shape);
        split_files = strsplit(files{i},'.');
        imwrite(img/max(max(img)),[split_files{1} '.jpg']);
        waitbar(i/length(files),h,['正在生成第' num2str(i) '/' num2str(length(files)) '张CWD时频图像']);
    end
elseif mode == 3
    for i = 1:length(files)
        s = load(files{i}).('s');
        img = CTFD(s);
        split_files = strsplit(files{i},'.');
        imwrite(img,[split_files{1} '.jpg']);
        waitbar(i/length(files),h,['正在生成第' num2str(i) '/' num2str(length(files)) '张CTFD时频图像']);
    end
end
if delete_mat
    waitbar(1,h,'正在删除.mat文件');
    for i = 1:length(files)
        delete(files{i});
    end
end
waitbar(1,h,'完成');
close(h)