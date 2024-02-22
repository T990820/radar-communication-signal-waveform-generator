% Function: getDataset
% Author:   Nongkai Tian
% Date:     2024/1/25 19:32
% Instruments:
%   对生成的时域信号进行预处理
%   mode
%       mode(1) = 1：将paths下所有文件夹内的.mat文件格式的时域信号文件转为.npy文件格式。如果.mat文件存储的是复数则先取其幅值再保存成.npy文件
%       mode(2) = 1：将paths中所有文件夹中的文件按照相同的策略去打乱重排，重排之后的文件名格式是"数字索引_调制方式"
%       mode(3) = 1：将paths下所有文件夹内的灰度图像转为彩色图像
%   file_path_length：文件名中数字索引部分的长度，应当等于每个文件夹中文件数量的位数
%   paths：待处理的所有文件夹路径，这些文件夹应包含相同数量的文件
function [] = Preproccess(mode,paths,file_path_length)
if mode(1) == 1
    h = waitbar(0,'Initializing','name','Generating .npy files ...');
    index = 0;
    file_nums = 0;
    for path = paths
        file_nums = file_nums + length(getAllFiles(path)); % 获取全部文件的数量
    end
    for path = paths
        all_files = getAllFiles(path);
        for i = 1:length(all_files)
            splited_filename = strsplit(all_files{i},'.'); % 将文件名以'.'为分割符进行分解，得到的元胞数组中第2个元素是文件类型
            if strcmp(splited_filename{2},"mat")
                s = load(all_files{i});
                s = s.s;
                if ~isreal(s)
                    s = abs(s);
                end
                writeNPY(s,[splited_filename{1} '.npy']);
                delete(all_files{i});
            end
            index = index + 1;
            waitbar(index/file_nums,h,['generating ' num2str(index) '/' num2str(file_nums) ' .npy files']);
        end
    end
    delete(h);close all;
end
if mode(2) == 1
    %% 判断paths中若干文件夹内文件的数量是否相等，如不相等则报错退出
    file_nums = [];
    for path = paths
        file_nums = [file_nums length(dir(path))];
    end
    if ~isequal(file_nums,file_nums(1)*ones(1,length(file_nums)))
        error('paths中若干文件夹内文件个数不相等');
    end
    %% 按相同的策略随机打乱每个文件夹内的文件
    file_index = randperm(length(dir(paths(1)))-2);
    for path = paths
        all_files = dir(path);
        for i = 3:length(all_files)
            splited_filename1 = strsplit(all_files(i).name,'_'); % 将文件名以'_'为分割符进行分解，得到的元胞数组中第一个元素是调制方式，第二个是索引
            splited_filename2 = strsplit(all_files(i).name,'.'); % 将文件名以'.'为分割符进行分解，得到的元胞数组中第2个元素是文件类型
            file_name = sprintf(['%0' num2str(file_path_length) 's'],num2str(file_index(i-2))); % 新文件名的数字部分。应该是数字部分在前，调制方式在后
            source_path = [all_files(i).folder '\' all_files(i).name];
            target_path = [all_files(i).folder '\' file_name '_' splited_filename1{1} '.' splited_filename2{2}];
            copyfile(source_path, target_path);
            delete([all_files(i).folder '\' all_files(i).name]);
        end
    end
end
if mode(3) == 1
    for path = paths
        all_files = getAllFiles(path);
        for i = 1:length(all_files)
            I = imread(all_files{i});
            C = ind2rgb(I,parula(256));
            imwrite(C, all_files{i});
        end
    end
end