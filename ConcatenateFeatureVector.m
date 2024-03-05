function ConcatenateFeatureVector(paths,destination_folder)
%% 判断paths中若干文件夹内文件的数量是否相等，如不相等则报错退出
file_nums = [];
for path = paths
    if exist(path)
        all_files_folders = dir(path); % 判断path下有无文件夹
        for i = 3:length(all_files_folders)
            if all_files_folders(i).isdir == 1
                error([convertStringsToChars(path) '中存在文件夹！']);
            end
        end
        file_nums = [file_nums length(dir(path))]; % 统计path路径下文件的个数
    else
        error([convertStringsToChars(path) '不存在！']);
    end
end
if ~isequal(file_nums,file_nums(1)*ones(1,length(file_nums)))
    error('paths中若干文件夹内文件个数不相等');
end
all_files = getAllFiles(path); % 获取所有.mat文件的文件名
for i = 1:size(all_files,1)
    detail_path = strsplit(all_files{i},'\'); % 获取all_files中第i个文件的详细路径
    filename = detail_path(end);
    all_files{i} = filename;
end
mkdir(destination_folder);
for i = 1:length(all_files)
    s = [];
    for path = paths
        feature_vector = load([convertStringsToChars(path) '\' cell2mat(all_files{i})]);
        fields = fieldnames(feature_vector); % 获取结构体中所有的字段名称
        s = [s load([convertStringsToChars(path) '\' cell2mat(all_files{i})]).(cell2mat(fields(1)))];
    end
    save([convertStringsToChars(destination_folder) '\' cell2mat(all_files{i})],'s');
end
end