% Function: getDataset
% Author:   Nongkai Tian
% Date:     2023/4/27 00:46
% Instruments: Get all images from strPath and its' subfolders, eg. all_files = getAllFiles('训练集/');
% 参考文献：
%   [1] https://blog.csdn.net/qq_36320710/article/details/116742409
function [files] = getAllFiles(inputdir)
if ~exist(inputdir,'file')
    error(['文件/文件夹' convertStringsToChars(inputdir) '不存在！']);
end
files=cell(0);
dirs = dir(inputdir);
for i=1:length(dirs)
    if strcmp(dirs(i).name,'.')==1||strcmp(dirs(i).name,'..')==1
        continue;
    else
        if(dirs(i).isdir==1)
            files=[files;getAllFiles([dirs(i).folder,'\',dirs(i).name])];
        else
            files = [files;[dirs(i).folder,'\',dirs(i).name]];
        end
    end
end
end