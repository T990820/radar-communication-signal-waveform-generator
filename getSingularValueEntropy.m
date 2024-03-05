function getSingularValueEntropy(paths)
h = waitbar(0,'Initializing','name','Starting SVD Analysis ...');
index = 0;
file_nums = 0; % paths中全部文件的数量
for path = paths % 获取全部文件的数量并判断文件/文件夹是否存在
    if ~exist(path,'file')
        error(['文件/文件夹' path '不存在！']);
    end
    file_nums = file_nums + length(getAllFiles(path));
end
for path = paths
    if ~isfolder(path) % 获取path文件夹下所有文件的绝对路径
        all_file_paths = path;
    else
        all_file_paths = getAllFiles(path);
    end
    for i = 1:size(all_file_paths,1)
        file_path = all_file_paths{i};
        I = double(imread(file_path))/255;
        [~,S,~] = svd(I);
        singular_value = diag(S);
        sorted_singular_value = sort(singular_value,'descend');
        sigma = sorted_singular_value(1:20);
        singular_value_entropy = -sum(sigma/sum(sigma).*log(sigma/sum(sigma))/log(2));
        filename_layerout = strsplit(file_path,'.');
        save([filename_layerout{1} '.mat'],'singular_value_entropy');
        delete(file_path);
        index = index + 1;
        waitbar(index/file_nums,h,['generating ' num2str(index) '/' num2str(file_nums) ' singular value entropy']);
    end
end
end