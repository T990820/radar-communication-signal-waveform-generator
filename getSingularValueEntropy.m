function getSingularValueEntropy(paths)
for path = paths
    if ~exist(path,'file')
        error(['文件/文件夹' path '不存在！']);
    end
    if ~isfolder(path)
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
    end
end
end