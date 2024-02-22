function data = getParameterDict(filename)
[~, ~, data] = xlsread(filename);
for i = 1:size(data,1)
    if isnan(data{i,1})
        data{i,1} = data{i-1,1};
    end
end
end