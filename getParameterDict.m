% 读取信号参数.xlsx中定义的参数
function data = getParameterDict(filename)
[~, ~, data] = xlsread(filename);
for i = 1:size(data,1)
    if exist('program_name')~=0
        if regexp(program_name, 'octave' )==1
            if length(data{i,1}) == 0
                data{i,1} = data{i-1,1};
            end
        end
    else
        if isnan(data{i,1})
            data{i,1} = data{i-1,1};
        end
    end
end
end
