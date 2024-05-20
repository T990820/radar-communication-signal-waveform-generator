% 读取信号参数.xlsx中
function param = getParameter(dict, modulation, key)
for i = 1:size(dict,1)
    if strcmp(dict{i,1}, modulation) && strcmp(dict{i,2}, key)
        param = dict{i,3};
        return
    end
end
error("索引modulation或key不存在！");
end

