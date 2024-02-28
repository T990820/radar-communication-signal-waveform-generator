% Instruments:
%   Mode
%       Mode(1) = 1: 生成盒维数
%       Mode(2) = 1: 生成信息维数
% 参考文献
%   [1] 曲志昱,毛校洁,侯长波.基于奇异值熵和分形维数的雷达信号识别[J].系统工程与电子技术,2018,40(02):303-307.
function getFractalDimension(mode, paths)
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
        fractal_dimension = zeros(1,length(mode));
        if mode(1) == 1
            s = load(file_path).('s');
            fractal_dimension(1) = getBoxDimension(s);
        end
        if mode(2) == 1
            s = load(file_path).('s');
            fractal_dimension(2) = getInformationDimension(s);
        end
        save(file_path,'fractal_dimension');
    end
end
end

function Db = getBoxDimension(X) % 参考自文献[1]的2.1节
N = length(X); % X是归一化后的fft模值，参考自我与毛校洁的微信聊天记录
delta = 1/N;
max_component = 0;
min_component = 0;
for i = 1:N-1
    max_component = max_component + max(X(i),X(i+1)) * delta;
    min_component = min_component + min(X(i),X(i+1)) * delta;
end
N_delta = N + (max_component - min_component) / delta / delta;
Db = -log(N_delta) / log(delta);
end

function DI = getInformationDimension(X) % 参考自文献[1]的2.2节
N = length(X);
Y = zeros(1,N-1);
for i = 1:N-1
    Y(i) = abs(X(i+1)-X(i));
end
W = sum(Y);
P = Y / W;
DI = sum(P.*log(P))/log(N);
end