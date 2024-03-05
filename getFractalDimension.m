% Instruments:
%   Mode
%       Mode(1) = 1: 生成盒维数
%       Mode(2) = 1: 生成信息维数
% 参考文献
%   [1] 曲志昱,毛校洁,侯长波.基于奇异值熵和分形维数的雷达信号识别[J].系统工程与电子技术,2018,40(02):303-307.
function getFractalDimension(mode, paths)
h = waitbar(0,'Initializing','name','分形维数分析');
global_index = 1;
samples_num = 0; % 总样本数量
for path = paths
    all_file_paths = getAllFiles(path);
    samples_num = samples_num + size(all_file_paths,1);
end
for path = paths
    all_file_paths = getAllFiles(path);
    for i = 1:size(all_file_paths,1)
        waitbar(global_index/samples_num,h,['正在生成第' num2str(i) '/' num2str(samples_num) '个样本的分形维数']);
        file_path = all_file_paths{i};
        fractal_dimension = zeros(1,length(mode));
        if mode(1) == 1
            s = load(file_path).('s');
            fractal_dimension(1) = getBoxDimension(s);
            if isnan(fractal_dimension(1))
                error([file_path '计算出的盒维数是NAN']);
            end
        end
        if mode(2) == 1
            s = load(file_path).('s');
            fractal_dimension(2) = getInformationDimension(s);
            if isnan(fractal_dimension(2))
                error([file_path '计算出的信息维数是NAN']);
            end
        end
        save(file_path,'fractal_dimension');
        global_index = global_index + 1;
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
P(P<1e-6)=1e-6; % 避免出现0×log(0)的情况，从而使结果为NaN
DI = sum(P.*log(P))/log(N);
end