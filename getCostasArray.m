% 名称：getCostasArray
% 作者：田秾恺
% 日期：2022/8/8 16:48
% 参数：dim(Costas序列长度)、arrs(所有可行的Costas序列构成的行向量组)
% 功能：获取所有长度为dim的Costas序列
% 注意：矩阵的离散自相关函数参考自余凡著《基于分布式架构的Costas序列搜索系统的设计与实现》P7
%      矩阵的离散自相关函数实例参考自姚建国著《Costas序列在雷达信号设计中的应用研究》P1
function [res] = getCostasArray(dim)
arrs = getFullPermutation(1:dim,dim);
res  = [];
for i = 1:size(arrs,1)
    matrix = zeros(dim,dim);
    for j = 1:dim
        matrix(dim-arrs(i,j)+1,j) = 1;
    end
    R = zeros(2*dim-1,2*dim-1);
    for j = 1-dim:dim-1     % 左右方向，与列数相关
        for k = 1-dim:dim-1 % 下上方向，与行数相关
            if j <=0 && k <= 0 % 左下
                for m = 1-k:size(matrix,1)     % 按行遍历
                    for n = 1:size(matrix,2)+j % 按列遍历
                        if matrix(m,n) == 1 && matrix(m+k,n-j) == 1
                            R(dim+j,dim+k) = R(dim+j,dim+k) + 1;
                        end
                    end
                end
            elseif j <=0 && k >= 0 % 左上
                for m = 1:size(matrix,1)-k
                    for n = 1:size(matrix,2)+j
                        if matrix(m,n) == 1 && matrix(m+k,n-j) == 1
                            R(dim+j,dim+k) = R(dim+j,dim+k) + 1;
                        end
                    end
                end
            elseif j >=0 && k <= 0 % 右下
                for m = 1-k:size(matrix,1)
                    for n = 1+j:size(matrix,2)
                        if matrix(m,n) == 1 && matrix(m+k,n-j) == 1
                            R(dim+j,dim+k) = R(dim+j,dim+k) + 1;
                        end
                    end
                end
            elseif j >=0 && k >= 0 % 右上
                for m = 1:size(matrix,1)-k
                    for n = 1+j:size(matrix,2)
                        if matrix(m,n) == 1 && matrix(m+k,n-j) == 1
                            R(dim+j,dim+k) = R(dim+j,dim+k) + 1;
                        end
                    end
                end
            end
        end
    end
    for j = 1:2*dim-1
        for k = 1:dim
            temp = R(j,k);
            R(j,k) = R(j,2*dim-k);
            R(j,2*dim-k) = temp;
        end
    end
    R = R';
    flag = true;
    for j = 1:2*dim-1
        for k = 1:2*dim-1
            if j ~= dim && k ~= dim && R(j,k) > 1
                flag = false;
            end
        end
    end
    if flag == true
        res = [res;arrs(i,:)];
    end
end

function [res] = getFullPermutation(V,N)
narginchk(2,2) ;
if isempty(V) || N == 0
    M = [] ;
elseif fix(N) ~= N || N < 1 || numel(N) ~= 1
    error('combn:negativeN','Second argument should be a positive integer') ;
elseif N==1
    M = V(:).';
else
    if nargout<2
        M = local_allcomb(V,N) ;
    else
        IND = local_allcomb(1:numel(V),N) ;
        M = V(IND) ;
    end
end
NewRows = size(M,2);
for i = 1:size(M,1)
    if length(M(i,:))-length(unique(M(i,:))) ~= 0
        M(i,:) = -1 * ones(1,size(M,2));
        NewRows = NewRows - 1;
    end
end
res = zeros(NewRows,size(M,2));
index = 1;res = [];
for i = 1:size(M,1)
    if sum(M(i,:) > 0)
        res(index,:) = M(i,:);
        index = index + 1;
    end
end

function Y = local_allcomb(X,N)
if N>1
    [Y{N:-1:1}] = ndgrid(X) ;
    Y = reshape(cat(N+1,Y{:}),[],N) ;
else
    Y = X(:);
end