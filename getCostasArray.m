% 名称：getCostasArray
% 作者：田秾恺
% 日期：2022/8/8 16:48
% 参数：dim(Costas序列长度)、arrs(所有可行的Costas序列构成的行向量组)
% 功能：获取所有长度为dim的Costas序列
% 注意：矩阵的离散自相关函数参考自余凡著《基于分布式架构的Costas序列搜索系统的设计与实现》P7
%      矩阵的离散自相关函数实例参考自姚建国著《Costas序列在雷达信号设计中的应用研究》P1
function [] = getCostasArray(dim)
arrs = perms(1:dim);
Costas  = [];
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
        Costas = [Costas;arrs(i,:)];
    end
end
save(['COSTAS n=' num2str(dim) '.mat'],'Costas')