function [res] = Build_Exponen_Input(input_X,q,alpha,mode)
% 构建 指数 输入
% input_X: 输入信号向量
%  指数输入，没有受阶数的影响
% 简化方式(mode)： 输入 * 指数（输入）
% q：裁剪延时拍频项：自拍频，延迟过长的项不需要。(q为数)
% 需在外部函数处进行 自适应参数 alpha的更新
if nargin < 4
    mode='exponen_add';
end

res = [];
u=1;

if strcmp(mode,'exponen_add')

    for i = 1:q
        for j = 1:length(input_X)-i
            res = [res exp(alpha*(input_X(j+i)+input_X(j)))];
            index_mat_nonlinear1(u, :) = [j, j+i];
            u=u+1;
        end
    end

elseif strcmp(mode,'exponen')

    for i = 1:length(input_X)
        res = [res exp(alpha*input_X(i))];
        index_mat_nonlinear1(u) = i;
        u=u+1;
    end

elseif strcmp(mode,'exponen_multiply')

    for i = 1:q
        for j = 1:length(input_X)-i
            res = [res input_X(j+i)*exp(alpha*(input_X(j)))];
            index_mat_nonlinear1(u, :) = [j, j+i];
            u=u+1;
        end
    end
elseif strcmp(mode,'exponen_multiply_add')
    % 用于更新系数
    % x1*x2*exp(x1*x2)
    for i = 1:q
        for j = 1:length(input_X)-i
            res = [res (input_X(j+i)+input_X(j))*exp(alpha*(input_X(j)+input_X(j+i)))];
            index_mat_nonlinear1(u, :) = [j, j+i];
            u=u+1;
        end
    end

end

end


