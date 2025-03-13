function [res] = Build_Hermit_Input(input_X,order)
% 构建 指数 输入
% input_X: 输入信号向量
% order: Hermit 级数
% 非线性系数中不存在自拍频项

res = [];
u=1;% 二阶一次交叉项系数
uu=1;% 三阶一次交叉项系数

vz=1;% 三阶二，一次交叉项系数
zv=1;% 三阶一，二次交叉项系数

z=1;% 第一阶的系数
v=1;% 第二阶的系数
l=1;% 第三阶的系数
if order == 2
    % 第0阶
    res=[res 1];
    % 第1阶
    for i = 1:length(input_X)
        res = [res 2*input_X(i) ];
        index_linear1(z, :) = i;
        z=z+1;
    end
    % 第2阶
    for i = 1:length(input_X)
        res = [res 4*input_X(i).^2-2];
        index_nonlinear2(v, :) = i;
        v=v+1;
    end



    % 第2阶一次的交叉项
    for i = 1:length(input_X)-1
        for j = i+1:length(input_X)
            res = [res 4*input_X(i)*input_X(j)];
            index_mat_nonlinear1(u, :) = [i, j];
            u=u+1;
        end
    end

elseif order == 3

    % 第0阶
    res=[res 1];
    % 第1阶
    for i = 1:length(input_X)
        res = [res 2*input_X(i) ];
        index_linear1(z, :) = i;
        z=z+1;
    end
    % 第2阶
    for i = 1:length(input_X)
        res = [res 4*input_X(i).^2-2];
        index_nonlinear2(v, :) = i;
        v=v+1;
    end

    % 第3阶
    for i = 1:length(input_X)
        res = [res 8*input_X(i).^3-12*input_X(i)];
        index_nonlinear3(l, :) = i;
        l=l+1;
    end

    % 第3阶二,一次的交叉项
    for i = 1:length(input_X)
        for j = i+1:length(input_X)
            res = [res (4*input_X(i).^2-2)*(2*input_X(j))];
            index_mat_nonlinear2(vz, :) = [i, j];
            vz=vz+1;
        end
    end

    % 第3阶一,二次的交叉项
    for i = 1:length(input_X)
        for j = i+1:length(input_X)
            res = [res (4*input_X(j).^2-2)*(2*input_X(i))];
            index_mat_nonlinear2(zv, :) = [i, j];
            zv=zv+1;
        end
    end


    % 第3阶的一次项的交叉项
    for i = 1:length(input_X)
        for j = i+1:length(input_X)
            for m = j+1:length(input_X)
                res = [res 8*input_X(i) * input_X(j) * input_X(m)];
                index_mat_nonlinear3(uu, :) = [i, j,m];
                uu=uu+1;
            end
        end
    end
end

end





