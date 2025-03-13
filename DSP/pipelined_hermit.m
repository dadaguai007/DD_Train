function [s,e,w]=pipelined_hermit(input,dn,k,tap,ref,sps,step_len1,step_len2)

% K个pipeline
% 每组一个q值
q=zeros(1,k);
% 抽头向量
w =zeros(k,1);

% 数据的长度
L1=length(input);
% 信号输出长度
n = 2*floor(L1/sps/2);

% hermit多项式的输出长度
L=length(Build_Hermit_Input([zeros(1,tap),1],3));
% 每个管道块的抽头向量B
vector_b=zeros(1,L);

% 参考抽头 % 参考抽头应该存在两种方法进行定义
input = cat(2,input(ref+1:end),dn(end-ref+1:end));

% 输出
s = zeros(n, 1);


% 待输入区间的长度及向量确定
L2=k-1+tap;
x1=zeros(1,L2);


for idx = 1:n-1

    %一阶前馈输入
    x1 = cat(2,x1(sps+1:end),input(sps*idx-sps+1:1:sps*idx));
    % 每个管道进行一定的延迟
    for i=1:k
        input_x=input(tap*(i-1)+1:tap*i);
        input_hermit=cat(2,input_x,q(i));
        out_hermit(i,:)=Build_Hermit_Input(input_hermit,3);
        % 输出与抽头相乘
        alpha=out_hermit(i,:).*vector_b;
        q(i)=sum(alpha);
    end


    % 输出为单值(1*N N*1)。
    s(idx)=q*w;
    % 误差
    e(idx) = dn(idx) - s(idx);

    % 每个管道块抽头更新
    for i=1:k
        out(i,:)=w(i)*out_hermit(i,:);
    end
    % 分子求和
    term=sum(out);
    % 迭代式
    term_k=term/(norm(term));

    % 管道更新
    vector_b=vector_b+e(idx)*step_len1*term_k;

    % 抽头更新
    w = w +  e(idx) * step_len2 * q.';

end


end