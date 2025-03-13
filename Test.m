clc;clear;close all;
addpath("D:\001-处理中\相干算法\optical_communication\")
addpath("Data\SOA1\")
addpath("Fncs\")
addpath("DSP\")

load('ref_d2.mat')
load('Data\SOA1\20240526_pam4_28G_150mA_SOAin_-2_btb\ROP--10.0.mat')

fs = 80e9; % sampling rate
sps = 2; % sample per symbol in dsp
osr = 4; % oversampling rate (compared with baud rate)

fb = 28e9; % 信号的波特率 baud rate

% 参考序列（重要）
ref=xx;
ref2 = repelem(ref,sps);
ref_seq = repmat(ref,1,1000);
ref_seq = ref_seq(:);
label=ref_seq;

%数字时钟恢复
data = resample(data,osr*fb,fs);  %四倍上采样，为了满足后续时钟恢复频率
f_clk_rec = cr(data.',fb,osr,0); % 注意：f_clk_rec是对fb的估计，不是对fb*osr的估计！
t_cur = 1/osr/fb*(0:length(data)-1); % 注意：插值前data1的采样率是osr*fb
t_new = 0:1/f_clk_rec/osr:t_cur(end); % 注意：插值后data2的采样率是osr*f_clk_rec
data1 = interp1(t_cur,data,t_new,'spline');
data1 = resample(data1,sps,osr);
data2 = std(ref)*normalization(data1.');
x=pnorm(data1);
% 同步
[xn,ref2,ff] = sync(x,ref2);


EQ=struct();
EQ.u=0.0015;
EQ.k1=31;
EQ.k2=15;
EQ.ref=8;
EQ.sps=2;
EQ.lamda=0.9999;

% FFE_LMS
[sigRx_E,en,w] = FFE_LMS(EQ, xn, label);

figure;hold on
plot(xn(1:1e5),'.')
plot(sigRx_E(1:1e5),'k.')
plot(label(1:1e5),'.')
legend('接收信号','均衡后信号','发送信号')

figure
stem(w(:))
title("均衡器抽头")

figure;
semilogy(abs(en(1:1e4)).^2)
xlabel("迭代次数")
ylabel("误差")



% 解码，计算误码
% 重新量化
A1=[-2 0 2];
% 参考序列
[~,label1] = quantiz(label,A1,[-3,-1,1,3]);
label1=decoding(label1,'gray');
% 接收序列
[~,I] = quantiz(sigRx_E,A1,[-3,-1,1,3]);
I=decoding(I,'gray');
ncut=1e5;
[ber,num,error_location] = CalcBER(I(ncut:end),label1(ncut:end)); %计算误码率
fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);