clc;clear;close all;
current_date = date;
disp(current_date);

addpath("D:\001-处理中\相干算法\optical_communication\")
addpath("D:\PhD\IM_DD_DSP\Data\DML-NRZ-50G\")
addpath("D:\PhD\IM_DD_DSP\Data\DML-NRZ-50G\252mA")
addpath("Fncs\")
addpath("DSP\")

% 读取参考序列
load ref32768 % can be replaced by prbs.m
ref = ref32768;

fs = 80e9; % sampling rate
sps = 2; % sample per symbol in dsp
osr = 4; % oversampling rate (compared with baud rate)
fb = 50e9; % 信号的波特率 baud rate

% 参考序列
ref2 = repelem(ref,sps);
ref_seq = repmat(ref,1,1000);
ref_seq = ref_seq(:);
label=ref_seq;

% 解码为0,1码，用于解误码
ref_seq=0.5*(ref_seq+1);
% 文件存储
datapath='Output\DML_NRZ_50G_BTB_ffe';
if ~exist(datapath,'dir')
    mkdir(datapath);
end

% 文件存储
% Vpp = [4.0,5.1,6.1,7.1,8.1,9.1,10.1,11.1,12.0,13.0,14.0,15.0];
Vpp=4:1:15;
Vpp=flip(Vpp);
BER=zeros(length(Vpp),1);
V_PP=zeros(length(Vpp),1);

foldername = '252mA';

for powID = 1:length(Vpp)
    BER_NUM=zeros(length(Vpp),1);
    for fileID = 1:1
        filename = sprintf('%s\\12dB-at252mA@NRZ_0km_-%1.1fdBm-PRBS215-ch4-50G-%d.mat',...
            foldername,Vpp(powID),fileID);
        load(filename);
        data = data';
        %数字时钟恢复
        Time_recovery;
        % 均衡
        PAM_Equ;
        % 解码
        NRZ_Decode;
        BER_NUM(fileID)=ber;
        close all;
    end
    
    ber=mean(BER_NUM);
    % 数据存储
    name=strcat('BER_',num2str(Vpp(powID)),'Vpp');
    %     save(sprintf('%s\\%s.mat',datapath,name),'ber');
    BER(powID)=ber;
    V_PP(powID)=-Vpp(powID);
end


berplot = BERPlot_David();

berplot.interval=2;
berplot.flagThreshold=1;
berplot.flagRedraw=1;
berplot.flagAddLegend=0;
berplot.orderInterpolation=1;
berplot.plot(V_PP,BER);