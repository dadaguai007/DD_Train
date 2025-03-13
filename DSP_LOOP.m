classdef DSP_LOOP < handle

    properties
        % configuration of LOOP

        % 步长
        u=0.01;
        step=0.001;
        % 前馈抽头
        k1=31;
        k2=31;
        k3=31;

        % 反馈抽头
        d1=15;
        d2=15;
        d3=15;

        % 参考抽头，由第一个前馈抽头决定最大值
        ref=8;

        lamda=0.9999;

        Deocde_type='PAM';
        Config; % 存储图形的配置参数
        flagWhiteMarkerFace = 0; %控制标记的填充颜色
        % auto zoom configuration
        flagAutoZoom = 1;% 标志是否自动缩放图像。
        % minor tick
        flagMinorTick = 1;% 标志是否显示次要刻度
        % legend
        flagAddLegend = 1;% 标志是否添加图例


    end



    methods
        function obj = DSP_LOOP()
            % initialize the figure configuration
            Config.AxisLineWidth = 1.25;
            Config.AxisFontSize = 8;
            Config.LineWidth =  1.5;
            Config.MarkerSize = 8;
            Config.FontSize = 12;
            Config.FigureSize = [420 350];
            %Config.AxisLineWidth：坐标轴线的宽度。
            %Config.AxisFontSize：坐标轴上标签的字体大小。
            %Config.LineWidth：曲线的线宽。
            % Config.MarkerSize：标记的大小。
            % Config.FontSize：图像中的字体大小。
            % Config.FigureSize：绘图窗口的大小，以像素为单位。
            obj.Config = Config;
        end

        function signal=choose_eq(obj,xn,label,type)

            switch type
                case 'ffe_lms'
                    EQ.u=obj.u;
                    EQ.k1=obj.k1;
                    EQ.ref=obj.ref;
                    EQ.sps=2;
                    [signal,~,~] = FFE_LMS(EQ, xn, label);
                case 'ffe_rls'
                    EQ.u=obj.u;
                    EQ.k1=obj.k1;
                    EQ.ref=obj.ref;
                    EQ.lamda=obj.lamda;
                    EQ.sps=2;
                    [signal,~,~] = FFE_RLS(EQ, xn, label);
                case 'dfe_lms'
                    EQ.u=obj.u;
                    EQ.k1=obj.k1;
                    EQ.k2=obj.d1;
                    EQ.ref=obj.ref;
                    EQ.sps=2;
                    [signal,~,~] = DFE_LMS(EQ, xn, label);
                case 'dfe_rls'
                    EQ.u=obj.u;
                    EQ.k1=obj.k1;
                    EQ.k2=obj.d1;
                    EQ.ref=obj.ref;
                    EQ.lamda=obj.lamda;
                    EQ.sps=2;
                    [signal,~,~] = DFE_RLS(EQ, xn, label);

                case 'volterra_ffe_lms'
                    taps_list=[obj.k1,obj.k2,obj.k3,0,0,0];
                    [signal,~,~]=volterra_dfe_lms(xn,label,2,obj.ref,taps_list,obj.u);
                case 'volterra_dfe_lms'
                    taps_list=[obj.k1,obj.k2,obj.k3,obj.d1,0,0];
                    [signal,~,~]=volterra_dfe_lms(xn,label,2,obj.ref,taps_list,obj.u);
                case 'volterra_dfe_rls'
                    taps_list=[obj.k1,obj.k2,obj.k3,obj.d1,0,0];
                    [signal,~,~]=volterra_dfe_rls(xn,label,2,obj.ref,taps_list,obj.lamda);
            end


        end

        function ber=NRZ_Decode(obj,sigRx_E,ref_seq)
            offset=0;
            % normalize
            sigRx_E = sigRx_E - mean(sigRx_E);
            sigRx_E = sigRx_E./max(sigRx_E);
            % make decision and convert into 0,1 sequence
            out = sign(sigRx_E-offset);
            out = 0.5*(out+1);

            ncut=1e5;
            [ber,~,~] = CalcBER(out(ncut:end),ref_seq(ncut:end)); %计算误码率
        end

        function ber=PAM_Decode(obj,sigRx_E,label)
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
            [ber,~,~] = CalcBER(I(ncut:end),label1(ncut:end)); %计算误码率
        end

        % 适配步长循环方式
        % 一般选取参考抽头为中心处
        function Loop_step(obj,xn,label,type)
            Axis=0:obj.step:obj.u;
            for index=1:length(Axis)
                u_axis=Axis(index);
                obj.u=u_axis;
                signal=choose_eq(obj,xn,label,type);

                switch obj.Deocde_type
                    case 'NRZ'
                        ber=NRZ_Decode(obj,signal,label);
                    case 'PAM'
                        ber=PAM_Decode(obj,signal,label);
                end
                BER(index)=ber;
            end
            %绘图
            figure(90)
            plot(Axis,BER,'Linewidth',obj.Config.LineWidth);
            set(gcf,'Position', obj.Config.FigureSize);

        end

        % 适配DFE，FFE抽头循环方式

        function M=Loop_one(obj,xn,label,type)
            % 存储元胞组
            M={};
            k_size=obj.k1;
            for k=1:k_size
                for Ref=1:k
                    obj.ref=Ref;
                    obj.k1=k;
                    signal=choose_eq(obj,xn,label,type);
                    switch obj.Deocde_type
                        case 'NRZ'
                            ber=NRZ_Decode(obj,signal,label);
                        case 'PAM'
                            ber=PAM_Decode(obj,signal,label);
                    end
                    % 存储BER
                    M{k,Ref}=ber;
                end
            end
        end


        function N=Loop_fed_one(obj,xn,label,type)
            N={};
            d_size=obj.d1;
            for k=1:d_size
                obj.d1=k;
                M=Loop_one(obj,xn,label,type);
                % 存储BER
                N{k}=M;
            end

        end

        % Volterra的循环方式

        function N=Loop_two(obj,xn,label,type)
            N={};
            k2_size=obj.k2;
            for k=1:k2_size
                obj.k2=k;
                M=Loop_one(obj,xn,label,type);
                % 存储BER
                N{k}=M;
            end

        end

    end





end