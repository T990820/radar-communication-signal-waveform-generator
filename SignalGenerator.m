classdef SignalGenerator
    properties(Access = private)
        fs        % 采样频率
        Points    % 采样点数
        TimeWidth % 采样时长
        CodeWidth % 码元宽度
        V         % 电平
        P         % 功率
    end
    methods(Access = public)
        function obj = SignalGenerator(fs,V,n,TimeWidth) % 构造函数
            obj.fs        = fs;        % sampling frequency
            obj.V         = V;         % 信号电平
            obj.Points    = n;         % 默认采样点数
            obj.TimeWidth = TimeWidth; % 时域宽度
        end
        function [s,noise] = generateNS(obj,time_width, n, f0,snr,SignalType)
            t = linspace(0,time_width,n);
            if SignalType == 1
                s = obj.V*cos(2*pi*f0*t);
            else
                s = obj.V*cos(2*pi*f0*t);
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateEQFM(obj,f_valley,n,band_width,time_width,snr,SignalType)
            k = 8*band_width/3/time_width/time_width;
            t = linspace(0,time_width,n);
            T = time_width/2;
            if SignalType == 1
                s = obj.V*cos(2*pi*f_valley*t+pi*k*(t-T).^3);
            else
                s = obj.V*cos(2*pi*f_valley*t+pi*k*(t-T).^3);
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateLFM(obj,snr,f0,band_width,n,SignalType)
            time_width = n / obj.fs; % 时域宽度
            t = linspace(0,time_width,n);
            if SignalType == 1
                s = obj.V*cos(2*pi*f0*t+pi*(band_width/time_width)*t.*t);
            else
                s = obj.V*cos(2*pi*f0*t+pi*(band_width/time_width)*t.*t);
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateVTFM(obj,snr,f0,band_width,SignalType)
            % 频率上升沿
            time_width = obj.Points / obj.fs / 2; % 时域宽度
            t = linspace(0, time_width, obj.Points / 2);
            if SignalType == 1
                s1 = obj.V*cos(2*pi*f0*t+pi*(band_width/time_width)*t.*t); % 频率上升沿的时域信号
            else
                s1 = obj.V*cos(2*pi*f0*t+pi*(band_width/time_width)*t.*t);
                s1 = HilbertTransfer(s1)/sqrt(2);
            end
            fmax = f0 + band_width;
            % 频率下降沿
            if SignalType == 1
                s2 = obj.V*cos(2*pi*fmax*t-pi*(band_width/time_width)*t.*t); % 频率上升沿的时域信号
            else
                s2 = obj.V*cos(2*pi*fmax*t-pi*(band_width/time_width)*t.*t);
                s2 = HilbertTransfer(s2)/sqrt(2);
            end
            s = [s1 s2];
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateBPSK(obj,Code,f,t,snr,SignalType)
            s = [];
            for j = 1:length(Code)
                if Code(j) == 0
                    if SignalType == 1
                        s=[s obj.V*cos(2*pi*f*t-pi/3)];
                    else
                        s=[s HilbertTransfer(obj.V*cos(2*pi*f*t-pi/3))/sqrt(2)];
                    end
                else
                    if SignalType == 1
                        s=[s obj.V*cos(2*pi*f*t+2*pi/3)];
                    else
                        s=[s HilbertTransfer(obj.V*cos(2*pi*f*t+2*pi/3))/sqrt(2)];
                    end
                end
            end
            s = obj.zero_fill_or_transaction(s);
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateQPSK(obj,snr,Code,SignalType)
            f = (0.1+0.1*rand(1))*obj.fs;
            code_width = obj.TimeWidth / length(Code) * 2;
            t = linspace(0,obj.TimeWidth/(obj.Points/obj.fs/code_width),code_width*obj.fs);
            s = [];
            for j = 1:obj.Points/obj.fs/code_width
                if Code(2*j-1) == 0 && Code(2*j) == 0
                    if SignalType == 1
                        s = [s obj.V*cos(2*pi*f*t)];
                    else
                        s = [s HilbertTransfer(obj.V*cos(2*pi*f*t))/sqrt(2)];
                    end
                elseif Code(2*j-1) == 0 && Code(2*j) == 1
                    if SignalType == 1
                        s = [s obj.V*cos(2*pi*f*t+pi/2)];
                    else
                        s = [s HilbertTransfer(obj.V*cos(2*pi*f*t+pi/2))/sqrt(2)];
                    end
                elseif Code(2*j-1) == 1 && Code(2*j) == 0
                    if SignalType == 1
                        s = [s obj.V*cos(2*pi*f*t+pi)];
                    else
                        s = [s HilbertTransfer(obj.V*cos(2*pi*f*t+pi))/sqrt(2)];
                    end
                elseif Code(2*j-1) == 1 && Code(2*j) == 1
                    if SignalType == 1
                        s = [s obj.V*cos(2*pi*f*t+3*pi/2)];
                    else
                        s = [s HilbertTransfer(obj.V*cos(2*pi*f*t+3*pi/2))/sqrt(2)];
                    end
                end
            end
            s = obj.zero_fill_or_transaction(s);
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateBFSK(obj,f,Code,snr,SignalType)
            s = zeros(1,obj.Points);
            t = linspace(0,obj.TimeWidth/length(Code),obj.Points/length(Code));
            code_width = obj.Points/length(Code); % 每个码元的采样点数
            for j = 1:length(Code)
                if SignalType == 1
                    s((j-1)*code_width+1:j*code_width) = obj.V*cos(2*pi*(Code(j)+1)*f*t);
                else
                    s((j-1)*code_width+1:j*code_width) = HilbertTransfer(obj.V*cos(2*pi*(Code(j)+1)*f*t))/sqrt(2);
                end
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateSFM(obj,ff,fmin,fmax,snr,SignalType)
            t = linspace(0,obj.Points/obj.fs,obj.Points);
            if SignalType == 1
                s = obj.V*vco(sin(2*pi*ff*t),[fmin,fmax],obj.fs);
            else
                s = obj.V*vco(sin(2*pi*ff*t),[fmin,fmax],obj.fs);
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateCOSTAS(obj,COSTAS_sequence_length,f_min,f_max,snr,SignalType)
            if ~exist(['COSTAS n=' num2str(COSTAS_sequence_length) '.mat'],'file')
                getCostasArray(COSTAS_sequence_length);
            end
            CostasArrs = load(['COSTAS n=' num2str(COSTAS_sequence_length) '.mat']).Costas;
            index = randi(size(CostasArrs,1),1);
            arr = CostasArrs(index,:);
            fstep = (f_max - f_min)/(COSTAS_sequence_length - 1);
            f = f_min+fstep*(arr-1);
            t = linspace(0,(obj.Points/COSTAS_sequence_length-1)/obj.fs,obj.Points/COSTAS_sequence_length);
            s = zeros(1,obj.Points);
            code_width = obj.Points/COSTAS_sequence_length; % 每个码元采的点数
            for j=1:COSTAS_sequence_length
                if SignalType == 1
                    s((j-1)*code_width+1:j*code_width) = obj.V*cos(2*pi*f(j)*t);
                else
                    s((j-1)*code_width+1:j*code_width) = HilbertTransfer(obj.V*cos(2*pi*f(j)*t))/sqrt(2);
                end
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateTANFM(obj,snr,SignalType)
            a = 2.5;
            b = atan(a);
            T = 2*obj.TimeWidth;
            t = linspace(0,obj.Points/obj.fs,obj.Points);
            if SignalType == 1
                s = obj.V*sin(2*pi*obj.fs/(8+4*rand())*tan(2*b*t/T)/2/tan(b).*t);
            else
                s = obj.V*sin(2*pi*obj.fs/(4+4*rand())*tan(2*b*t/T)/2/tan(b).*t);
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateFrank(obj,N,f,n,time_width,snr,SignalType)
            FrankCode = zeros(N,N);
            for j = 1:N
                for k = 1:N
                    FrankCode(j,k) = (j-1)*(k-1);
                end
            end
            s = [];
            for j = 1:N
                for k = 1:N
                    t = linspace(0,time_width/N/N,fix(n/N/N));
                    if SignalType == 1
                        s = [s obj.V*cos(2*pi*f*t+2*pi*FrankCode(j,k)/N)];
                    else
                        s = [s HilbertTransfer(obj.V*cos(2*pi*f*t+2*pi*FrankCode(j,k)/N))/sqrt(2)];
                    end
                end
            end
            s = obj.zero_fill_or_transaction(s);
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateP1(obj,snr,N,SignalType)
            P1Code = zeros(N,N);
            for j = 1:N
                for k = 1:N
                    P1Code(j,k) = -pi/N*(N+1-2*k)*(N*(k-1)+(j-1));
                end
            end
            if SignalType == 2
                f = (0.1+0.1*rand(1))*obj.fs; % 采样点数较少时频率不应过高
            else
                f = (0.05+0.1*rand(1))*obj.fs;
            end
            s = [];
            for j = 1:N
                for k = 1:N
                    t = linspace((N*(j-1)+k-1)*obj.TimeWidth/N/N,(N*(j-1)+k-1)*obj.TimeWidth/N/N+obj.TimeWidth/N/N,fix(obj.Points/N/N));
                    if SignalType == 1
                        s = [s cos(2*pi*f*t+P1Code(j,k))];
                    else
                        s = [s HilbertTransfer(cos(2*pi*f*t+P1Code(j,k)))/sqrt(2)];
                    end
                end
            end
            s = obj.zero_fill_or_transaction(s);
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateP2(obj,snr,SignalType)
            % N必须是偶数
            N = floor(6+3*rand());
            while mod(N,2)~=0
                N = floor(6+3*rand());
            end
            P2Code = zeros(N,N);
            for j = 1:N
                for k = 1:N
                    P2Code(j,k) = pi/N*(N+1-2*k)*(N+1-2*j);
                end
            end
            if SignalType == 2
                f = (0.05+0.1*rand(1))*obj.fs;
            else
                f = (0.05+0.1*rand(1))*obj.fs;
            end
            s = [];
            for j = 1:N
                for k = 1:N
                    t = linspace((N*(j-1)+k-1)*obj.TimeWidth/N/N,(N*(j-1)+k-1)*obj.TimeWidth/N/N+obj.TimeWidth/N/N,fix(obj.Points/N/N));
                    if SignalType == 1
                        s = [s cos(2*pi*f*t+P2Code(j,k))];
                    else
                        s = [s HilbertTransfer(cos(2*pi*f*t+P2Code(j,k)))/sqrt(2)];
                    end
                end
            end
            s = obj.zero_fill_or_transaction(s);
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateP3(obj,snr,SignalType)
            if SignalType == 2
                f = (0.1+0.1*rand(1))*obj.fs;
                N=64;
            else
                f = (0.05+0.1*rand(1))*obj.fs;
                N=128;
            end
            for k = 1:N
                phase(k)=(pi/N)*(k-1)*(k-1); %Compute the P3 Phase
            end
            s = [];
            for j = 1:N
                t = linspace((j-1)*obj.TimeWidth/N,j*obj.TimeWidth/N,fix(obj.Points/N));
                if SignalType == 1
                    s = [s cos(2*pi*f*t+phase(j))];
                else
                    s = [s HilbertTransfer(cos(2*pi*f*t+phase(j)))/sqrt(2)];
                end
            end
            s = obj.zero_fill_or_transaction(s);
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateP4(obj,snr,SignalType)
            N=64;
            for k = 1:N
                phase(k)=((pi/N)*(k-1)^2)-(pi*(k-1)); % P4 Phase
            end
            if SignalType == 2
                f = (0.1+0.1*rand(1))*obj.fs;
            else
                f = (0.05+0.1*rand(1))*obj.fs;
            end
            s = [];
            for j = 1:N
                t = linspace((j-1)*obj.TimeWidth/N,j*obj.TimeWidth/N,fix(obj.Points/N));
                if SignalType == 1
                    s = [s cos(2*pi*f*t+phase(j))];
                else
                    s = [s HilbertTransfer(cos(2*pi*f*t+phase(j)))/sqrt(2)];
                end
            end
            s = obj.zero_fill_or_transaction(s);
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateT1(obj,snr,f,k,n,SignalType)
            time_width = n / obj.fs;
            t = linspace(0,time_width,n);
            j = fix(k*t/time_width);
            phase_state_num = 2;
            phi = mod(2*pi/phase_state_num*fix((k*t-j*time_width).*j*phase_state_num/time_width),2*pi);
            s = [];
            j = 1;
            for i = t
                s = [s cos(2*pi*f*i+phi(j))];
                j = j + 1;
            end
            if SignalType == 2
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateT2(obj,snr,f,k,n,SignalType)
            time_width = n / obj.fs;
            t = linspace(0,time_width,n);
            j = fix(k*t/time_width);
            phase_state_num = 2;
            phi = mod(2*pi/phase_state_num*floor((k*t-j*time_width).*(2*j-k+1)/time_width*phase_state_num/2),2*pi);
            s = [];
            j = 1;
            for i = t
                s = [s cos(2*pi*f*i+phi(j))];
                j = j + 1;
            end
            if SignalType == 2
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateT3(obj,snr,f,n,SignalType)
            time_width = n / obj.fs;
            deltaF = 0.05*obj.fs;
            phase_state_num = 2;
            t = linspace(0,time_width,n);
            phi = mod(2*pi/phase_state_num*floor(phase_state_num*deltaF*t.^2/2/time_width),2*pi);
            s = [];
            j = 1;
            for i = t
                s = [s cos(2*pi*f*i+phi(j))];
                j = j + 1;
            end
            if SignalType == 2
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateT4(obj,snr,f,n,SignalType)
            time_width = n / obj.fs;
            t = linspace(0,time_width,n);
            deltaF = 0.1*obj.fs;
            phase_state_num = 2;
            phi = mod(2*pi/phase_state_num*floor(phase_state_num*deltaF*t.^2/2/time_width-phase_state_num*deltaF*t/2),2*pi);
            s = [];
            j = 1;
            for i = t
                s = [s cos(2*pi*f*i+phi(j))];
                j = j + 1;
            end
            if SignalType == 2
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateLFM_SFM(obj,snr,SignalType)
            syms k t F B A C f0
            Dphi = 2*pi*f0+2*pi*(k*t+C)*(A+B*cos(2*pi*F*t));
            phi = int(Dphi,t);
            phi = matlabFunction(phi);
            t = linspace(obj.TimeWidth,2*obj.TimeWidth,obj.Points);
            TT = 5+rand()*5; % 正弦频率分量在TimeWidth的时间长度内震荡5到10个周期
            ff = TT/obj.TimeWidth;
            freqWidth = (0.1+0.3*rand(1))*obj.fs; % f ranges from 0.1*fs to 0.4*fs
            A = 3;
            B = 0.5;
            C = 30;
            f0 = 0.5*obj.fs-A*C; % initial frequency
            k = -1*freqWidth/A/obj.TimeWidth/2;
            if SignalType == 1
                s = obj.V*cos(phi(A,B,C,ff,f0,k,t));
            else
                s = obj.V*cos(phi(A,B,C,ff,f0,k,t));
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateEQFM_SFM(obj,snr,SignalType)
            syms f0 k t T A B F
            Dphi = 2*pi*(f0+3/2*k*(t-T)^2)*(A+B*cos(2*pi*F*t));
            phi = int(Dphi,t);
            phi = matlabFunction(phi);
            t = linspace(0,obj.TimeWidth,obj.Points);
            T = obj.TimeWidth/2;
            freqWidth = (0.1+0.3*rand(1))*obj.fs; % f ranges from 0.1*fs to 0.4*fs
            k = 2*freqWidth/3/T/T;
            A = 3;
            B = 0.5;
            f0 = rand(1)*0.1*obj.fs/A/2;
            TT = 6+rand()*6; % 正弦频率分量在TimeWidth的时间长度内震荡6到12个周期
            ff = TT/obj.TimeWidth;
            if SignalType == 1
                s = obj.V*cos(phi(A,B,ff,T,f0,k,t));
            else
                s = obj.V*cos(phi(A,B,ff,T,f0,k,t));
                s = HilbertTransfer(s)/sqrt(2);
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateDDC_MASK(obj,snr,E,M,symbol_rate,SignalType) % 生成下变频后的MASK信号
            assert(SignalType==2,'下变频后应当是复指数信号');
            symbol_period = 1 / symbol_rate;                    % 符号周期
            symbol_nums = fix(obj.Points/obj.fs/symbol_period); % 符号数量
            t = linspace(0,obj.Points/symbol_nums/obj.fs,obj.Points/symbol_nums);
            d = sqrt(3*E/(M^2-1));
            s = [];
            for k = 1:symbol_nums
                rand_sequence = randperm(M);
                code = rand_sequence(1);
                ak = (2*code-1-M)*d;
                if SignalType == 1
                    s = [s sqrt(E)*ak*ones(1,length(t))];
                elseif SignalType == 2
                    s = [s sqrt(E)*ak*ones(1,length(t))];
                end
            end
            s = obj.zero_fill_or_transaction(s);
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateDDC_MPSK(obj,snr,E,M,symbol_rate,SignalType) % 生成下变频后的MASK信号
            assert(SignalType==2,'下变频后应当是复指数信号');
            symbol_period = 1 / symbol_rate;                    % 符号周期
            symbol_nums = fix(obj.Points/obj.fs/symbol_period); % 符号数量
            t = linspace(0,obj.Points/symbol_nums/obj.fs,obj.Points/symbol_nums);
            d = sqrt(3*E/(M^2-1));
            s = [];
            for k = 1:symbol_nums
                rand_sequence = randperm(M);
                m = rand_sequence(1);
                phik = 2*pi*(m-1)/M;
                if SignalType == 1
                    s = [s sqrt(E)*cos(phik)*ones(1,length(t))];
                elseif SignalType == 2
                    s = [s sqrt(E)*exp(1i*phik)*ones(1,length(t))];
                end
            end
            s = obj.zero_fill_or_transaction(s);
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s,noise] = generateDDC_MFSK(obj,snr,E,M,deltaf,symbol_rate,SignalType) % 生成下变频后的MASK信号
            assert(SignalType==2,'下变频后应当是复指数信号');
            symbol_period = 1 / symbol_rate;                    % 符号周期
            symbol_nums = fix(obj.Points/obj.fs/symbol_period); % 符号数量
            t = linspace(0,obj.Points/symbol_nums/obj.fs,obj.Points/symbol_nums);
            d = sqrt(3*E/(M^2-1));
            s = [];
            for k = 1:symbol_nums
                rand_sequence = randperm(M);
                m = rand_sequence(1);
                omegak = (2*m-1-M)*deltaf*2*pi;
                if SignalType == 1
                    s = [s sqrt(E)*cos(omegak*t).*ones(1,length(t))];
                elseif SignalType == 2
                    s = [s sqrt(E)*exp(1i*omegak*t).*ones(1,length(t))];
                end
            end
            s = obj.zero_fill_or_transaction(s);
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [s, noise] = generateMQAM(obj,snr,SignalType,M,code,fc)
            if mod(length(code), length(dec2bin(M-1))) ~= 0
                error([num2str(M) 'QAM需要' num2str(length(dec2bin(M-1))) '个二进制数来表示一个符号，而出入的码元数量不能被' num2str(length(dec2bin(M-1))) '整除']);
            end
            constellation = qammod(0:M-1,M);
            symbol_num = length(code)/length(dec2bin(M-1)); % 符号数量
            symbol_width = obj.Points / obj.fs / symbol_num; % 每个符号的时域宽度
            s = [];
            for i = 1:symbol_num
                t = linspace(0, symbol_width, obj.Points / symbol_num);
                symbol = bin2dec(num2str(code((i-1)*symbol_num+1:i*symbol_num)));
                s_temp = real(exp(1i*2*pi*fc*t)*constellation(symbol)); 
                if SignalType == 2
                    s_temp = HilbertTransfer(s_temp);
                end
                s = [s s_temp];
            end
            if isinf(snr)
                noise = zeros(1,length(s));
                return
            end
            noise = obj.getNoise(s,snr,SignalType);
            s = s+noise;
        end
        function [P] = getP(obj, s) % 获取信号s的平均功率
            E = 0; % 总能量
            for i = 1:length(s)
                E = E + abs(s(i))^2/obj.fs;
            end
            P = E/(length(s)/obj.fs);
        end
        function [snr]     = getSNR(obj,s,noise) % 获取s+noise的信噪比
            Ps = obj.getP(s);
            Pnoise = obj.getP(noise);
            snr = 10*log(Ps/Pnoise)/log(10);
        end
        function []        = checkSNR(obj,s,noise,snr)
            if abs(obj.getSNR(s,noise)-snr)>0.5
                warning("信噪比偏差过大，请重新检查代码");
            end
        end
        function [res]     = zero_fill_or_transaction(obj,s) % 信号s的长度大于采样点数时截断，小于采样点数时补0
            if length(s) > obj.Points
                res = s(1:obj.Points);
                return;
            else
                res = [s zeros(1,obj.Points-length(s))];
            end
        end
        function noise = getNoise(obj,s,snr,SignalType)
            p = obj.getP(s);
            Pnoise = p/10^(snr/10);
            dBW = 10*log(Pnoise)/log(10);
            noise = wgn(1,length(s),dBW);
            if SignalType == 2
                noise = HilbertTransfer(noise)/sqrt(2);
            end
            obj.checkSNR(s,noise,snr);
        end
    end
end