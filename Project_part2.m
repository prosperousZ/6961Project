clear all

k = 2048;%the total numbver of subcaiiers
L = 200;%nmber of zero-padded symbols
numberOfTap = L+1;%number of taps
Fc = 24000;%carrier freq
B = 8000;%bandwidth
samplingRate = 256000;%sampleing-rate
W = 24;%number of ZP-OFDM symbols in one packet
Lambda = 24;%the oversampling factor
Y = load('benchmark_rece_data_174623_1472.mat');
Y = struct2array(Y);
%remove noise
Ypb = bandpass(Y,[-1000+Fc,8000 + Fc],samplingRate);

%underwater paras
c = 1500;
v = 1.03;
a = v/c;

%maxX = max(Ypb, [], "all");

%plot(Ypb);
Ttx = 8.2695;
Trx = (2117580 - 700) * (1/samplingRate);% sample
a_hat = Ttx/Trx -1;
YPB_re = resample(Ypb, round((1+a_hat) * 1e5),1e5);
Ls = 192;
Ms = 256;
Lp = 24;
N = Lp*Ls - 1;
h = Ls * fir1(N,1/Ms,kaiser(N+1 , 7.8562));
YPB_re_hat = upfirdn(YPB_re,h,Ls,Ms);
%*******This part is synchronization************
pilot= load('pilot_signal_for_synchronization.mat');
pilot = struct2array(pilot);
correlation = xcorr(YPB_re_hat,pilot);
figure(2)
plot(abs(correlation));

%starting point of OFDM symbol
n0 = 292726;
%n0_withGap = n0 + 2400;
YPB_re_hat = YPB_re_hat(n0:end);

vector= load('itc_1007_compfilter.mat');
vector = struct2array(vector);
figure (3)
plot(vector)
%convolve signal with the "vector"
YPB_re_hat_vector = conv(YPB_re_hat, vector);
delay = 50;
YPB_re_hat_vector = YPB_re_hat_vector(delay+1:end);

fs = 192000;
ts = 1/fs;
n_passbandToBaseBand = [0:length(YPB_re_hat_vector)-1].';
YBB_I = real(YPB_re_hat_vector .* 2.*cos(2.*pi.*Fc.*ts .*n_passbandToBaseBand ));
YBB_Q = -imag(YPB_re_hat_vector .* 2.*sin(2.*pi.*Fc.*ts .*n_passbandToBaseBand ));

YBB = YBB_I + 1j.*YBB_Q;
figure(4)
plot(YBB)

%***********************step 8, match filtering******************
beta = 0.125;
delay_R = 100;
%removed delay
R=rcosdesign(beta, 2*delay_R, Lambda, 'sqrt');
YBB_est = conv(YBB,R);
%This step is to remove delay after conv
YBB_est = YBB_est(2401:end);
%figure(5)
%plot(YBB_est);

%declare n, i, and z_m_1
n = [1:(k+L)*Lambda].';
i = [1:k+L].';
z_m_1 = zeros(k,1);

%power over null subcarriers
p = zeros((2400-2200)/1,(2-(-2))/0.1);
for n_hat_0_1 = [2200:1:2400]
    for epsilon_1 = [-2:0.1:2]
        %CFO compensation:
        YBB_hat_1 = YBB_est(n_hat_0_1 + n-1) .* exp((-1j .* 2.*pi.*epsilon_1).*(n-1+n_hat_0_1) .*ts);
        %Down-sampling
        YBB_hat_1 = YBB_hat_1(i*Lambda);
        for m = 1:k
            z_m_1(m) = sum(YBB_hat_1(i) .*exp(-1j * (2.*pi.*(m-1) .* (i-1))./k));
            %for q = 1:length(i)
             %   z_m_1(m) = z_m_1(m) + YBB_hat_1(q) * exp(-1j*((2*pi*(m-1)*(q-1))/k));
            %end
        end
        %for w = 1:k
         %   p(n_hat_0_1-2200+1,(epsilon_1-(-2))/0.1+1) =  sum()
        %end
        %p(n_hat_0_1-2200+1,(epsilon_1-(-2))/0.1+1) = sum(z_m_1(i))
    end
end



