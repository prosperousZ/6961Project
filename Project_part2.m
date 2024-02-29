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
Trx =
a_hat = Ttx/Trx -1;
YPB_re = resample(Ypb, round((1+a_hat) * 1e5),1e5);
Ls = 192;
Ms = 256;
Lp = 24;
N = Lp*Ls - 1;
h = Ls * fir1(N,1/Ms,kaiser(N+1 , 7.8562));
YPB_re_hat = upfirdn(YPB_re,h,Ls,Ms);

