%Course Project Part 1 task 2
%Author: Katelyn, Faith, Haoze

clear all
close all
% Parameters
M=1024;   %# of subcarriers
N=4;    %Number of OFDM-QPSK SYMBOLS symbols
L=5;    %taps
%Delcare channel h
h=[0.227, 0.46, 0.688, 0.46, 0.227];
%h=[1 0 0 0 0];
%QPSK Symbol Map
QPSKmap= [1+1j, -1+1j, -1-1j, 1-1j];
s = zeros(M,1);
y_freq = zeros(M,1);
y=zeros(M+(L-1),1);

Fs=192e3;
Fd=8e3;
beta=0.125;
delay=100;

%Store 4 symbols here
d_original = zeros(M,N);

%result symbol
d_result =  zeros(M,N);

%**************iteration start*********************
for k=1:N
%Generate 1 OFDM SYMBOL
  for j = 1 : M
     randomIndex = randi(length(QPSKmap), 1);
     s(j) = QPSKmap(randomIndex);%random choose from QPSK
     d_original(j,k) = s(j);
  end
  
    %Inverse FFT
  x=ifft(s);
  
  tem = zeros(L-1,1);

  s_cp=[tem; x;tem];
  %Y Output after going through channel impulse response of data w/o cp

  %upsample
  lambda=Fs/Fd;
  xup=upsample(s_cp, lambda);

  %Square Root cosine Pulse
   sr_cos=rcosdesign(beta, 2*delay, lambda, 'sqrt');
   %xBB=conv(xup, sr_cos);
   figure(), impz(sr_cos)
  
   %Pass
 

  for i=L:M+2*(L-1)
      y(i-(L-1))=h(1)*s_cp(i)+h(2)*s_cp(i-1)+h(3)*s_cp(i-2)+h(4)*s_cp(i-3)+h(5)*s_cp(i-4);
  end
  
 
  %step 5 from lecture 11 instructor' note which is similar to fft(y), but
  %not the same
  for ii = 1:M
      sum = 0;
      for kk = 1: M+(L-1)
          sum = sum + y(kk) .* exp(-(1j*2*pi*(kk-1)*(ii-1))/M);
      end
      y_freq(ii) = sum;
  end
  
  %step 6
  h_freq = fft(h, M).';
  d_freq=y_freq./h_freq;
  for w = 1 : M
    
     d_result(w,k) = y_freq(w);
  end

end
%***********iteration end***************************
%decide if bits are 1 or 0
