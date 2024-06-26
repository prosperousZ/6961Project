%Course Project Part 1
%Author: Katelyn, Faith, Haoze

clear all
close all
% Parameters
M=16;   %# of subcarriers
N=4;    %Number of OFDM-QPSK SYMBOLS symbols
L=5;    %taps
%Delcare channel h
h=[0.227, 0.46, 0.688, 0.46, 0.227];
%QPSK Symbol Map
QPSKmap= [1+1j, -1+1j, -1-1j, 1-1j];
s = zeros(M,1);

%Store 4 symbols here
d_original = zeros(M,4);

%result symbol
d_result =  zeros(M,4);

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

  cp=x(length(s)-(L-1)+1:end);
  s_cp=[cp; x];



  %Y Output after going through channel impulse response of data w/o cp
  y=zeros(M,1);
  for i=L:M+L-1
      y(i-(L-1))=h(1)*s_cp(i)+h(2)*s_cp(i-1)+h(3)*s_cp(i-2)+h(4)*s_cp(i-3)+h(5)*s_cp(i-4);
  end

  y_freq=fft(y);
  
  h_freq = fft(h, M).';
  d_freq=y_freq./h_freq;
  for w = 1 : M
    
     d_result(w,k) = d_freq(w);
  end

end
%***********iteration end***************************
%decide if bits are 1 or 0
