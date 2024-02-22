%Course Project Part 1 task 2
%Author: Katelyn, Faith, Haoze

clear all
close all
% Parameters
M=2048;   %# of subcarriers
N=4;    %Number of OFDM-QPSK SYMBOLS symbols
L=200;    %taps
B=8e3;
fc=24e3;

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
d_original = zeros(M*N,1);

%result symbol
d_result =  zeros(M*N, 1);

%Passband symbols
xPB_syms= [];

%Creates Chirp
   t1=0.05;
   t2=0.2;
   f_0=fc-4000;
   B_ch=B/t1;
   t=[0:(t1*Fs)-1]/Fs;
   f_t=f_0+B_ch*t;
   ch=cos(2*pi*f_t.*t).';
   ch_zp=zeros(t2*Fs, 1);
   chirp=[ch; ch_zp; ch; ch_zp; ch; ch_zp; ch; ch_zp];
 

%**************iteration start*********************
for k=0:N-1
%Generate 1 OFDM SYMBOL
  for j = 1 : M
     randomIndex = randi(length(QPSKmap), 1);
     s(j) = QPSKmap(randomIndex);%random choose from QPSK
     d_original(j+k*M) = s(j);
  end
  
  %Inverse FFT
  x=ifft(s);
  
  %Add zp
  tem = zeros(L-1,1);
  s_cp=[tem; x;tem];

  %upsample
  lambda=Fs/Fd;
  xup=upsample(s_cp, lambda);
  t=[0:length(xup)-1]/Fs;
  figure(), plot(t, abs(xup))

  %Square Root cosine Pulse
   sr_cos=rcosdesign(beta, 2*delay, lambda, 'sqrt');
   xBB=conv(xup, sr_cos);
   
   %RCos Plot
   t=[0:length(sr_cos)-1]/Fs;         %Time Domain Plot
   figure(), plot(t, abs(sr_cos))
    
   SR_cos=fftshift(fft(sr_cos));        %Freq Domain Plot
   f=[-length(SR_cos)/2:(length(SR_cos))/2-1]*Fs;
   figure(), plot(f, abs(SR_cos))

   %Data Plot
   t=[0:length(xBB)-1]/Fs;                       %Time Domain Plot
   figure(), plot(t, abs(xBB))
    
   XBB=fftshift(fft(xBB));                     %Freq Domain Plot
   f=[-length(xBB)/2:(length(xBB))/2-1]/length(xBB)*Fs;
   figure(), plot(f, abs(XBB))
  
   %Passband
   Ts=1/B;
   ts=Ts/lambda;
   n=[0:length(xBB)-1].';
   xPB=real(exp(1j*2*pi*fc*n*ts).*xBB);
   len_pb=length(xPB);
   xPB_syms(1+k*len_pb:k*len_pb+len_pb)=xPB(1:len_pb);
    
   %Convert Signal back to baseband
   xBB_i=xPB.*sin(2*pi*fc*n*ts)*2;
   xBB_r=xPB.*cos(2*pi*fc*n*ts)*2;
   xBB=xBB_r+1i*xBB_i;

   %Match Filter
   xOUF=conv(xBB, sr_cos);

   %Down Sample
   xBBd=xOUF(delay*lambda+1:lambda:end-delay*lambda);

   %fft to convert symbols back into frequency domain
    for ii = 1:M
         sum = 0;
         for kk = 1: M+(L-1)
             sum = sum + xBBd(kk) .* exp(-(1j*2*pi*(kk-1)*(ii-1))/M);
         end
         y_freq(ii) = sum;
    end
   

   d_result(1+k*M:k*M+M) = y_freq(1:M);
end
    
   %Add Chirp to OFDM Symbol
   xPB_ch=[chirp; xPB_syms; chirp];
   
   %Plot of chirp spectrogram
   Nfft=4096;
   figure(), spectrogram(xPB_ch, Nfft, Nfft*3/4, Nfft, Fs, 'yaxis' )
   colormap jet
   set(gca, 'Fontsize', 12, 'FontWeight', 'bold')
   title('Transmission Spectrogram')
%***********iteration end***************************
%decide if bits are 1 or 0
