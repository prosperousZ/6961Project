clear all

N = 4;
M = 16;
L = 5;

h = [0.227, 0.46, 0.688,0.46,0.227];
s = [1+1i, 1-1i, -1+1i,-1-1i]/sqrt(2);


symbol_cp = zeros(M+L-1,1);% create symbol array

for i = 5 : M+L-1
    randomIndex = randi(length(s ), 1);
     symbol_cp(i) = s(randomIndex);%random choose from QPSK
end

%add cp
cp = L-1;
for j = 1:cp
    symbol_cp(j) = symbol_cp(M+j);
end




