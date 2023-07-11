clc
close all
clearvars

%% Script parameters:
M =64;                                 % Modulation order 
m = log2(M);                             % Binary word size
N = m*10000;                             % Number of bits
EbNo = linspace(1,100,100);     % EbNo noise for AWGN channel
f = 1E9;
%% Simulate information source:
rng default;                             % Init seed of the random gen
txBit = randi([0 1],N,1);                % Generate array of random bits

%% Modulator:
txInteger = bit2int(txBit,m);            % Convert group of m bits
txSymbol = qammod(txInteger,M,'gray');   % Generate symbols

%% Channel:
SNR = EbNo + 10*log10(m);  
% Add to EbNo the signal 

for i = 1:length(SNR)
    



    rxSymbol(i,:) = awgn(txSymbol,SNR(i),'measured'); 
    rxInteger(:,i) = qamdemod(rxSymbol(i,:),M,'gray')';
    rxBit(:,i) = int2bit(rxInteger(:,i),m);
end
i = 0;
%% constaliation 


figure_1 = scatterplot(rxSymbol(0.20*length(EbNo),:),1,0,'g.');
hold on;grid on;
scatterplot(txSymbol,1,0,'r*',figure_1);


%% BER
error = zeros(length(rxBit(1,:)),1);
for i = 1:1:length(rxBit(1,:)) % For each SNR level
    for q = 1:length(rxBit(:,i))
        %disp(rxBit(q,i))
        if rxBit(q,i) ~= txBit(q)
           error(i) = error(i)+1;
        end
    end
end


figure_2 = figure;
semilogy(SNR,(error./length(txBit)));
grid on;
xlabel("SNR in DB")
ylabel("Bit Error Rate")
    
