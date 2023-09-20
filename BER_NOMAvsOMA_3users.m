%% NOMA vs OMA BER comparison for 3 user scenario using OFDM blocks
%% Clear the workspace, command window, and close all figures
clc;
clear;
close all;
%%
% Define user distances from the base station (customize accordingly)
userDistances = [80, 40, 20];  % User distances in meters
% Define simulation parameters
numSymbols = 32;               % Number of symbols
numSubcarriers = 128;          % Number of subcarriers
CPlength = 16;                 % Cyclic Prefix length
channelLength = 16;            % Channel length
SNRdBRange = -10:10:80;        % SNR values in dB
numBlocks = 1024;              % Number of blocks
meanSquare = 1;                % Mean square value
signalPower = 1;               % Power of BPSK signal
pathLossExponent = 2;          % Path loss exponent
% Channel modeling for all users
h = cell(1, 3);                % Cell array to store channel responses
% User data and modulated signal for all user
dataUser = cell(1, 3);         % Cell array to store user data (NOMA)
signalUser = cell(1, 3);       % Cell array to store user signal (NOMA)
%% NOMA Simulation
% Arrays to store BER results for NOMA
BERUser_NOMA = zeros(3, length(SNRdBRange));  % Rows: Users, Columns: SNR values
for SNRIdx = 1:length(SNRdBRange)
    SNRdB = SNRdBRange(SNRIdx);
    SNR = 10^(SNRdB/10);
    noiseVar = meanSquare * signalPower / SNR;
    BERPerBlock = zeros(3, numBlocks); % Array to store BER per block
    for i = 1:numBlocks
        % Transmitter
        for user = 1:3
            dataUser{user} = randi([0 1], 1, numSymbols);
            signalUser{user} = 2 * dataUser{user} - 1;  % BPSK modulation
            channelGain = sqrt(userDistances(user)^-pathLossExponent);
            h{user} = channelGain * (randn(1, channelLength) + sqrt(-1) * randn(1, channelLength)) / sqrt(2); % Rayleigh Channel
        end
        % Superposition coded signal
        signal = sqrt(0.7) * signalUser{1} + sqrt(0.2) * signalUser{2} + sqrt(0.1) * signalUser{3};
        timeDomainSignal = ifft(signal, numSymbols);  % Inverse Discrete Fourier Transform (IDFT)
        % Add cyclic prefix
        signalWithCP = [timeDomainSignal(end-CPlength+1:end), timeDomainSignal];
        for user = 1:3
            if user == 1
                n = sqrt(noiseVar/2) * (randn(1, numSymbols+CPlength) + sqrt(-1) * randn(1, numSymbols+CPlength)); % AWGN noise
                receivedSignal = cconv(signalWithCP, h{user}, numSymbols+CPlength) + n;
                receivedSignal = receivedSignal(CPlength+1:end);
                receivedSymbol = fft(receivedSignal);
                channelFreqResponse = fft(h{user}, numSymbols);
                estimatedSymbol = receivedSymbol ./ channelFreqResponse;
                estimatedBits = (estimatedSymbol > 0);
                BERPerBlock(user, i) = sum(dataUser{user} ~= estimatedBits) / numSymbols;
            elseif user == 2
                % User 2 Receiver (Decode the user 1 data)
                n = sqrt(noiseVar/2) * (randn(1, numSymbols+CPlength) + sqrt(-1) * randn(1, numSymbols+CPlength)); % AWGN noise
                receivedSignal = cconv(signalWithCP, h{user}, numSymbols+CPlength) + n;
                receivedSignal = receivedSignal(CPlength+1:end);
                receivedSymbol = fft(receivedSignal);
                channelFreqResponse = fft(h{user}, numSymbols);
                estimatedSymbol2 = receivedSymbol ./ channelFreqResponse;
                estimatedBits21 = (estimatedSymbol2 > 0);
                % User 2 performs SIC (Successive Interference Cancellation) on user 1
                estimatedSymbol2_1 = estimatedSymbol2 - sqrt(0.7) * (2 * estimatedBits21 - 1);
                estimatedBits = (estimatedSymbol2_1 > 0);
                BERPerBlock(user, i) = sum(dataUser{user} ~= estimatedBits) / numSymbols;
            else
                % User 3 Receiver (Decode the user 1 data)
                n = sqrt(noiseVar/2) * (randn(1, numSymbols+CPlength) + sqrt(-1) * randn(1, numSymbols+CPlength)); % AWGN noise
                receivedSignal = cconv(signalWithCP, h{user}, numSymbols+CPlength) + n;
                receivedSignal = receivedSignal(CPlength+1:end);
                receivedSymbol = fft(receivedSignal);
                channelFreqResponse = fft(h{user}, numSymbols);
                estimatedSymbol3 = receivedSymbol ./ channelFreqResponse;
                estimatedBits31 = (estimatedSymbol3 > 0);
                % User 3 Receiver (Decode the user 2 data by performing SIC on user 1)
                estimatedSymbol3_1 = estimatedSymbol3 - sqrt(0.7) * (2 * estimatedBits31 - 1);
                estimatedBits32 = (estimatedSymbol3_1 > 0);
                % User 3 performs SIC (Successive Interference Cancellation) on user 1 and user 2
                estimatedSymbol3_1_2 = estimatedSymbol3 - sqrt(0.7) * (2 * estimatedBits31 - 1) - sqrt(0.2) * (2 * estimatedBits32 - 1);
                estimatedBits = (estimatedSymbol3_1_2 > 0);
                BERPerBlock(user, i) = sum(dataUser{user} ~= estimatedBits) / numSymbols;
            end
        end
    end
    % Calculate BER for the current User and store the result
    for user = 1:3
        numErrors = sum(BERPerBlock(user, :));
        calculatedBER = numErrors / (numBlocks * numSymbols);
        BERUser_NOMA(user, SNRIdx) = calculatedBER;
    end
end
%% OMA Simulation
% Arrays to store BER results for OMA
BERUser_OMA = zeros(3, length(SNRdBRange));  % Rows: Users, Columns: SNR values
for user = 1:3
    for SNRIdx = 1:length(SNRdBRange)
        SNRdB = SNRdBRange(SNRIdx);
        SNR = 10^(SNRdB/10);
        noiseVar = meanSquare * signalPower / SNR;
        BERPerBlock = zeros(1, numBlocks); % Array to store BER per block
        for i = 1:numBlocks
            % Transmitter for current User in OMA
            dataUser = randi([0 1], 1, numSymbols);
            signalUser = 2 * dataUser - 1;  % BPSK modulation
            % OMA: Separate signals for current User
            signalUserTimeDomain = ifft(signalUser, numSymbols);  % Inverse Discrete Fourier Transform (IDFT)
            % Add cyclic prefix for current User
            signalWithCPUser = [signalUserTimeDomain(end-CPlength+1:end), signalUserTimeDomain];
            % Channel modeling (with distance-dependent channel gain) for current User in OMA
            channelGain = sqrt(userDistances(user)^-pathLossExponent);
            h = channelGain * (randn(1, channelLength) + sqrt(-1) * randn(1, channelLength)) / sqrt(2); % Rayleigh Channel
            n = sqrt(noiseVar/2) * (randn(1, numSymbols+CPlength) + sqrt(-1) * randn(1, numSymbols+CPlength)); % AWGN noise
            % Receiver 
            receivedSignal = cconv(signalWithCPUser, h, numSymbols+CPlength) + n;
            receivedSignal = receivedSignal(CPlength+1:end);
            receivedSymbol = fft(receivedSignal);
            channelFreqResponse = fft(h, numSymbols);
            estimatedSymbol = receivedSymbol ./ channelFreqResponse;
            estimatedBits = (estimatedSymbol > 0);
            BERPerBlock(i) = sum(dataUser ~= estimatedBits) / numSymbols;
        end
        % Calculate BER for the current User in OMA and store the result
        numErrors = sum(BERPerBlock);
        calculatedBER = numErrors / (numBlocks * numSymbols);
        BERUser_OMA(user, SNRIdx) = calculatedBER;
    end
end
%% Plot BER vs SNR
figure
for user = 1:3
    semilogy(SNRdBRange, BERUser_NOMA(user, :), '*-', 'LineWidth', 1.5); hold on;
    semilogy(SNRdBRange, BERUser_OMA(user, :), 'o--', 'LineWidth', 1.5); hold on;
end
title('NOMA vs OMA BER', 'Interpreter', 'Latex');
xlabel('SNR (dB)', 'Interpreter', 'Latex');
ylabel('BER', 'Interpreter', 'Latex');
legend('NOMA User 1 (weaker)', 'OMA User 1 (weaker)', 'NOMA User 2 (mid)', 'OMA User 2 (mid)', ...
        'NOMA User 3 (stronger)', 'OMA User 3 (stronger)', 'Interpreter', 'Latex');
ylim([10^-5 10^-1]);
grid on;
