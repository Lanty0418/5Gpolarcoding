%% 5G New Radio Polar Coding
%
% This example highlights the new polar channel coding technique chosen for
% 5G New Radio (NR) communications system. Of the two main types of code
% constructions specified by 3GPP, this example models the CRC-Aided Polar
% (CA-Polar) coding scheme. This example describes the main components of
% the polar coding scheme with individual components for code construction,
% encoding and decoding along-with rate-matching. It models a polar-coded
% QPSK-modulated link over AWGN and presents Block-Error-Rate results for
% different message lengths and code rates for the coding scheme.

% Copyright 2018-2021 The MathWorks, Inc.


s = rng(100);       % Seed the RNG for repeatability

%%
% Specify the code parameters used for a simulation. 

% Code parameters
K = 54;             % Message length in bits, including CRC, K > 30
E = 124;            % Rate matched output length, E <= 8192

%EbNo = 0.8;         % EbNo in dB
L = 8;              % List length, a power of two, [1 2 4 8]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numFrames = 10000;     % Number of frames to simulate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linkDir = 'DL';     % Link direction: downlink ('DL') OR uplink ('UL')


if strcmpi(linkDir,'DL')
    % Downlink scenario (K >= 36, including CRC bits)
    crcLen = 24;      % Number of CRC bits for DL, Section 5.1, [6]
    poly = '24C';     % CRC polynomial
    nPC = 0;          % Number of parity check bits, Section 5.3.1.2, [6]
    nMax = 9;         % Maximum value of n, for 2^n, Section 7.3.3, [6]
    iIL = true;       % Interleave input, Section 5.3.1.1, [6]
    iBIL = false;     % Interleave coded bits, Section 5.4.1.3, [6]
else
    % Uplink scenario (K > 30, including CRC bits)
    crcLen = 11;      
    poly = '11';
    nPC = 0;          
    nMax = 10;        
    iIL = false;      
    iBIL = true;      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_vec=[56 43 54 164 164];
E_vec=[864 180 124 240 184];
snrVec=-12:1:8;  % or snrVec=-12:2:8;  is ok
blockErrRateAllLayers = zeros(length(K_vec), length(snrVec));
bitErrRateAllLayers = zeros(length(K_vec), length(snrVec));
legendLabels = cell(1, length(K_vec));  % 新增一個儲存圖例標籤的變數
  
for idx = 1:length(K_vec)
    K = K_vec(idx);
    E = E_vec(idx);
    R=K_vec(idx)/E_vec(idx);

    legendLabels{idx} = ['K/E = ' num2str(K/E)];  % 設定圖例標籤

for snrIdx=1:length(snrVec)
snrdB=snrVec(snrIdx);  
noiseVar = 1./(10.^(snrdB/10)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Channel
chan = comm.AWGNChannel('NoiseMethod','Variance','Variance',noiseVar);


% Error meter
ber = comm.ErrorRate;

%% Frame Processing Loop
numferr = 0;
for i = 1:numFrames

    % Generate a random message
    msg = randi([0 1],K-crcLen,1);
    
    % Attach CRC
    msgcrc = nrCRCEncode(msg,poly);
    
    % Polar encode
    encOut = nrPolarEncode(msgcrc,E,nMax,iIL);
    N = length(encOut);
    
    % Rate match
    modIn = nrRateMatchPolar(encOut,K,E,iBIL);
    
    % Modulate
    modOut = nrSymbolModulate(modIn,'QPSK');
    
    % Add White Gaussian noise
    rSig = chan(modOut);
    
    % Soft demodulate
    rxLLR = nrSymbolDemodulate(rSig,'QPSK',noiseVar);
    
    % Rate recover
    decIn = nrRateRecoverPolar(rxLLR,K,N,iBIL);
    
    % Polar decode
    decBits = nrPolarDecode(decIn,K,E,L,nMax,iIL,crcLen); 
    
    % Compare msg and decoded bits
    errStats = ber(double(decBits(1:K-crcLen)), msg);
    numferr = numferr + any(decBits(1:K-crcLen)~=msg);

end
blockErrRateAllLayers(idx, snrIdx) = numferr / numFrames;
        bitErrRateAllLayers(idx, snrIdx) = errStats(1);
disp(['Block Error Rate: ' num2str(numferr/numFrames) ...
      ', Bit Error Rate: ' num2str(errStats(1)) ...
      ', at SNR = ' num2str(snrdB) ' dB'])
end
end
% 用圖呈現 Block Error Rate
figure;
for idx = 1:length(K_vec)
    plot(snrVec, blockErrRateAllLayers(idx, :), '-o', 'DisplayName', legendLabels{idx});  % 使用新的圖例標籤
    hold on;
end
title('Block Error Rate vs SNR ');
xlabel('SNR (dB)');
ylabel('Block Error Rate');
legend('show');
set(gca, 'YScale', 'log'); % 將 y 軸設為對數尺度
ylim([1e-4, 1]); % 設定 y 軸範圍

% 用圖呈現 Bit Error Rate
figure;
for idx = 1:length(K_vec)
    plot(snrVec, bitErrRateAllLayers(idx, :), '-o', 'DisplayName', legendLabels{idx});  % 使用新的圖例標籤
    hold on;
end
title('Bit Error Rate vs SNR ');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
legend('show');
set(gca, 'YScale', 'log'); % 將 y 軸設為對數尺度
ylim([1e-5, 1]); % 設定 y 軸範圍

rng(s); % 恢復 RNG