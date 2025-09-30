clc;
clear;
close all;
addpath(genpath(pwd));
rng(1993); % For repeatable results

%%%%%*** Waveform Configuration ***%%%%%
% Create a format configuration object for a 1-by-1 HT transmission
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth = 'CBW20'; % 20 MHz channel bandwidth
cfgHT.NumTransmitAntennas = 1; % 1 transmit antennas
cfgHT.NumSpaceTimeStreams = 1; % 1 space-time streams
cfgHT.PSDULength = 2000; % PSDU length in bytes % 64
cfgHT.MCS = 0; % 1 spatial streams, BPSK rate-1/2
cfgHT.ChannelCoding = 'BCC'; % BCC channel coding

fs = wlanSampleRate(cfgHT); % Get the baseband sampling rate
ofdmInfo = wlanHTOFDMInfo('HT-Data',cfgHT); % Get the OFDM info
ind = wlanFieldIndices(cfgHT); % Indices for accessing each field within the time-domain packet

%%%%%*** Simulation Parameters ***%%%%%
snr = 0:5:15; %10:2:30;
global numTags;
numTags = 3;
global seqLenForEstChannel;
seqLenForEstChannel = 40;
preambleForEstChannel = survey_Mecha_funcGeneratePreamble(seqLenForEstChannel,numTags);

maxNumPackets = 2000; % The maximum number of packets at an SNR point

S = numel(snr); % 返回数组snr中元素的个数
numBitErrs = zeros(S,numTags);
berEst = zeros(S,numTags);

numDataSubcarrier = 52;
numPilotSubcarrier = 4;
estH = ones(numDataSubcarrier+numPilotSubcarrier,numTags);


for i = 1:S
    disp(['SNR: ',num2str(snr(i)),' dB...']);
    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = i;
    RandStream.setGlobalStream(stream);
    
    % Loop to simulate multiple packets
    n = 1; % Index of packet transmitted
    while  n<=maxNumPackets
        disp(['snr: ',num2str(snr(i)),' dB -> ','n: ',num2str(n),'-th packet']);
        %%%%%*** TX side ***%%%%%
        % Generate a packet waveform
        txPSDU = randi([0 1],cfgHT.PSDULength*8,1); % PSDULength in bytes
        tx = wlanWaveformGenerator(txPSDU,cfgHT);
        tx = [tx; zeros(15,cfgHT.NumTransmitAntennas)]; % Add trailing zeros to allow for channel filter delay
        
        exSig = [];
        H_TX_Tags = [];
        %%%%%*** TX-Tags backscatter channel & AWGN
        for idx_0 = 1:numTags
            bxCoeffForTxTag_real = -1+(1+1)*rand(1,1);
            bxCoeffForTxTag_imag = -1+(1+1)*rand(1,1);
            bxCoeffForTxTag_real = bxCoeffForTxTag_real*0.1;
            bxCoeffForTxTag_imag = bxCoeffForTxTag_imag*0.1;
            bxCoeffForTxTag = bxCoeffForTxTag_real+1i*bxCoeffForTxTag_imag;
            tmp_exSig = tx.*bxCoeffForTxTag;
            exSig = [exSig,tmp_exSig];
            H_TX_Tags = [H_TX_Tags,bxCoeffForTxTag];
        end
        
        %%%%%*** Tags side ***%%%%%
        % Backscatter at the tag
        temp = ceil((cfgHT.PSDULength*8+16+6)/26);
        numSymForPsdu = 0;
        numSymForTailPad = 0;
        if mod(temp,2) == 1
            numSymForPsdu = (numel(tx)-720-15-80-80-80)/80;
            numSymForTailPad = 2;
        else
            numSymForPsdu = (numel(tx)-720-15-80-80)/80;
            numSymForTailPad = 1;
        end
        numTagData = numSymForPsdu; % modulate one tag data per one symbol
        
        % Initial tags data
        tagData = zeros(numTagData,numTags);
        numPayload = numTagData-seqLenForEstChannel*numTags;
        actualPayloadBits = zeros(numPayload,numTags);
        for idx_1 = 1:numTags
            payload = randi([0,1],numPayload,1);
            actualPayloadBits(:,idx_1) = payload;
            tagData(:,idx_1) = [preambleForEstChannel(:,idx_1);payload];
        end
        
        % backscatter operation
        for idx_2 = 1:numTags
            bxSig{idx_2} = survey_Mecha_funcBackscatter(exSig(:,idx_2),tagData(:,idx_2),1);
        end
        
        %%%%%***** Backscatter channel & AWGN ***%%%%%
        H_Tags_RX = [];
        for idx_3 = 1:numTags
            bxCoeffForTagRx_real = -1+(1+1)*rand(1,1);
            bxCoeffForTagRx_imag = -1+(1+1)*rand(1,1);
            bxCoeffForTagRx_real = bxCoeffForTagRx_real*0.01;
            bxCoeffForTagRx_imag = bxCoeffForTagRx_imag*0.01;
            bxCoeffForTagRx = bxCoeffForTagRx_real+1i*bxCoeffForTagRx_imag;
            bxSig{idx_3} = bxSig{idx_3}.*bxCoeffForTagRx; % backscatter channel
            H_Tags_RX = [H_Tags_RX,bxCoeffForTagRx];
        end
        actualH = H_TX_Tags.*H_Tags_RX;
        
        
           
        %%%%%*** RX side ***%%%%%
        rx = complex(zeros(length(bxSig{1}),1));
        for idx_4 = 1:numTags
            rx = rx + bxSig{idx_4};
        end
%         rxFromTags = rx;
        [rxFromTags,~,~] = func_awgn(rx,snr(i),'measured'); % Received signal from Tags-RX channel
        
        
        
        [~,ofdmSymDerived] = survey_Mecha_funcOFDMSymDerived(txPSDU,cfgHT);
        [cfgOFDM,dataInd,pilotInd] = wlan.internal.wlanGetOFDMConfig(cfgHT.ChannelBandwidth, cfgHT.GuardInterval, 'HT', cfgHT.NumSpaceTimeStreams);
        ofdmDataDerived = ofdmSymDerived(cfgOFDM.DataIndices,:,:);
        ofdmPilotsDerived = ofdmSymDerived(cfgOFDM.PilotIndices,:,:);
        [ofdmDemodData,ofdmDemodPilots] = survey_Mecha_funcReceiver(rxFromTags(ind.HTData(1):ind.HTData(2)),cfgHT,1);
        ofdmDerived = [ofdmDataDerived;ofdmPilotsDerived];
        ofdmDemod = [ofdmDemodData;ofdmDemodPilots];
        
        % LS estimator
        for idx_5 = 1:numTags
            A = ofdmDemod(:,1+(1+(idx_5-1)*seqLenForEstChannel:idx_5*seqLenForEstChannel));
            B = ofdmDerived(:,1+(1+(idx_5-1)*seqLenForEstChannel:idx_5*seqLenForEstChannel));
            tmp_LL = size(A,1);
            for idx_LL = 1:tmp_LL
                tmp_H_real = funcLSEstimator(B(idx_LL,:)',real(A(idx_LL,:))');
                tmp_H_imag = funcLSEstimator(B(idx_LL,:)',imag(A(idx_LL,:))');
                tmp_H = tmp_H_real + 1i*tmp_H_imag;
                estH(idx_LL,idx_5) = tmp_H;
            end
        end
        
        payload_ofdmDemod = ofdmDemod(:,1+seqLenForEstChannel*numTags+1:end-numSymForTailPad); 
        payload_ofdmDerived = ofdmDerived(:,1+seqLenForEstChannel*numTags+1:end-numSymForTailPad);
        
        demodPayloadBits = survey_Mecha_funcDemd(estH,payload_ofdmDerived,payload_ofdmDemod);
        
        % calculate the number of bits
        for tt1 = 1:numTags
            numBitErrs(i,tt1) = numBitErrs(i,tt1) + biterr(actualPayloadBits(:,tt1),demodPayloadBits(:,tt1));
        end
        n = n+1;
        
    end
    % calculate bit error rate
    for tt2 = 1:numTags
        berEst(i,tt2) = numBitErrs(i,tt2)/(numPayload*maxNumPackets);
    end
    
end

aaa = 1;


