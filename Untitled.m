% configuration d'un équipement utilisateur
% The channel sampling rate depends on the FFT size used in the OFDM
% modulator. This can be obtained using the function lteSLSCFDMAInfo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ueConfig = struct('SidelinkMode','V2X');% Release 14 V2X mode
ueConfig.NSLRB = 1;                    % 10MHz bandwidth
ueConfig.DuplexMode = 'FDD';            % Duplex mode
ueConfig.CyclicPrefixSL = 'Normal';     % Cyclic prefix length
ueConfig.Modulation = '16QAM'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDSCH configuration for transmission mode 10 (TM10)
enb = struct;
enb.NDLRB = 6;
enb.CellRefP = 2;
enb.DuplexMode = 'FDD';
enb.CFI = 3;
enb.CyclicPrefix = 'Normal';
enb.NFrame = 0;
enb.NCellID = 0;

enb.PDSCH.TxScheme = 'Port7-14';
enb.PDSCH.NLayers = 1;
enb.PDSCH.RNTI = 1;
enb.PDSCH.NSCID = 0;
enb.PDSCH.Modulation = {'16QAM'};
enb.PDSCH.Rho = 0;
enb.PDSCH.RV = 0;
enb.PDSCH.NTurboDecIts = 5;
enb.PDSCH.PRBSet = (0:enb.NDLRB-1).';
enb.PDSCH.NTxAnts = 4;
enb.PDSCH.W = lteCSICodebook(enb.PDSCH.NLayers,enb.PDSCH.NTxAnts,0).';
enb.PDSCH.AltCodebook4Tx = 'Off';
enb.PDSCH.CSI = 'On';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SimCSIPeriod = [5 1];
SimCSIRS.CSIRSConfig = [0 5];                       % CSI-RS configuration
SimCSIRS.CSIRSPeriod = SimCSIPeriod; % CSI-RS period
SimCSIRS.NCSIID = [10 16];                          % CSI-RS scrambling identity
SimCSIRS.CSIRefP = enb.PDSCH.NTxAnts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SimCSIProcess.CSIRSResource  = [1 2 1 2]; % CSI-RS resources index (1 based)
SimCSIProcess.CSIIMResource  = [1 1 2 3]; % CSI-IM resources index (1 based)
SimCSIProcess.CSIMode        = {'PUCCH 1-1','PUSCH 3-1','PUSCH 3-1','PUSCH 3-1'}; % CSI reporting modes
SimCSIProcess.PMIMode        = {'Wideband' ,'Wideband' ,'Wideband','Wideband'};   % PMI reporting modes
SimCSIProcess.CodebookSubset = {'0x0000000000000001','000001','0x0000000000000001','000001'}; % Codebook subset restrictions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fading Channel and SNR Configuration
snrTP = [10 8]; % SNR of received transmission from TP1 and TP2
Noc = -174;      % dBm/15kHz average power spectral density

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel between TP1 and UE
ofdmInfo = lteOFDMInfo(enb);
chcfg = struct;
chcfg.DelayProfile = 'EPA';
chcfg.NRxAnts = 2;
chcfg.DopplerFreq = 5;
chcfg.MIMOCorrelation = 'Low';
chcfg.SamplingRate = ofdmInfo.SamplingRate;
chcfg.InitPhase = 'Random';
chcfg.ModelType = 'GMEDS';
chcfg.NTerms = 16;
chcfg.NormalizeTxAnts = 'On';
chcfg.NormalizePathGains = 'On';
chcfg.Seed = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%configuration channel for vehicular users
chcfg(2) = struct;                    % Channel config structure
chcfg(2).Seed = 6;                    % Channel seed / 6 par defaut
chcfg(2).NRxAnts = 2;                 % Number of receive antennas / 2par défaut
chcfg(2).DelayProfile ='EVA';         % Delay profile
chcfg(2).DopplerFreq = 500;           % Doppler frequency in Hz/ 500 par défaut
chcfg(2).MIMOCorrelation = 'High';     % Multi-antenna correlation/ 0 par défaut
chcfg(2).NTerms = 16;                 % Oscillators used in fading model
chcfg(2).ModelType = 'GMEDS';         % Rayleigh fading model type
chcfg(2).InitPhase = 'Random';        % Random initial phases
chcfg(2).NormalizePathGains = 'On';   % Normalize delay profile power
chcfg(2).NormalizeTxAnts = 'On';      % Normalize for transmit antennas
%channel1.Seed = 1;
channel1.SamplingRate = ofdmInfo.SamplingRate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rxWaveformSize = [channel1.SamplingRate*1e-3+15 channel1.NRxAnts];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%configuration channel for pedestrian users
% chcfg = struct;
% chcfg.DelayProfile = 'EPA';
% chcfg.NRxAnts = 2;
% chcfg.DopplerFreq = 5;
% chcfg.MIMOCorrelation = 'Low';
% %chcfg.SamplingRate = ofdmInfo.SamplingRate;
% chcfg.InitPhase = 'Random';
% chcfg.ModelType = 'GMEDS';
% chcfg.NTerms = 16;
% chcfg.NormalizeTxAnts = 'On';
% chcfg.NormalizePathGains = 'On';
% chcfg.Seed = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify UE-specific parameters
muNumLayers = [1 1 1 1 1 1];  % Number of layers for a maximum of 6 users
muNumRxAnts = [1 1 1 1 1 1];  % Number of receive antennas for a maximum of 6 users
muCodeRate = [0.5 0.5 0.5 0.5 0.5 0.5]; % Code rate for a maximum of 6 users
muModulation = {'16QAM';'16QAM';'16QAM';'16QAM';'16QAM';'16QAM'}; % Modulation for a maximum of 6 users 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ceccsi.FreqWindow = 1;
ceccsi.TimeWindow = 2;
ceccsi.InterpType = 'cubic';
ceccsi.PilotAverage = 'UserDefined';
ceccsi.InterpWinSize = 1;
ceccsi.InterpWindow = 'Causal';
ceccsi.Reference = 'CSIRS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Channel estimation configuration
% cec.PilotAverage = 'UserDefined'; % Type of pilot averaging 
% cec.FreqWindow = 13;              % Frequency averaging windows in REs
% cec.TimeWindow = 1;               % Time averaging windows in REs
% cec.InterpType = 'cubic';         % Interpolation type
% cec.Reference = 'Antennas';       % Reference for channel estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % DM-RS estimation
cecdmrs = ceccsi;
cecdmrs.Reference = 'DMRS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NUsers = 6;
NSLRB = 3;
NULRB = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enb.CSIRefP = SimCSIRS.CSIRefP;
    enb.CSIRSConfig = SimCSIRS.CSIRSConfig;
    enb.CSIRSPeriod = SimCSIRS.CSIRSPeriod;
    enb.NCSIID = SimCSIRS.NCSIID;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ZP CSI-RS resource for TP1
zp1 = '0000000000000000';                          % Assume all CSI configurations unmuted
zp1(SimCSIIM.ZeroPowerCSIRSConfig([1 2])+1) = '1'; % Mute CSI-IM #0,1 (background & TP2 interference)
zp1(SimCSIRS.CSIRSConfig(2)+1) = '1';              % Mute CSI-RS #1 (TP2 transmission)

% Add ZP CSI-RS resource to TP1 parameters
enb.ZeroPowerCSIRSConfig = zp1;
enb.ZeroPowerCSIRSPeriod = SimCSIPeriod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ZP CSI-RS resource for TP2
zp2 = '0000000000000000';                          % Assume all CSI configurations unmuted
zp2(SimCSIIM.ZeroPowerCSIRSConfig([1 3])+1) = '1'; % Mute CSI-IM #0,2 (background & TP1 interference)
zp2(SimCSIRS.CSIRSConfig(1)+1) = '1';              % Mute CSI-RS #0 (TP1 transmission)

% Add ZP CSI-RS resource to TP2 parameters
enb.ZeroPowerCSIRSConfig = zp2;
enb.ZeroPowerCSIRSPeriod = SimCSIPeriod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%UE Parameterization for CSI Estimation and Reporting
numCSIRS = numel(SimCSIRS.CSIRSConfig);
csirs = repmat(enb,numCSIRS,1);
for idx = 1:numCSIRS
    csirs(idx).CSIRefP = SimCSIRS.CSIRefP(idx);
    csirs(idx).CSIRSConfig = SimCSIRS.CSIRSConfig(idx);
    csirs(idx).CSIRSPeriod = SimCSIRS.CSIRSPeriod{idx};
    csirs(idx).NCSIID = SimCSIRS.NCSIID(idx);
end

numCSIIM = numel(SimCSIIM.ZeroPowerCSIRSConfig);
csiim = repmat(enb,numCSIIM,1);
for idx = 1:numCSIIM
    csiim(idx).ZeroPowerCSIRSConfig = SimCSIIM.ZeroPowerCSIRSConfig(idx);
    csiim(idx).ZeroPowerCSIRSPeriod = SimCSIIM.ZeroPowerCSIRSPeriod{idx};
    csiim(idx).CSIRSPeriod = 'Off';
end

SINRs = [1.3 1.3 2.3 3.7 5 6.8 9.2 10.9 13 14.8 17.1 18.9 21 23.9 24.3];
numCSIProcesses = numel(SimCSIProcess.CSIRSResource);
process = repmat(enb,numCSIProcesses,1);
for idx = 1:numCSIProcesses
    % Index CSI-RS and CSI-IM resources used by the process
    process(idx).CSIRSIdx = SimCSIProcess.CSIRSResource(idx);
    process(idx).CSIIMIdx = SimCSIProcess.CSIIMResource(idx);
 % Reporting configuration
    process(idx).PDSCH.CSIMode = SimCSIProcess.CSIMode{idx};
    process(idx).PDSCH.PMIMode = SimCSIProcess.PMIMode{idx};
    process(idx).PDSCH.CodebookSubset = SimCSIProcess.CodebookSubset{idx};
    process(idx).PDSCH.SINRs90pc = SINRs;

    % CSI-RS configuration for CSI estimation
    process(idx).CSIRefP = SimCSIRS.CSIRefP(process(idx).CSIRSIdx);
    process(idx).CSIRSConfig = SimCSIRS.CSIRSConfig(process(idx).CSIRSIdx);
    process(idx).CSIRSPeriod = SimCSIRS.CSIRSPeriod(process(idx).CSIRSIdx);
    process(idx).NCSIID = SimCSIRS.NCSIID(process(idx).CSIRSIdx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation Setup
% Convert to linear and Watts
nocLin = 10.^(Noc/10)*(1e-3); % linear in Watts
% Take into account number of antennas and FFT (OFDM) scaling
No = sqrt(nocLin/(2*double(ofdmInfo.Nfft)));
NocW = 10.^((Noc-30)/10); % convert to W/15kHz

% SINR = Es/Noc TS 36.101 Sec. 8.1.1
NocTot = NocW; % W/15kHz
snrLin = 10.^(snrTP/10);
Es = snrLin*NocTot; % W/15kHz

% Amplitude scaling factors
K = sqrt(Es);

% Set the random number seed
rng('default');

% Initialize containers which store the channel estimates interference
% measurements for CSI-RS and CSI-IM resources
csirshest = cell(numCSIRS,1);
csiimnest = zeros(numCSIIM,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize a buffer to store CQI reports for each process. Size buffer
% for subband CQI reporting.
numCSIReports = ceil(totSubframes/SimCSIPeriod(1));
pmiInfo = ltePMIInfo(process(1),setfield(process(1).PDSCH,'PMIMode','Subband')); %#ok<SFLD>
cqiBuffer = ones(numCSIReports,pmiInfo.NSubbands+1,numCSIProcesses);

% Initialize buffers to store the CRC, BER and TP selected
crcBuffer = cell(totSubframes,1);
berBuffer =  zeros(totSubframes,2);
tpBuffer = zeros(totSubframes,1);


csiReportIdx = 1; % Index of CSI report
lastOffset = 0;   % Initialize overall frame timing offset
frameOffset = 0;  % Initialize frame timing offset

%Simulation Loop

for nsf = 0:totSubframes-1
    % Initialize UE receive waveform with noise
    rxWaveform = No*complex(randn(rxWaveformSize),randn(rxWaveformSize));

    % Select PDSCH TP based on the highest reported wideband CQI
    bufferIdx = mod(csiReportIdx-2,numCSIReports)+1; % Buffered CQI to use
    if dpsOperation
        widebandCQI = permute(cqiBuffer(:,1,:),[1 3 2]);
        if (widebandCQI(bufferIdx,3)>=widebandCQI(bufferIdx,4))
            pdschTransmissionPoint = 1; % TP1
        else
            pdschTransmissionPoint = 2; % TP2
        end
    else
         pdschTransmissionPoint = 1; %%%%%%%%%%%%%#ok<UNRCH> % TP1
    end
% Generate TP1 and TP2 transmissions.
    % Each transmission contains cell-specific reference signal,
    % synchronizing signals, CSI-RS and TM10 OCNG. One transmission
    % contains the PDSCH for the UE.
   
        % Update subframe number for each transmission
        enb.NSubframe = nsf;

        % Turn off the CSI resource not transmitted by this TP
        tpenb = enb;
        tpenb.CSIRSPeriod{mod(1,2)+1} = 'Off';

        % Blank subframe; CRS, PSS and SSS
        sf = lteResourceGrid(tpenb,tpenb.PDSCH.NTxAnts);
        crsInd = lteCellRSIndices(tpenb); % Cell-specific reference signal
        sf(crsInd) = lteCellRS(tpenb);
        pssInd = ltePSSIndices(tpenb);    % Primary synchronizing signal
        sf(pssInd) = ltePSS(tpenb);
        sssInd = lteSSSIndices(tpenb);    % Secondary synchronizing signal
        sf(sssInd) = lteSSS(tpenb);

        % CSI-RS resource
        csitp = tpenb;
        csitp.ZeroPowerCSIRSPeriod = 'Off';
        csiInd = lteCSIRSIndices(csitp);
        sf(csiInd) = lteCSIRS(csitp);

        % TM10 OCNG transmitted apart from subframes containing PSS/SSS/PBCH
        if isempty(pssInd)
            tpenb.PDSCH.RNTI = 0;

            % Add OCNG at DMRS locations
            oncngInd = lteDMRSIndices(tpenb,tpenb.PDSCH);
            ocngSym = lteDMRS(tpenb,tpenb.PDSCH);
            sf(oncngInd) = ocngSym;

            % Add OCNG for PDSCH symbols
            [oncngInd,ocngInfo] = ltePDSCHIndices(tpenb,tpenb.PDSCH,tpenb.PDSCH.PRBSet);
            ocngSym = ltePDSCH(tpenb,tpenb.PDSCH,randi([0 1],ocngInfo.G,1));
            sf(oncngInd) = ocngSym;
        end
        % PDSCH and DMRS Transmission to UE from either TP1 or TP2.
        % The transport block size and modulation scheme for transmission
        % are selected using the reported CQI for the appropriate CSI
        % process. The PDSCH must be mapped around the ZP CSI-RS configured
        % for the UE and the ZP CSI-RS of the TP. Only transmit PDSCH for
        % subframes not containing CSI resources or PSS/SSS as per TS36.101
        % Table 9.3.6.1-1.
        if ~isempty(pssInd)||~isempty(csiInd)
            tbs = 0; % Transport block size is 0 as no PDSCH transmitted
        elseif (pdschTransmissionPoint==1)
            tpBuffer(nsf+1) = pdschTransmissionPoint;
            % Get relevant configuration for transmission point and CQI to
            % use.
            txenb = enb;
            cqi = cqiBuffer(bufferIdx,:,1+2);

            % Select subband for PDSCH transmission.
            % The PDSCH is transmitted in the highest differential CQI
            % subband. Subbands less than full size are excluded.
            partialSubband = (pmiInfo.k*pmiInfo.NSubbands>txenb.NDLRB);
            [~,idx] = max(cqi(2:(end-partialSubband))); % Maximum differential
            sbCandidates = find(cqi(2:(end-partialSubband))==cqi(idx+1));
            sb = sbCandidates(randi([1 numel(sbCandidates)],1,1));
            cqi = cqi(1)+cqi(1+sb); % Calculate SB PMI from wideband and differential
            txenb.PDSCH.PRBSet = ((sb-1)*pmiInfo.k+(0:(pmiInfo.k-1))).'; % PRB allocation for subband

            % Select MCS according to CQI using TS36.101 Table A.4-1 CSI
            % RMC RC.12 FDD (MCS.13), which defines the relationship
            % between CQI indices and MCS indices
            imcsTable = [-1 0 0 1 3 5 7 10 12 14 17 19 21 22 24 25];
            imcs = imcsTable(cqi+1);

            % Determine TBS and modulation order, fixed RI (1) and PMI (0).
            % Generate PDSCH for only a non-zero transport block size.
            [itbs,modulation] = lteMCS(imcs);
            tbs = double(lteTBS(size(txenb.PDSCH.PRBSet,1),itbs));
            if any(tbs)
                if ~iscell(modulation)
                    modulation = {modulation};
                end
                txenb.PDSCH.NLayers = 1;
                txenb.PDSCH.Modulation = modulation;
                txenb.PDSCH.W = lteCSICodebook(txenb.PDSCH.NLayers,txenb.PDSCH.NTxAnts,0).';

                % PDSCH mapping
                [pdschInd,pdschInfo] = ltePDSCHIndices(txenb,txenb.PDSCH,txenb.PDSCH.PRBSet);

                % Generate DL-SCH data
                txtrblk = arrayfun(@(x)randi([0 1],x,1),tbs,'UniformOutput',false);
                cw = lteDLSCH(txenb,txenb.PDSCH,pdschInfo.G,txtrblk);

                % Generate PDSCH symbols with cell identity of serving cell
                % for correct scrambling
                txenb.NCellID = enb.NCellID;
                pdschSym = ltePDSCH(txenb,txenb.PDSCH,cw);

                % Create UE specific DMRS configuration to allow for
                % scrambling code to change depending on transmission point
                txenb.NCellID = enb.NCellID;
                if pdschTransmissionPoint == 2
                    txenb.NCellID = enb(2).NCellID;
                end
                dmrsInd = lteDMRSIndices(txenb,txenb.PDSCH);
                dmrsSym = lteDMRS(txenb,txenb.PDSCH);

                % Map PDSCH and DMRS
                sf(pdschInd) = pdschSym;
                sf(dmrsInd) = dmrsSym;
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Throughput Results
% The throughput results for all users are displayed in the MATLAB(R) command
% window after the simulation for each SNR point is completed. They are also 
% captured in output arrays |simThroughput| and |maxThroughput|. 

% legendString = cell(NUsers,1);
% figure;
% for userIdx = 1:NUsers
%     plot(SNRIn, simThroughput(:,userIdx)*100./maxThroughput(:,userIdx),'*-.');
%     hold on;
%     legendString{userIdx} = strcat('UE-' ,num2str(userIdx), ': ', ...
%         num2str(muNumLayers(userIdx)), ' layer(s), ' ,num2str(NTxAnts), ...
%         ' TxAnt(s), ', num2str(muNumRxAnts(userIdx)), ' RxAnt(s)');
% end
% grid on;
% xlabel('SNR (dB)');
% ylabel('Throughput (%)');
% legend(legendString,'Location','SouthEast');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%