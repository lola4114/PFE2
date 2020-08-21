% ue and enb caracteristics
dpsOperation = true; % Enable DPS {true,false}
totSubframes = 150;  % Number of subframes to simulate

% Transmission point 1 cell-wide settings
enb = struct;
enb.NDLRB = 6;
enb.CellRefP = 2;
enb.DuplexMode = 'FDD';
enb.CFI = 3;
enb.CyclicPrefix = 'Normal';
enb.NFrame = 0;
enb.NCellID = 0;

% PDSCH configuration for transmission mode 10 (TM10)
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

enb = repmat(enb,2,1);
enb(2).NCellID = 6;
enb(2).PDSCH.NTxAnts = 2;
enb(2).PDSCH.W = lteCSICodebook(enb(2).PDSCH.NLayers,enb(2).PDSCH.NTxAnts,0).';

SimCSIPeriod = [5 1]; % [Tcsi-rs Dcsi-rs]

% CSI-RS resource: {CSI-RS #0, CSI-RS #1}
SimCSIRS.CSIRSConfig = [0 5];                       % CSI-RS configuration
SimCSIRS.CSIRSPeriod = {SimCSIPeriod,SimCSIPeriod}; % CSI-RS period
SimCSIRS.NCSIID = [10 16];                          % CSI-RS scrambling identity
SimCSIRS.CSIRefP = [enb(1).PDSCH.NTxAnts enb(2).PDSCH.NTxAnts];

SimCSIProcess.CSIRSResource  = [1 2 1 2]; % CSI-RS resources index (1 based)
SimCSIProcess.CSIIMResource  = [1 1 2 3]; % CSI-IM resources index (1 based)
SimCSIProcess.CSIMode        = {'PUCCH 1-1','PUSCH 3-1','PUSCH 3-1','PUSCH 3-1'}; % CSI reporting modes
SimCSIProcess.PMIMode        = {'Wideband' ,'Wideband' ,'Wideband','Wideband'};   % PMI reporting modes
SimCSIProcess.CodebookSubset = {'0x0000000000000001','000001','0x0000000000000001','000001'}; % Codebook subset restrictions

%Fading Channel and SNR Configuration
snrTP = [11 8]; % SNR of received transmission from TP1 and TP2
Noc = -98;      % dBm/15kHz average power spectral density

% Channel between TP1 and UE
ofdmInfo = lteOFDMInfo(enb(1));
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
chcfg.Seed = 1;

% Channel between TP2 and UE
chcfg = repmat(chcfg,2,1);
chcfg(2).Seed = 2;

% Calculate the size of the received waveform
rxWaveformSize = [chcfg(1).SamplingRate*1e-3+15 chcfg(1).NRxAnts];

%Channel Estimation Configuration
% CSI-RS estimation
ceccsi.FreqWindow = 1;
ceccsi.TimeWindow = 2;
ceccsi.InterpType = 'cubic';
ceccsi.PilotAverage = 'UserDefined';
ceccsi.InterpWinSize = 1;
ceccsi.InterpWindow = 'Causal';
ceccsi.Reference = 'CSIRS';

% DM-RS estimation
cecdmrs = ceccsi;
cecdmrs.Reference = 'DMRS';

%Transmission Point Parameterization for CSI Estimation and Reporting
% Set TP1 and TP2 CSI-RS configuration
for enbIdx = 1:2
    enb(enbIdx).CSIRefP = SimCSIRS.CSIRefP;
    enb(enbIdx).CSIRSConfig = SimCSIRS.CSIRSConfig;
    enb(enbIdx).CSIRSPeriod = SimCSIRS.CSIRSPeriod;
    enb(enbIdx).NCSIID = SimCSIRS.NCSIID;
end

% ZP CSI-RS resource for TP1
zp1 = '0000000000000000';                          % Assume all CSI configurations unmuted
zp1(SimCSIIM.ZeroPowerCSIRSConfig([1 2])+1) = '1'; % Mute CSI-IM #0,1 (background & TP2 interference)
zp1(SimCSIRS.CSIRSConfig(2)+1) = '1';              % Mute CSI-RS #1 (TP2 transmission)

% Add ZP CSI-RS resource to TP1 parameters
enb(1).ZeroPowerCSIRSConfig = zp1;
enb(1).ZeroPowerCSIRSPeriod = SimCSIPeriod;

% ZP CSI-RS resource for TP2
zp2 = '0000000000000000';                          % Assume all CSI configurations unmuted
zp2(SimCSIIM.ZeroPowerCSIRSConfig([1 3])+1) = '1'; % Mute CSI-IM #0,2 (background & TP1 interference)
zp2(SimCSIRS.CSIRSConfig(1)+1) = '1';              % Mute CSI-RS #0 (TP1 transmission)

% Add ZP CSI-RS resource to TP2 parameters
enb(2).ZeroPowerCSIRSConfig = zp2;
enb(2).ZeroPowerCSIRSPeriod = SimCSIPeriod;

%UE Parameterization for CSI Estimation and Reporting
numCSIRS = numel(SimCSIRS.CSIRSConfig);
csirs = repmat(enb(1),numCSIRS,1);
for idx = 1:numCSIRS
    csirs(idx).CSIRefP = SimCSIRS.CSIRefP(idx);
    csirs(idx).CSIRSConfig = SimCSIRS.CSIRSConfig(idx);
    csirs(idx).CSIRSPeriod = SimCSIRS.CSIRSPeriod{idx};
    csirs(idx).NCSIID = SimCSIRS.NCSIID(idx);
end

numCSIIM = numel(SimCSIIM.ZeroPowerCSIRSConfig);
csiim = repmat(enb(1),numCSIIM,1);
for idx = 1:numCSIIM
    csiim(idx).ZeroPowerCSIRSConfig = SimCSIIM.ZeroPowerCSIRSConfig(idx);
    csiim(idx).ZeroPowerCSIRSPeriod = SimCSIIM.ZeroPowerCSIRSPeriod{idx};
    csiim(idx).CSIRSPeriod = 'Off';
end

SINRs = [1.3 1.3 2.3 3.7 5 6.8 9.2 10.9 13 14.8 17.1 18.9 21 23.9 24.3];
numCSIProcesses = numel(SimCSIProcess.CSIRSResource);
process = repmat(enb(1),numCSIProcesses,1);
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
        pdschTransmissionPoint = 1; %#ok<UNRCH> % TP1
    end

    % Generate TP1 and TP2 transmissions.
    % Each transmission contains cell-specific reference signal,
    % synchronizing signals, CSI-RS and TM10 OCNG. One transmission
    % contains the PDSCH for the UE.
    for enbIdx = 1:2
        % Update subframe number for each transmission
        enb(enbIdx).NSubframe = nsf;

        % Turn off the CSI resource not transmitted by this TP
        tpenb = enb(enbIdx);
        tpenb.CSIRSPeriod{mod(enbIdx,2)+1} = 'Off';

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
        elseif (pdschTransmissionPoint==enbIdx)
            tpBuffer(nsf+1) = pdschTransmissionPoint;
            % Get relevant configuration for transmission point and CQI to
            % use.
            txenb = enb(enbIdx);
            cqi = cqiBuffer(bufferIdx,:,enbIdx+2);

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
                txenb.NCellID = enb(1).NCellID;
                pdschSym = ltePDSCH(txenb,txenb.PDSCH,cw);

                % Create UE specific DMRS configuration to allow for
                % scrambling code to change depending on transmission point
                txenb.NCellID = enb(1).NCellID;
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

        % OFDM modulate, pass through a fading channel, scale for SNR and
        % add to existing receive waveform
        txWaveform = lteOFDMModulate(tpenb,sf);
        txWaveform = [txWaveform; zeros(15,size(txWaveform,2))]; %#ok<AGROW>
        chcfg(enbIdx).InitTime = nsf/1e3;
        rxWaveform = rxWaveform + K(enbIdx)*lteFadingChannel(chcfg(enbIdx),txWaveform);
    end

    % Receiver Synchronization and OFDM demodulation
    % Synchronize using the PSS/SSS of the serving cell (TP1) and OFDM
    % demodulate
    if ~isempty(pssInd)
        frameOffset = lteDLFrameOffset(enb(1),rxWaveform);
        if (frameOffset > 25)
            frameOffset = lastOffset;
        end
        lastOffset = frameOffset;
    end
    rxWaveform = rxWaveform(1+frameOffset:end,:);
    rxsf = lteOFDMDemodulate(enb(1),rxWaveform);

    % Calculate CSI-RS Estimates
    % Generate channel estimates for CSI-RS resources configured at the UE.
    for idx = 1:numCSIRS
        % Calculate channel estimate
        csirs(idx).NSubframe = nsf;
        csirshest{idx} = lteDLChannelEstimate(csirs(idx),csirs(idx).PDSCH,ceccsi,rxsf);
    end

    % Calculate interference using CSI-IM
    % For each CSI-IM resource calculate the energy in resource elements.
    % This is the noise+interference estimate.
    for idx = 1:numCSIIM
        % Calculate noise and interference estimate
        csiim(idx).NSubframe = nsf;
        imIndices = lteCSIRSIndices(csiim(idx));
        imSym = lteExtractResources(imIndices,rxsf);
        csiimnest(idx) = mean(abs(imSym(:)).^2);
    end

    % CSI Process Reporting
    % When estimated CSI resource elements are not zero calculate CSI
    if ~isempty(csiInd)
        % For each CSI process calculate CSI feedback
        for idx = 1:numCSIProcesses
            % Update subframe number
            process(idx).NSubframe = nsf;

            % Extract CSI-RS estimate and CSI-IM estimate for process
            hest = csirshest{process(idx).CSIRSIdx};
            nest = csiimnest(process(idx).CSIIMIdx);

            % Calculate CQI/PMI/RI, condition CQI on PMI/RI selection
            [ri,PMISet] = lteRISelect(process(idx),process(idx).PDSCH,hest,nest);
            process(idx).PDSCH.PMISet = PMISet;
            process(idx).PDSCH.NLayers = ri;
            process(idx).PDSCH.NCodewords = min(ri,2);
            [cqi,sinrs] = lteCQISelect(process(idx),process(idx).PDSCH,hest,nest);
            cqiBuffer(csiReportIdx,1:numel(cqi),idx) = cqi;
        end

        % New CSI report
        csiReportIdx = csiReportIdx+1;
    end

    % PDSCH Demodulation
    % If the transport block size is not 0, a PDSCH exists to decode
    if any(tbs)
        % Estimate channel using DMRS. To use the correct scrambling
        % sequence for the DMRS use the configuration for the active TP.
        [dmrshest,dmrsnest] = lteDLChannelEstimate(txenb,txenb.PDSCH,cecdmrs,rxsf);

        % Extract PDSCH symbols from received grid and channel estimate
        [sym,symhest] = lteExtractResources(pdschInd,rxsf,dmrshest);

        % Scale the received symbols by the PDSCH power factor Rho and
        % decode the PDSCH with the CellID for the serving cell
        sym = sym*(10^(-txenb.PDSCH.Rho/20));
        txenb.NCellID = enb(1).NCellID;
        [cws,recsym] = ltePDSCHDecode(txenb,txenb.PDSCH,sym,symhest,dmrsnest);
        % Scale cws by 1/K(1) to avoid numerical issues with DLSCH decoding
        cws = cellfun(@(x) x*(1/K(1)),cws,'UniformOutput',false);
        [trblk,crc] = lteDLSCHDecode(txenb,txenb.PDSCH,tbs,cws,[]);

        % Store CRC and BER
        crcBuffer{nsf+1} = double(crc);
        berBuffer(nsf+1,:) = [sum(trblk{1}~=txtrblk{1}),numel(trblk{1})];
    end
end

%Conclusion and Results

%hCoMPResults(totSubframes,SimCSIPeriod,crcBuffer,berBuffer,cqiBuffer,tpBuffer);
