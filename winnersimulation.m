function res = winnersimulation()

%winner+B1 & winner+C2


addpath('C:\Users\PC\OneDrive\Documents\MATLAB\Examples\R2019a\lte');
%addpath('"C:\Program Files\MATLAB\R2019a\bin\"');

%s = rng(21);    % For repeatability
%v = 60;     %absolute speed in km/h
AA(1) = winner2.AntennaArray('UCA', 16, 0.3);
AA(2) = winner2.AntennaArray('UCA', 12, 0.3);
AA(3) = winner2.AntennaArray('UCA', 8,  0.3);
AA(4) = winner2.AntennaArray('UCA', 4,  0.05);

BSIdx    = {1}; % Index in antenna array inventory vector
%N = input ('add number of pedestrian user equipments : \n'); 
%V = input ('add number of vehicle user equipments : \n');
%N = 2;
%V = 4;

%for i = 1:1:N
   
    
        %pedestrian
%Resultat(i)   = 4;     % Index in antenna array inventory vector


%end

%vehicles
%for j = N+1:1:V
%Resultat(j) = 3;
%end 
    MSIdx    = [4 4 3 3 3 3]; 
 %cellularLinks = N;
numLinks = 6;               % Number of links
range    = 433*250; 
%seed = 5;
% Layout range (meters)
cfgLayout = winner2.layoutparset(MSIdx, BSIdx, numLinks, AA, range);
%cfgLayout.CenterFrequency = 6.00e9;   % carrier frequency in Hz
%cfgLayout.StreetWidth = 20;        %street width specified on 5GAA.org


    
   
    % 11 for C2 
%cfgLayout.PropagConditionVector (i)= 0;  % 0 for NLOS
    
    cfgLayout.Pairing = [1 1 1 1 1 1  ;
        2 3 4 5 6 7 ]; 
    
  
cfgLayout.ScenarioVector = [11 11 3 3 3 3 ];     % 3 for B1 
cfgLayout.PropagConditionVector=[0 0 1 1 1 1 ];  % 0 for NLOS
    
numBSSect = sum(cfgLayout.NofSect);
numMS = length(MSIdx);
    
cfgLayout.Stations(1).Pos(1:2)=[125,216.5];
cfgLayout.Stations(2).Pos(1:2)=[7,216.5];
cfgLayout.Stations(3).Pos(1:2)=[243,216.5];
cfgLayout.Stations(4).Pos(1:2)=[125,432];
%cfgLayout.Stations(5).Pos(1:2)=[130,428]; 
cfgLayout.Stations(5).Pos(1:2)=[248,230];           
%cfgLayout.Stations(7).Pos(1:2)=[247,220];       
cfgLayout.Stations(6).Pos(1:2)=[0,220];   
cfgLayout.Stations(7).Pos(1:2)=[125,5];       
%cfgLayout.Stations(10).Pos(1:2)=[130,1];        
%cfgLayout.Stations(11).Pos(1:2)=[4,200.5];             
          


% Randomly draw MS velocity
  cfgLayout.Stations(2).Velocity = rand(3,1) -0.5;
    cfgLayout.Stations(3).Velocity = rand(3,1)-0.5;
for i = 4:numMS
   cfgLayout.Stations(i).Velocity = 16 + 1.*rand(3,1);
end
%cfgLayout.Stations.Pos(1:2)=vehicle_distribution(v);
BSPos = cell2mat({cfgLayout.Stations(1:numBSSect).Pos});   
MSPos = cell2mat({cfgLayout.Stations(numBSSect+1:end).Pos}); 

scrsz = get(groot, 'ScreenSize');
figSize = min(scrsz([3,4]))/2.3;
figure('Position', [scrsz(3)*.5-figSize/2,scrsz(4)*.7-figSize/2,figSize,figSize]);
hold on; grid on;
hBS = plot(BSPos(1,:), BSPos(2,:), 'or');   % Plot BS
hMS = plot(MSPos(1,:), MSPos(2,:), 'xr');   % Plot MS
for linkIdx = 1:numLinks  % Plot links
    pairStn = cfgLayout.Pairing(:,linkIdx);
    pairPos = cell2mat({cfgLayout.Stations(pairStn).Pos});
    plot(pairPos(1,:), pairPos(2,:), '-b');
end
xlim([0 250]); ylim([0 433]);
xlabel('X Position (meters)'); ylabel('Y Position (meters)')
legend([hBS, hMS], 'BS', 'MS', 'location', 'northwest');

frameLen = 150;   % Number of samples to be generated

%cfgWim = winner2.wimparset;

%cfgWim.RandomSeed = 10; % For repeatability

cfgWim = winner2.wimparset;
cfgWim.NumTimeSamples= frameLen;
cfgWim.SampleDensity = 20;
cfgWim.IntraClusterDsUsed  = 'yes';
cfgWim.CenterFrequency     = 6e9;
cfgWim.UniformTimeSampling = 'no';
cfgWim.ShadowingModelUsed  = 'yes';
cfgWim.PathLossModelUsed   = 'yes';
cfgWim.RandomSeed          = 31415926;  % For repeatability


 chanCoef = winner2.wim(cfgWim,cfgLayout) ;

WINNERChan = comm.WINNER2Channel(cfgWim, cfgLayout);
chanInfo = info(WINNERChan);
disp(chanCoef);

txSig = cellfun(@(x) [ones(1,x);zeros(frameLen-1,x)], ...
   num2cell(chanInfo.NumBSElements)', 'UniformOutput', false);
%txsig = awgn (txSig, 10, 'measured');
% Pass impulse signal through each link
rxSig = WINNERChan(txSig);
 [~,pathgains] = step(WINNERChan, txSig);
figure('Position', [scrsz(3)*.3-figSize/2,scrsz(4)*.25-figSize/2,figSize,figSize]);
hold on; 
for linkIdx = 1:numLinks
    delay = chanInfo.ChannelFilterDelay(linkIdx);
    stem(((0:(frameLen-1))-delay)/chanInfo.SampleRate(linkIdx), ...
        abs(rxSig{linkIdx}(:,1)));
   
end
maxX = max((cell2mat(cellfun(@(x) find(abs(x) < 1e-8, 1, 'first'), ...
    rxSig.', 'UniformOutput', false)) - chanInfo.ChannelFilterDelay)./ ...
    chanInfo.SampleRate);
minX = -max(chanInfo.ChannelFilterDelay./chanInfo.SampleRate);
xlim([minX, maxX]);
xlabel('Time (s)'); ylabel('Magnitude');
legend('Link 1', 'Link 2', 'Link 3', 'Link 4', 'Link 5', 'Link 6');
title('Impulse Response at First Receive Antenna');

SA = dsp.SpectrumAnalyzer( ...
   'Name',         'Frequency response', ...
    'SpectrumType', 'Power density', ...
    'SampleRate',   chanInfo.SampleRate(3), ...
    'Position',     [scrsz(3)*.7-figSize/2,scrsz(4)*.25-figSize/2,figSize,figSize], ...
    'Title',        'Frequency Response', ...
    'ShowLegend',   true, ...
    'ChannelNames', {'Link 3','Link 4'});

SA(cell2mat(cellfun(@(x) x(:,1), rxSig(3:4,1)', 'UniformOutput', false)));

%rng(s); % Restore RNG

%for i = 1 : 1 : 6
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% pathgains part
    
   % x = [ones(1,numLinks); zeros(frameLen-1,numLinks)];
%[y,pathgains] = step(WINNERChan, txSig);
disp(pathgains);
res  = pathgains;


