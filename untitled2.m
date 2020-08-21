% configuration du canal
 chcfg.DelayProfile = 'EPA';
   chcfg.NRxAnts = 1;
   chcfg.DopplerFreq = 5;
   chcfg.MIMOCorrelation = 'Low';
   chcfg.Seed = 1;
   chcfg.InitPhase = 'Random';
   chcfg.ModelType = 'GMEDS';
   chcfg.NTerms = 16;
   chcfg.NormalizeTxAnts = 'On';
   chcfg.NormalizePathGains = 'On';


%configuration de l'eNB
enb.NDLRB = 9;
enb.CyclicPrefix = 'Normal';
enb.PHICHDuration = 'Normal';
enb.CFI = 3;
enb.Ng = 'Sixth';
enb.CellRefP = 1;
enb.NCellID = 1;
enb.NSubframe = 0;
enb.DuplexMode = 'FDD';

%configuration d'un équipement utilisateur