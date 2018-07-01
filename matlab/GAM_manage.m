% GAM_manage.m
%  
% Called From: prepFATE.m % J. Stewart
  
% Description: Goal: to get all data in a format that is ready to import
% into R for SEM analysis. 
%  
% Outside Functions Called: 
% sem_mbari.m % J. Stewart
% sem_orca.m % J. Stewart

% uses variables created in c_makeSPUE.m, c_importPFEL.m, c_filterPFEL.m

% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 19-Oct-2011 13:59:15  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : GAM_manage.m 


semTrixMulti = []; % this is for MBARI analyses
semTrixMultiFull = []; % this is for combo analyses--comment/uncomment datasets to compare. Possiblilities: 
                        % mbari, juv_rockfish, orcawale
%% prepare for MBARI

prep_mbari % J. Stewart % this imports Hake/Myctophids from VARS and makes PPUE_sess and PPUE_fulldive. Also, OMZdepth and upwelling
datasetN = 'mbari';

Squid = DSER; % squid from c_makeER.m
DepthOMZ = depthOMZ_mbari; % OMZ depthfrom sem_mbari.m 
Temp25m = temp25m_mbari; % temperature at 25±5 m.
Prey = PreyDSER; % prey: hake,loligo,myctophidae
DayL = semDaylight_hr; % day length (hours)
UIwin = nansum(semUItrix,2); % upwelling from sem_mbari
NOI = NOItrix; % El Nino index NOI
DepthStrat = strat_mbari; % stratification
LatN = latmean;

% make the MBARI matrix for SEM in R 
% Create the matrix and then remove nans, as they can't be modeled and will make SEM fussy
semTrix = [unidivenum Squid DepthOMZ Temp25m UIwin NOI DayL DepthStrat LatN Prey SqMonth SqYear ContourDistKM ROV_ID sum(SqdepthbinscountDownSess)' DurDownSess];  

semTrix(find(isnan(semTrix(:,3))),:) = []; % get rid of any of the OMZdepths that are nan--clutters the sem
semTrix(find(isnan(semTrix(:,4))),:) = []; % get rid of any of the Temp25m that are nan--clutters the sem
semTrix(find(isnan(semTrix(:,8))),:) = []; % get rid of any of the stratification that are nan--improper dives.

%%% Uncomment to save a .txt file of MBARI data to run in R: FINAL VERSION
% outfile = [DirFATEdata 'gamTrix_' datasetN '.txt'];
% fout = fopen(outfile,'w');
% header = 'UniDiveNum Squid DepthOMZ Temp25m UIwin NOI DayL Strat LatN Hake Loligo Myct SqMonth SqYear ContourDistKM ROV_ID SquidCount Effort'; % write header to file
% fprintf(fout,'%s \n',header);
% fprintf(fout,'%4.1f\t%3.4f\t%3.4f\t%4.4f\t%4.0f\t%1.3f\t%2.3f\t%3.0f\t%2.3f\t%4.3f\t%4.3f\t%4.3f\t%2.0f\t%4.0f\t%4.3f\t%1.0f\t%4.3f\t%4.3f\n', semTrix.');

%semTrixMultiFull = [semTrixMultiFull; Squid DepthOMZ Temp25m UIwin NOI DayL DepthStrat LatN Prey repmat(1,length(LatN),1)];  % add ChlA!

% to add to historic table:
mean(semTrix(:,3))
std(semTrix(:,3))

d = depthOLZ_mbari;
d(find(isnan(d)),:) = [];
mean(d)
std(d)

%% prepare from ORCAWALE

prep_orca % by J. Stewart
datasetN = 'orca';

UniDiveNum = CastNumUni; % unique dive ID
Squid = ORSqEvents; % squid from sem_orca: count
DepthOMZ = depthOMZ_orca; 
Temp25m = temp25m_orca;
Prey = nan(length(DepthOMZ),3);
DayL = OR_Daylight_hr; 
UIwin = nansum(ORsemUItrix,2);
NOI = NOIorca; 
DepthStrat = strat_orca;
LatN = OrcaLat;
DistShelfKM = ContourDistKMorca;

semTrixMultiFull = [semTrixMultiFull; UniDiveNum Squid DepthOMZ Temp25m UIwin NOI DayL DepthStrat LatN DistShelfKM repmat(2,length(LatN),1)];  

%% prepare from JuvRockfish Cruise

prep_jrock % J. Stewart
datasetN = 'jrock';

UniDiveNum = CruiseHaulCTDUni;
Squid = JRockSqEventsJ;% JRockTrix(:,1); % CPUE
DepthOMZ = depthOMZ_jrock; % JRockTrix(:,2) this is filtered to remove the shelf
Temp25m = temp25m_jrock; % JRockTrix(:,3);
DayL = JRock_Daylight_hr; %JRockTrix(:,6); 
UIwin = nansum(JRocksemUItrix,2); %JRockTrix(:,4);
NOI = NOIjrock; %JRockTrix(:,5); 
DepthStrat = strat_jrock;% JRockTrix(:,7); % stratification
LatN = LatJuvCruiseHaul; % JRockTrix(:,8);
DistShelfKM = ContourDistKMjuv;

semTrixMultiFull = [semTrixMultiFull; UniDiveNum Squid DepthOMZ Temp25m UIwin NOI DayL DepthStrat LatN DistShelfKM repmat(3,length(LatN),1)]; 

%% NWFSC FRAM

prep_nwfscfram %%% need to continue this, it's working if I can get sem_calcs to go


prep_jrock % J. Stewart
datasetN = 'fram';

UniDiveNum = CruiseHaulCTDUni;
Squid = JRockSqEventsJ;% JRockTrix(:,1); % CPUE
DepthOMZ = depthOMZ_jrock; % JRockTrix(:,2) this is filtered to remove the shelf
Temp25m = temp25m_jrock; % JRockTrix(:,3);
DayL = JRock_Daylight_hr; %JRockTrix(:,6); 
UIwin = nansum(JRocksemUItrix,2); %JRockTrix(:,4);
NOI = NOIjrock; %JRockTrix(:,5); 
DepthStrat = strat_jrock;% JRockTrix(:,7); % stratification
LatN = LatJuvCruiseHaul; % JRockTrix(:,8);
DistShelfKM = ContourDistKMjuv;

semTrixMultiFull = [semTrixMultiFull; UniDiveNum Squid DepthOMZ Temp25m UIwin NOI DayL DepthStrat LatN DistShelfKM repmat(3,length(LatN),1)]; 

%% prepare from Brodeur

prep_saip % J. Stewart
datasetN = 'saip';

Squid = SAIP_hauldata(:,7);
DepthOMZ = SAIP_hauldata(:,4); % SAIP: these are all NaNs. All of them. 
Temp25m = SAIP_hauldata(:,5); % temperature at 25±5 m.
Prey = nan(length(DepthOMZ),3); 
DayL = SAIP_Daylight_hr; % problem with no day
UIwin = nansum(semUItrix,2); % problem here because no day data associated with dates
NOI = NOIsaip; % El Nino index NOI. monthly.
DepthStrat = SAIP_hauldata(:,6); % stratification
LatN = SAIP_hauldata(:,3);

% data aren't really up for comparison across groups: no oxygen data, no dates

%% prepare from canada. 
sem_canadian % J. Stewart
datasetN = 'canada';

% Canadian data is not in a format that can easily be processed, and no OMZ data. 

%% MultigroupFull save

semTrixMultiFull2 = semTrixMultiFull;

semTrixMultiFull2(find(isnan(semTrixMultiFull2(:,3))),:) = []; % get rid of any of the Temp25m that are nan--clutters the sem
semTrixMultiFull2(find(isnan(semTrixMultiFull2(:,4))),:) = []; % get rid of any of the Temp25m that are nan--clutters the sem
semTrixMultiFull2(find(isnan(semTrixMultiFull2(:,8))),:) = []; % get rid of any of the stratification that are nan--improper dives.

%%% uncomment to save the file
% outfile = [DirFATEdata 'semTrixMultiFull.txt'];
% fout = fopen(outfile,'w');
% header = 'UniDiveNum Squid DepthOMZ Temp25m UIwin NOI DayL Strat LatN DistShelfKM DataGroup'; % write header to file
% fprintf(fout,'%s \n',header);
% fprintf(fout,'%3.2f\t%3.4f\t%3.4f\t%4.4f\t%4.0f\t%1.3f\t%2.3f\t%3.0f\t%2.3f\t%4.4f\t%1.0f\n', semTrixMultiFull2.');

% save CTDfile_multi, which is made in prep_orca.m and prep_jrock.m, to calculate variables in R
%%% uncomment to save file
% outfile = [DirFATEdata 'CTDMultiFull.txt'];
% fout = fopen(outfile,'w');
% header = 'DiveNum Depth Oxyg Temper DataGroup'; % write header to file
% fprintf(fout,'%s \n',header);
% fprintf(fout,'%3.2f\t%4.3f\t%3.3f\t%2.3f\t%1.0f\n', CTDfile_multi.');




%% 
  
disp('Completed GAM_manage.m') 
% ===== EOF [GAM_manage.m] ======  
