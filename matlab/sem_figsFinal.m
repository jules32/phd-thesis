% sem_figsFinal.m
%  
% Called From: 
  
% Description: 
%  
% Outside Functions Called: 
%  
% AUTHOR    : J. Stewart jules32@gmail.com 
% DATE      : 28-Feb-2013 15:57:13  
% Revision  : 1.00  
% DEVELOPED : 7.9.0.529 (R2009b) OSX 
% FILENAME  : sem_figsFinal.m 

Mfilename = mfilename; 
  
%% semTrix = [1unidivenum 2Squid 3DepthOMZ 4Temp25m 5UIwin 6NOI 7DayL 8DepthStrat 9LatN 10Prey 
% 13SqMonth 14SqYear 15ContourDistKM 16ROV_ID 17sum(SqdepthbinscountDownSess)' 18DurDownSess' 19PDOtrix
% 20NPGOtrix 21lonmean];  % add ChlA!

semTrixSqYrUni = unique(semTrix(:,14));
semTrixSqMoUni = unique(semTrix(:,13));
  
%% Fig 2: bubble plot REDO WITH ONLY APPROPRIATE DATA

% Squid SPUE Monthly Table
% Sq and Prey month-year ones were made in sem_prey: semTrixSqYMTot, semTrixMyYMTot etc
ScatterMatrixA = [];
for y = 1:length(semTrixSqYrUni)
    ylog = semTrix(:,14) == semTrixSqYrUni(y);
    for m = 1:length(semTrixSqMoUni)
        mlog = semTrixSqMoUni(m) == semTrix(:,13);
        my = ylog + mlog;
        mylog = my == 2;
        SPDtemp = semTrix(mylog,2);
        ScatterMatrixA = [ScatterMatrixA; semTrixSqYrUni(y) semTrixSqMoUni(m) length(SPDtemp)...
            length(SPDtemp(SPDtemp>0)) (length(SPDtemp(SPDtemp>0))/length(SPDtemp))*100]; 
    end
end

ScatterMatrix = ScatterMatrixA(ScatterMatrixA(:,1)>1996,:);
ScatterMatrix(isnan(ScatterMatrix(:,5)),:) = [];

figure; hold on;
colormap(flipud(gray))
scatter(ScatterMatrix(:,1), ScatterMatrix(:,2), ScatterMatrix(:,3)*10, 'k')
scatter(ScatterMatrix(ScatterMatrix(:,5)>0,1), ScatterMatrix(ScatterMatrix(:,5)>0,2), ...
    ScatterMatrix(ScatterMatrix(:,5)>0,3)*10, ScatterMatrix(ScatterMatrix(:,5)>0,5), 'filled')
set(gca, 'XLim', [min(unique(ScatterMatrix(:,1)))-1 max(unique(ScatterMatrix(:,1)))])
xlabel('SqYear ', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('SqMonth ', 'FontSize', 14, 'FontWeight', 'bold')
% title('ROV dives and Humboldt squid occurrence',...
%     'FontSize', 16, 'FontWeight', 'bold')
fixfig
FigName = ['Dive_Squid_BubbleplotFiltered.pdf'];
annotate_JS(Mfilename, gcf, FigName)
% print('-dpdf', [DirFATEsem '/_figs/' FigName]);

%% do it without the crazy way. 

yrcutoff = 3; % this separates the year at March: so a year goes from March-Feb

% semTrix = [1unidivenum 2Squid 3DepthOMZ 4Temp25m 5UIwin 6NOI 7DayL 8DepthStrat 9LatN 10Prey 
% 13SqMonth 14SqYear 15ContourDistKM 16ROV_ID 17sum(SqdepthbinscountDownSess)' 18DurDownSess' 19PDOtrix
% 20NPGOtrix 21lonmean];  % add ChlA!

% this is all with dateCTDV, which is all data, but since it is dealing
% with OMZ values it will self-select. 
dateCTDSqYr = dateCTDV(:,1);
YrInd = find(dateCTDV(:,2) < yrcutoff); 
dateCTDSqYr(YrInd) = dateCTDSqYr(YrInd)-1;
dateCTDSqMonth = dateCTDV(:,2); % and change Jan-Feb to from 1,2 to 13,14 just to make it different
MInd1 = find(dateCTDV(:,2) == 1);
MInd2 = find(dateCTDV(:,2) == 2);
dateCTDSqMonth(MInd1) = 13;
dateCTDSqMonth(MInd2) = 14;

OMZmean_stdY = nan(length(semTrixSqYrUni),2);
OMZmedian_q_nY = nan(length(semTrixSqYrUni),3);
OLZmedian = nan(length(semTrixSqYrUni),1);

Sq_median_q_nY = nan(length(semTrixSqYrUni),4);
Ha_median_q_nY = nan(length(semTrixSqYrUni),4);
Lo_median_q_nY = nan(length(semTrixSqYrUni),4);
My_median_q_nY = nan(length(semTrixSqYrUni),4);

% values by year
for i=1:length(semTrixSqYrUni)
    
    % for OMZ--only using Midwater 1 dives between 36.5, 37 and 121.75,122.3
    plog = dateCTDSqYr == semTrixSqYrUni(i);
    depthp = depth(plog);
    oxygp = oxyg(plog);
    latp = lat(plog);
    lonp = lon(plog);
    omzdef = oxygp >= 0.49 & oxygp < 0.51 & depthp < 700 & latp >=36.5 & latp < 37 & lonp >-122.3 & lonp <121.75;
    olzdef = oxygp >= 1.49 & oxygp < 1.51 & depthp < 700 & latp >=36.5 & latp < 37 & lonp >-122.3 & lonp <121.75;

    
    OMZmean_stdY(i,1) = mean(depthp(omzdef));
    OMZmean_stdY(i,2) = std(depthp(omzdef));
    OMZmedian_q_nY(i,1) = median(depthp(omzdef));
    OMZmedian_q_nY(i,2) = quantile(depthp(omzdef), 0.25);
    OMZmedian_q_nY(i,3) = quantile(depthp(omzdef), 0.75);
    OMZmedian_q_nY(i,4) = length(depthp(omzdef));
    OLZmedian(i) = median(depthp(olzdef));
    
    % for squid
    Slog = sqcountdatesallV(:,1) == semTrixSqYrUni(i);
    
    Sq_median_q_nY(i,1) = nanmedian(sqcountdepthsall(Slog));
    Sq_median_q_nY(i,2) = quantile(sqcountdepthsall(Slog), 0.25);
    Sq_median_q_nY(i,3) = quantile(sqcountdepthsall(Slog),0.75);
    Sq_median_q_nY(i,4) = length(sqcountdepthsall(Slog));
    
    % for hake
    Hlog = hacountdepthsall(:,3) == semTrixSqYrUni(i);
    
    Ha_median_q_nY(i,1) = nanmedian(hacountdepthsall(Hlog,1));
    Ha_median_q_nY(i,2) = quantile(hacountdepthsall(Hlog,1), 0.25);
    Ha_median_q_nY(i,3) = quantile(hacountdepthsall(Hlog,1),0.75);
    Ha_median_q_nY(i,4) = length(hacountdepthsall(Hlog,1));
  
  % for loligo
    Llog = locountdepthsall(:,3) == semTrixSqYrUni(i);
    
    Lo_median_q_nY(i,1) = nanmedian(locountdepthsall(Llog,1));
    Lo_median_q_nY(i,2) = quantile(locountdepthsall(Llog,1), 0.25);
    Lo_median_q_nY(i,3) = quantile(locountdepthsall(Llog,1),0.75);
    Lo_median_q_nY(i,4) = length(locountdepthsall(Llog,1));
    
    % for myct
    Mlog = mycountdepthsall(:,3) == semTrixSqYrUni(i);
    
    My_median_q_nY(i,1) = nanmedian(mycountdepthsall(Mlog,1));
    My_median_q_nY(i,2) = quantile(mycountdepthsall(Mlog,1), 0.25);
    My_median_q_nY(i,3) = quantile(mycountdepthsall(Mlog,1),0.75);
    My_median_q_nY(i,4) = length(mycountdepthsall(Mlog,1));
    
end


% values by unidivenum

Sq_median_dnum = nan(length(unidivenum),1);
Lo_median_dnum = nan(length(unidivenum),1);
Ha_median_dnum = nan(length(unidivenum),1);
My_median_dnum = nan(length(unidivenum),1);

% for i=1:length(unidivenum)
%     
%     % for squid
%     Slog = sqcountdivenumall == unidivenum(i);
%     Sq_median_dnum(i) = nanmedian(sqcountdepthsall(Slog));
%     
%     % for loligo
%     Llog = locountdepthsall(:,2) == unidivenum(i);
%     Lo_median_dnum(i) = nanmedian(locountdepthsall(Llog,1));
%    
%    % for hake
%     Hlog = hacountdepthsall(:,2) == unidivenum(i);
%     Ha_median_dnum(i) = nanmedian(hacountdepthsall(Hlog,1));
%    
%     % for myct
%     Mlog = mycountdepthsall(:,2) == unidivenum(i);
%     My_median_dnum(i) = nanmedian(mycountdepthsall(Mlog,1));
% 
% end

%% Now those figures

% figure; hold on
% plot(semTrixSqYrUni, Sq_median_q_nY(:,1), 'k')
% plot(semTrixSqYrUni, Sq_median_q_nY(:,2), '--k')
% plot(semTrixSqYrUni, Sq_median_q_nY(:,3), '--k')
% plot(semTrixSqYrUni, My_median_q_nY(:,1), 'b')
% plot(semTrixSqYrUni, My_median_q_nY(:,2), '--b')
% plot(semTrixSqYrUni, My_median_q_nY(:,3), '--b')
% fixfig
% xlabel('Year ', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('Depth (m) ', 'FontSize', 14, 'FontWeight', 'bold')
% % title('Black = mean squid depths, grey = prey. dark shading = OMZ, light shading = OLZ')
% set(gca, 'YDir', 'reverse', 'XLim', [semTrixSqYrUni(2) max(semTrixSqYrUni)], 'YLim', [0 1000])
% FigName = ['DepthO2_sq_prey.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% % print('-dpdf', [DirFATEsem '_figs/' FigName]);
% 
% % need to make all of these figures and then put them over each other in Matlab
% AllnY = log([Sq_median_q_nY(:,4); Ha_median_q_nY(:,4); Lo_median_q_nY(:,4); My_median_q_nY(:,4)]);
% 
% figure; hold on % squid
% h(1) = area(semTrixSqYrUni, repmat(1000,length(OMZmedian_q_nY(:,1)),1));
% h(2) = area(semTrixSqYrUni, OMZmedian_q_nY(:,1));
% set(h(1),'FaceColor',[1 0.85 1]) 
% set(h(2),'FaceColor',[1 1 1])
% cm = makeColorMap([1,1,1], [1,0,0], length(semTrixSqYrUni)); %makeColorMap from MatlabCentral
% colormap(cm)
% scatter(semTrixSqYrUni, Sq_median_q_nY(:,1), Sq_median_q_nY(:,3)-Sq_median_q_nY(:,2), Sq_median_q_nY(:,4), 'filled')
% scatter(semTrixSqYrUni, Sq_median_q_nY(:,1), Sq_median_q_nY(:,3)-Sq_median_q_nY(:,2), [1,0,0])
% plot(semTrixSqYrUni, OMZmedian_q_nY(:,1), 'k-')
% set(gca, 'YLim', [0 1000], 'XLim', [min(semTrixSqYrUni) max(semTrixSqYrUni)], 'YDir', 'reverse')
% fixfig
% xlabel('Year ', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('Depth (m) ', 'FontSize', 14, 'FontWeight', 'bold')
% FigName = ['DepthO2_Sq.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% % print('-dpdf', [DirFATEsem '_figs/' FigName]);
% 
% figure; hold on % myct
% cm = makeColorMap([1,1,1], [0,0,1], length(semTrixSqYrUni)); %makeColorMap from MatlabCentral
% colormap(cm)
% scatter(semTrixSqYrUni, My_median_q_nY(:,1), My_median_q_nY(:,3)-My_median_q_nY(:,2), My_median_q_nY(:,4), 'filled')
% scatter(semTrixSqYrUni, My_median_q_nY(:,1), My_median_q_nY(:,3)-My_median_q_nY(:,2), [0,0,1])
% set(gca, 'YLim', [0 1000], 'XLim', [min(semTrixSqYrUni)-1 max(semTrixSqYrUni)], 'YDir', 'reverse')
% fixfig
% FigName = ['DepthO2_My.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% % print('-dpdf', [DirFATEsem '_figs/' FigName]);
% 
% figure; hold on % hake
% cm = makeColorMap([1,1,1], [0,1,0], length(semTrixSqYrUni)); %makeColorMap from MatlabCentral
% colormap(cm)
% scatter(semTrixSqYrUni, Ha_median_q_nY(:,1), Ha_median_q_nY(:,3)-Ha_median_q_nY(:,2), Ha_median_q_nY(:,4), 'filled')
% scatter(semTrixSqYrUni, Ha_median_q_nY(:,1), Ha_median_q_nY(:,3)-Ha_median_q_nY(:,2), [0,1,0])
% set(gca, 'YLim', [0 1000], 'XLim', [min(semTrixSqYrUni)-1 max(semTrixSqYrUni)], 'YDir', 'reverse')
% fixfig
% FigName = ['DepthO2_Ha.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% % print('-dpdf', [DirFATEsem '_figs/' FigName]);
% 
% 
% figure; hold on % loligo
% cm = makeColorMap([1,1,1], [0,1,1], length(semTrixSqYrUni)); %makeColorMap from MatlabCentral
% colormap(cm)
% scatter(semTrixSqYrUni, Lo_median_q_nY(:,1), Lo_median_q_nY(:,3)-Lo_median_q_nY(:,2)+1, Lo_median_q_nY(:,4), 'filled')
% scatter(semTrixSqYrUni, Lo_median_q_nY(:,1), Lo_median_q_nY(:,3)-Lo_median_q_nY(:,2)+1, [0,1,1])
% set(gca, 'YLim', [0 1000], 'XLim', [min(semTrixSqYrUni)-1 max(semTrixSqYrUni)], 'YDir', 'reverse')
% fixfig
% FigName = ['DepthO2_Lo.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% % print('-dpdf', [DirFATEsem '_figs/' FigName]);
% 
% 
% plot(semTrixSqYrUni, Sq_median_q_nY(:,2), '--k')
% plot(semTrixSqYrUni, Sq_median_q_nY(:,3), '--k')
% plot(semTrixSqYrUni, My_median_q_nY(:,1), 'b')
% plot(semTrixSqYrUni, My_median_q_nY(:,2), '--b')
% plot(semTrixSqYrUni, My_median_q_nY(:,3), '--b')
% fixfig
% xlabel('Year ', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('Depth (m) ', 'FontSize', 14, 'FontWeight', 'bold')
% % title('Black = mean squid depths, grey = prey. dark shading = OMZ, light shading = OLZ')
% set(gca, 'YDir', 'reverse', 'XLim', [semTrixSqYrUni(2) max(semTrixSqYrUni)], 'YLim', [0 1000])
% FigName = ['DepthO2_sq_prey.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% % print('-dpdf', [DirFATEsem '_figs/' FigName]);

%% Boxplot subplots ACTUALLY THEY ARE ERROBAR PLOTS but this is how you'd do it: %

semTrixSqYMTot = zeros(12, length(semTrixSqYrUni));
semTrixHaYMTot = zeros(12, length(semTrixSqYrUni));
semTrixLoYMTot = zeros(12, length(semTrixSqYrUni));
semTrixMyYMTot = zeros(12, length(semTrixSqYrUni));

for y = 1:length(semTrixSqYrUni)
    ylog = semTrix(:,14) == semTrixSqYrUni(y);
    for m = 1:length(semTrixSqMoUni)
        mlog = semTrixSqMoUni(m) == semTrix(:,13);
        my = ylog + mlog;
        mylog = my == 2;
       
        % matrices
        semTrixSqYMTot(m,y) = sum(semTrix(mylog,2));
        semTrixHaYMTot(m,y) = sum(semTrix(mylog,10));
        semTrixLoYMTot(m,y) = sum(semTrix(mylog,11));
        semTrixMyYMTot(m,y) = sum(semTrix(mylog,12));
        
    end
end


figure
subplot(2,2,1); hold on % loligo
LoYTot = sum(semTrixLoYMTot);
h(1) = area(semTrixSqYrUni, repmat(1000,length(OMZmedian_q_nY(:,1)),1)); 
h(2) = area(semTrixSqYrUni, OMZmedian_q_nY(:,1));
h(3)= area(semTrixSqYrUni, OLZmedian);
set(h(1),'FaceColor',[0.6 0.6 0.6]) 
set(h(2),'FaceColor',[0.8 0.8 0.8])
set(h(3),'FaceColor',[1 1 1])
scatter(semTrixSqYrUni, Lo_median_q_nY(:,1), LoYTot/max(LoYTot)*100+0.1, 'c', 'filled')
errorbar(semTrixSqYrUni, Lo_median_q_nY(:,1), abs(Lo_median_q_nY(:,2)-Lo_median_q_nY(:,1)), abs(Lo_median_q_nY(:,3)-Lo_median_q_nY(:,1)), 'c', 'linewidth', 2)
text(semTrixSqYrUni(2), 800, ['max=' num2str(max(LoYTot)) ', min=' num2str(min(LoYTot(LoYTot>0)))])
set(gca, 'XLim', [min(semTrixSqYrUni) max(semTrixSqYrUni)], 'YLim', [0 1000], 'YDir', 'reverse') 

subplot(2,2,2); hold on % hake
HaYTot = sum(semTrixHaYMTot);
h(1) = area(semTrixSqYrUni, repmat(1000,length(OMZmedian_q_nY(:,1)),1)); 
h(2) = area(semTrixSqYrUni, OMZmedian_q_nY(:,1));
h(3)= area(semTrixSqYrUni, OLZmedian);
set(h(1),'FaceColor',[0.6 0.6 0.6]) 
set(h(2),'FaceColor',[0.8 0.8 0.8])
set(h(3),'FaceColor',[1 1 1])
scatter(semTrixSqYrUni, Ha_median_q_nY(:,1), HaYTot/max(HaYTot)*100+0.1, 'b', 'filled')
% scatter(semTrixSqYrUni, Ha_median_q_nY(:,1), log(HaYTot+1)*10, 'b', 'filled')
errorbar(semTrixSqYrUni, Ha_median_q_nY(:,1), abs(Ha_median_q_nY(:,2)-Ha_median_q_nY(:,1)), abs(Ha_median_q_nY(:,3)-Ha_median_q_nY(:,1)), 'b', 'linewidth', 2)
text(semTrixSqYrUni(2), 800, ['max=' num2str(max(HaYTot)) ', min=' num2str(min(HaYTot(HaYTot>0)))])
set(gca, 'XLim', [min(semTrixSqYrUni) max(semTrixSqYrUni)], 'YLim', [0 1000], 'YDir', 'reverse') 

subplot(2,2,3); hold on % myct
MyYTot = sum(semTrixMyYMTot);
h(1) = area(semTrixSqYrUni, repmat(1000,length(OMZmedian_q_nY(:,1)),1)); 
h(2) = area(semTrixSqYrUni, OMZmedian_q_nY(:,1));
h(3)= area(semTrixSqYrUni, OLZmedian);
set(h(1),'FaceColor',[0.6 0.6 0.6]) 
set(h(2),'FaceColor',[0.8 0.8 0.8])
set(h(3),'FaceColor',[1 1 1])
scatter(semTrixSqYrUni, My_median_q_nY(:,1), MyYTot/max(MyYTot)*100+0.1, 'r', 'filled')
errorbar(semTrixSqYrUni, My_median_q_nY(:,1), abs(My_median_q_nY(:,2)-My_median_q_nY(:,1)), abs(My_median_q_nY(:,3)-My_median_q_nY(:,1)), 'r', 'linewidth', 2)
text(semTrixSqYrUni(2), 800, ['max=' num2str(max(MyYTot)) ', min=' num2str(min(MyYTot(MyYTot>0)))])
set(gca, 'XLim', [min(semTrixSqYrUni) max(semTrixSqYrUni)], 'YLim', [0 1000], 'YDir', 'reverse') 

subplot(2,2,4); hold on % humboldt
SqYTot = sum(semTrixSqYMTot);
h(1) = area(semTrixSqYrUni, repmat(1000,length(OMZmedian_q_nY(:,1)),1)); 
h(2) = area(semTrixSqYrUni, OMZmedian_q_nY(:,1));
h(3)= area(semTrixSqYrUni, OLZmedian);
set(h(1),'FaceColor',[0.6 0.6 0.6]) 
set(h(2),'FaceColor',[0.8 0.8 0.8])
set(h(3),'FaceColor',[1 1 1])
scatter(semTrixSqYrUni, Sq_median_q_nY(:,1), SqYTot/max(SqYTot)*100+0.1, 'k', 'filled')
errorbar(semTrixSqYrUni, Sq_median_q_nY(:,1), abs(Sq_median_q_nY(:,2)-Sq_median_q_nY(:,1)), abs(Sq_median_q_nY(:,3)-Sq_median_q_nY(:,1)), 'k', 'linewidth', 2)
text(semTrixSqYrUni(2), 800, ['max=' num2str(max(SqYTot)) ', min=' num2str(min(SqYTot(SqYTot>0)))])
set(gca, 'XLim', [min(semTrixSqYrUni) max(semTrixSqYrUni)], 'YLim', [0 1000], 'YDir', 'reverse') 

fixfig
FigName = ['DepthO2_Subplot.pdf'];
annotate_JS(Mfilename, gcf, FigName)
orient('landscape')
% print('-dpdf', [DirFATEsem '_figs/' FigName]);

%% This is the boxplot fig

% figure; 
% subplot(2,2,1); hold on % humboldt
% boxplot(locountdepthsall(:,1),locountdepthsall(:,3),'colors', 'k','symbol', '.w', 'outliersize', 1); % but this wouldn't work with overlay of abundance
% set(gca, 'YLim', [0 1000], 'YDir', 'reverse') 
% 
% subplot(2,2,2); hold on % humboldt
% boxplot(hacountdepthsall(:,1),hacountdepthsall(:,3),'colors', 'k','symbol', '.w', 'outliersize', 1); % but this wouldn't work with overlay of abundance
% set(gca, 'YLim', [0 1000], 'YDir', 'reverse') 
% 
% subplot(2,2,3); hold on % humboldt
% boxplot(mycountdepthsall(:,1),mycountdepthsall(:,3),'colors', 'k','symbol', '.w', 'outliersize', 1); % but this wouldn't work with overlay of abundance
% set(gca,  'YLim', [0 1000], 'YDir', 'reverse') 
% 
% subplot(2,2,4); hold on % humboldt
% boxplot([sqcountdepthsall; 100; 100],[sqcountdatesallV(:,1); 2000; 2001],'colors', 'k','symbol', '.w', 'outliersize', 1); % but this wouldn't work with overlay of abundance
% set(gca,  'YLim', [0 1000], 'YDir', 'reverse') 
% FigName = ['DepthO2_SubplotBoxplot.pdf'];
% annotate_JS(Mfilename, gcf, FigName)
% orient('landscape')
% % print('-dpdf', [DirFATEsem '_figs/' FigName]);
% 
% 
% 
% 
% 
% % and associated regression
% figure; hold on % by unidivenum
% plot(depthOMZ_mbari,My_median_dnum, 'xr')
% plot(depthOMZ_mbari,Sq_median_dnum, 'ok')
% plot(depthOMZ_mbari,Lo_median_dnum, 'sc')
% plot(depthOMZ_mbari,Ha_median_dnum, '.b')
% ylabel('Median organism depth (m) ', 'FontSize', 14, 'FontWeight', 'bold');
% xlabel('Median OMZ upper boundary (m) ', 'FontSize', 14, 'FontWeight', 'bold')
% set(gca,'XLim', [300 700], 'YLim', [0 1200], 'YDir', 'reverse')
% fixfig
% 
% % by year
% figure; hold on % by unidivenum
% plot(OMZmedian_q_nY(:,1),Sq_median_q_nY(:,1), 'xr')
% plot(OMZmedian_q_nY(:,1),My_median_q_nY(:,1), 'ok')
% plot(OMZmedian_q_nY(:,1),Lo_median_q_nY(:,1), 'sc')
% plot(OMZmedian_q_nY(:,1),Ha_median_q_nY(:,1), '.b')
% ylabel('Median organism depth (m) ', 'FontSize', 14, 'FontWeight', 'bold');
% xlabel('Median OMZ upper boundary (m) ', 'FontSize', 14, 'FontWeight', 'bold')
% set(gca,'XLim', [300 700], 'YLim', [0 1200], 'YDir', 'reverse')
% fixfig
% 



%% 
  
disp('Completed sem_figsFinal.m') 
% ===== EOF [sem_figsFinal.m] ======  
