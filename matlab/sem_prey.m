% sem_prey.m
%
% Called From: sem.prep

% Description: Need to do a lot of work in Text Wrangler for this.
% variable names often are labeled with shorthand: H=hake, L=loligo,
% M=myctophids
%
% Outside Functions Called:
%
% AUTHOR    : J. Stewart jules32@gmail.com
% DATE      : 16-Oct-2012 12:58:31
% Revision  : 1.00
% DEVELOPED : 7.9.0.529 (R2009b) OSX
% FILENAME  : sem_prey.m

Mfilename = mfilename;
dd = 0:50:4000; % depth ranges to work with

%% hake

prepneeded=0;
if prepneeded
    % hakeheaders.txt has added headers and None changed to NaN
    filename = '/Users/juliastewart/Dropbox/FATE_2012/DATA/ROVprey/hakeheaders.txt'; % change from hake/loligo/myctophidae
    
    fid = fopen(filename,'rt');
    temp = textscan(fid,'%f %f %s %f %f %f %f %f %f %s %s %s %f %f %s %s %f %s %f %s %s %s %s %s %s %f %s %s','delimiter','\t', 'headerlines', 1);
    fclose(fid);
    
    [PDiveNumber PDepth HConceptName PLatitude PLongitude PTemperature POxygen HX1 HX2 ...
        HRecordedDate HRovName HX3 HX4 VideoTime PAssociationsA PAssociationsB HX5 HShipName ...
        HX6 PAssociations HVideoArchiveName HChiefScientist HX7 HImage HX8 HNotes HObserver ...
        ] = temp{:}; % P stands for Prey
    
    clear HX1 HX2 HX3 HX4 HX5 HX6 HX7 HX8 HX9 PAssociationsA PAssociationsB
    PDepth(PDepth<0)=nan;
    
    % Dive Num
    DV = repmat(0.1,length(HRovName), 1); % so now any time it's a Ventana dive it has .1 after the DiveNumber,
    DV(strmatch('Tiburon', HRovName)) = 0.2; % and if it's a Tiburon it has .2
    DV(strmatch('Doc Ricketts', HRovName)) = 0.3; %and if it's a Doc Ricketts it has .3
    PDiveNumber = PDiveNumber + DV;
    
    HakeTrixA = [PDiveNumber PDepth PLatitude PLongitude PTemperature POxygen];
    
    % First match them to unidivenum
    
    HakeTrix = [];
    AssociationsS = [];
    HDateS = [];
    
    for i=1:length(unidivenum)
        PTemp = PDiveNumber == unidivenum(i);
        HakeTrix = [HakeTrix; HakeTrixA(PTemp,:)];
        AssociationsS = [AssociationsS; PAssociations(PTemp,:)];
        HDateS = [HDateS; HRecordedDate(PTemp,:)];
    end
    HakeTrix = [HakeTrix ones(size(HakeTrix,1),1)]; % this column of ones will turn into counts
    
    % not possible to deal with dates because there are weird ones with only time and some are >24h. So
    % assume that they are not seen on the downcast. 
    HDateC = regexp(HDateS, '\d{2}\:\d{2}', 'match');
%     HDateS(HDateC>24) = NaN;
    HDateN = datenum(HDateS);
    HDateV = datevec(HDateN);
    HDateV2 = datenum(0,0,0,HDateV(:,4),HDateV(:,5),HDateV(:,6));
    
    % after pull them all in and match with unidivenum and adding column for counts
    % put both HakeTrix and AssociationsS into excel saved as HakeFix.txt,
    % then text wrangler:
    % remove "maybe" and "incorrect"
    % go through and manually remove duplicate identity-refs
    % replace 999 with counts=5
    % replace school/schooling with counts=10
    % replace dense with counts=100 for myctophids, 50 for hake, 20 for market squid.
    % save as HakeFix2.txt
    % read HakeFix2.txt in and continue.
end

%% process cleaned prey data

Preyfiles = {'hake', 'loligo', 'myctophidae'};
PreyEvents = nan(length(unidivenum),length(Preyfiles));
PreyDSER = nan(length(unidivenum),length(Preyfiles));
PreyDepthMed = nan(length(unidivenum),length(Preyfiles));

% PreyDivesOnlydepthsDown = zeros(length(dd),length(unique(sqcountdivenumall))); % for depth
% PreyDivesOnlycountsDown = zeros(length(dd),length(unique(sqcountdivenumall)));

PreydepthbinscountDownH = nan(length(dd),length(unidivenum));
PreydepthbinscountDownL = nan(length(dd),length(unidivenum));
PreydepthbinscountDownM = nan(length(dd),length(unidivenum));

HaDivei = [];
LoDivei = [];
MyDivei = [];
HaDiveTrix = [];
LoDiveTrix = [];
MyDiveTrix = [];
PrcAll = nan(length(Preyfiles),2);

for p = 1:length(Preyfiles) % read in each prey file individually and process
    filename = ['/Users/juliastewart/Dropbox/FATE_2012/DATA/ROVprey/' Preyfiles{p} 'Fixed2.txt'];
    
    fid = fopen(filename,'rt');
    temp = textscan(fid,'%f %f %f %f %f %f %f %s','delimiter','\t', 'headerlines', 1);
    fclose(fid);
    
    [PDiveNumber PDepth PLat PLon PTemp POxy PCount PAssociations] = temp{:}; % P for Prey
    
    prc = [5 95];
    prc2 = prctile(PDepth, prc);
    PrcAll(p,:) = prc2;
    figure
    hist(PDepth, 1:25:2000)
    set(gca, 'XLim', [0 2000])
    xlabel('Depths (m)  ');
    ylabel('Counts ')
    title([Preyfiles(p); num2str(prc(2)-prc(1)) '% of depths between ' num2str(prc2(1)) ' and ' num2str(prc2(2))])
    fixfig
    FigName = ['DepthHistPrc_' Preyfiles{p} '.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
%     print('-dpdf', ['/Users/juliastewart/Dropbox/FATE_2012/workSEM/_figs/' FigName]);

    PEvents = zeros(length(unidivenum),1);
    PDSER = zeros(length(unidivenum),1);
    prcountdepthsall = [];
    PrDepth = [];
    
    pp = 0; % in order to index properly
    for i = 1:length(unidivenum)
        ii = PDiveNumber == unidivenum(i);
        oo = POxy < 10 & POxy > 0; % need to remove weird observations
        io = ii + oo;
        ilog = io == 2;
        ctd = divenum == unidivenum(i);
        idepth = PDepth(ilog);
        icount = PCount(ilog);
        ioxy = POxy(ilog);
        itemp = PTemp(ilog);
        dateind = dateCTD(ctd);
        depthCTD = depth(ctd);
        ymax = max(depthCTD);
        depthCTDmaxes(i) = ymax;
        ymaxindex = find(depthCTD == ymax);
        ymaxtime = dateind(ymaxindex(1)); % Take the first index encountered
        
        dateindV = datevec(dateind);
        dateindVSqYr = dateindV(:,1);
        YrInd = find(dateindV(:,2) < yrcutoff); 
        dateindVSqYr(YrInd) = dateindVSqYr(YrInd)-1;
        dateindVSqMonth = dateindV(:,2); % and change Jan-Feb to from 1,2 to 13,14 just to make it different from normal months
        MInd1 = find(dateindV(:,2) == 1);
        MInd2 = find(dateindV(:,2) == 2);
        dateindVSqMonth(MInd1) = 13;
        dateindVSqMonth(MInd2) = 14;

        if sum(icount) > 0
            prdepth = [];
            for u = 1:length(icount)
                stemp = repmat(idepth(u),icount(u),1);
                utemp = repmat(unidivenum(i),icount(u),1);
                ttemp = repmat(dateCTDstartSqYearV(i,1),icount(u),1);
                PrDepth = [PrDepth; stemp utemp ttemp];
            end
            prcountdepthsall = [prcountdepthsall; PrDepth];
            
            PreyDepthMed(i,p) = median(PrDepth(:,1)); % median prey depth per dive
        end
            
        if p == 1
            PreydepthbinscountDownH(:,i) = hist(idepth,dd); %(Sqtime<=ymaxtime)
            HaDivei = [PDiveNumber(ilog) idepth ioxy itemp icount ...
                repmat(dateindVSqYr(1),length(idepth),1) repmat(dateindVSqMonth(1),length(idepth),1)];
            hacountdepthsall = prcountdepthsall;
        elseif p == 2
            PreydepthbinscountDownL(:,i) = hist(idepth,dd); %(Sqtime<=ymaxtime)
            LoDivei = [PDiveNumber(ilog) idepth ioxy itemp icount ...
                repmat(dateindVSqYr(1),length(idepth),1) repmat(dateindVSqMonth(1),length(idepth),1)];
            locountdepthsall = prcountdepthsall;
        elseif p == 3
            PreydepthbinscountDownM(:,i) = hist(idepth,dd); %(Sqtime<=ymaxtime)
            MyDivei = [PDiveNumber(ilog) idepth ioxy itemp icount ...
                repmat(dateindVSqYr(1),length(idepth),1) repmat(dateindVSqMonth(1),length(idepth),1)];
            mycountdepthsall = prcountdepthsall;
        end
        
        if sum(PCount(ilog))>0
            PEvents(i) = sum(icount(idepth>prc2(1) & idepth<=prc2(2))); % total prey count for each dive, tailored for each prey depth range
            HDtemp = sum(PEvents(i))./sum(depthbinsdurDown(dd>prc2(1) & dd<prc2(2),i));
            
            % downcast
%             pp = pp+1;
%             PreyDivesOnlydepthsDown.ctd(:,pp) = hist(depthCTD(dateind<=ymaxtime),dd); %bindurDown;
%             PreyDivesOnlycountsDown.Preyfiles{p}(:,pp) = hist(idepth,dd); %(Sqtime<=ymaxtime)
            
            if isinf(HDtemp)
                PDSER(i) = NaN;
            else
                PDSER(i) = HDtemp;
            end
        end
        
        HaDiveTrix = [HaDiveTrix; HaDivei];
        LoDiveTrix = [LoDiveTrix; LoDivei];
        MyDiveTrix = [MyDiveTrix; MyDivei];
        
    end
    PreyEvents(:,p) = PEvents;
    PreyDSER(:,p) = PDSER;    
    
end

% Limit prey by their own depth ranges  
ADtotDownH = nansum(PreydepthbinscountDownH(dd>PrcAll(1,1) & dd<PrcAll(1,2),:),2)./nansum(depthbinsdurDown(dd>PrcAll(1,1) & dd<PrcAll(1,2),:),2);
ADtotDownL = nansum(PreydepthbinscountDownL(dd>PrcAll(2,1) & dd<PrcAll(2,2),:),2)./nansum(depthbinsdurDown(dd>PrcAll(2,1) & dd<PrcAll(2,2),:),2);
ADtotDownM = nansum(PreydepthbinscountDownM(dd>PrcAll(3,1) & dd<PrcAll(3,2),:),2)./nansum(depthbinsdurDown(dd>PrcAll(3,1) & dd<PrcAll(3,2),:),2);

% This was when we were limiting prey by Humboldt squid depth ranges  
% ADtotDownH = nansum(PreydepthbinscountDownH(sessbindef,:),2)./nansum(depthbinsdurDownSess,2);
% ADtotDownL = nansum(PreydepthbinscountDownL(sessbindef,:),2)./nansum(depthbinsdurDownSess,2);
% ADtotDownM = nansum(PreydepthbinscountDownM(sessbindef,:),2)./nansum(depthbinsdurDownSess,2);

% troubleshoot for 2010 loligos--although there are 30 unique dives from
% the text file, none show up in the count. 
% a = [3494; 3503; 3505; 3571; 181; 190; 3535; 157; 159; 3594; 3483; 3519; 3588; 214; 213; 215; 3493; 3600; 3545; 3599; 3477; 3504; 3502; 143; 152; 3580; 158; 3595; 212; 217];
% for i = 1:length(a)
%     aa = a(i) == unidivenum;
%     if sum(aa)>0
%         disp(a(i))
%     end
% end

save -mat MBARI_dataPrey_clean.mat HaDiveTrix LoDiveTrix MyDiveTrix PreyEvents PreyDSER ADtotDownH ADtotDownL ADtotDownM

% for Elliott/R
% cd(DirFATEsem)
% outfile = ['gamTrix_Hake_2013_02_22.txt'] % switch this and XXDiveTrix below to be what you want to save
% fout = fopen(outfile,'w');
% header = 'DiveNum Depth Oxy Temp Count SqYear SqMonth'; % write header to file
% fprintf(fout,'%s \n',header);
% fprintf(fout,'%4.1f\t%4.4f\t%2.3f\t%2.3f\t%4.0f\t%4.0f\t%2.0f\n', HaDiveTrix.');

% and for squid:
DiveDateVSqYr = DiveDateV(:,1);
YrInd = find(DiveDateV(:,2) < yrcutoff); 
DiveDateVSqYr(YrInd) = DiveDateVSqYr(YrInd)-1;
DiveDateVSqMonth = DiveDateV(:,2); % and change Jan-Feb to from 1,2 to 13,14 just to make it different from normal months
MInd1 = find(DiveDateV(:,2) == 1);
MInd2 = find(DiveDateV(:,2) == 2);
DiveDateVSqMonth(MInd1) = 13;
DiveDateVSqMonth(MInd2) = 14;
SqDiveTrix = [DiveNumber Depth Oxygen Temperature VarsCount DiveDateVSqYr DiveDateVSqMonth];

cd(DirFATEsem)
% outfile = ['gamTrix_Squid_2013_02_22.txt']
% fout = fopen(outfile,'w');
% header = 'DiveNum Depth Oxy Temp Count SqYear SqMonth'; % write header to file
% fprintf(fout,'%s \n',header);
% fprintf(fout,'%4.1f\t%4.4f\t%2.3f\t%2.3f\t%4.0f\t%4.0f\t%2.0f\n', SqDiveTrix.');

%% prey figures

figure
scatter(MyDiveTrix(:,3), MyDiveTrix(:,2), MyDiveTrix(:,5)*10) 
xlabel('oxygen (ml/l) ')
ylabel('sighting depth (m) ')
title('myctophid')
set(gca, 'YDir', 'reverse')

SqYMTot = zeros(12, length(SqYearUni));
HaYMTot = zeros(12, length(SqYearUni));
LoYMTot = zeros(12, length(SqYearUni));
MyYMTot = zeros(12, length(SqYearUni));

SqYMTotV = [];
HaYMTotV = [];
LoYMTotV = [];
MyYMTotV = [];
DateYMV = [];
DurYMV = [];
DurYMVprey = [];

for y = 1:length(SqYearUni)
    ylog = dateCTDstartSqYearV(:,1) == SqYearUni(y);
    for m = 1:length(SqMonthUni)
        mlog = SqMonthUni(m) == dateCTDstartSqYearV(:,2);
        my = ylog + mlog;
        mylog = my == 2;
        
        % matrices
        SqYMTot(m,y) = sum(DSER(mylog));
        HaYMTot(m,y) = sum(PreyDSER(mylog,1));
        LoYMTot(m,y) = sum(PreyDSER(mylog,2));
        MyYMTot(m,y) = sum(PreyDSER(mylog,3));
        
        % vectors
        DateYMV = [DateYMV; datenum(SqYearUni(y),SqMonthUni(m),1,0,0,0)];
        DurYMV = [DurYMV; sum(sum(depthbinsdurDownSess(:,mylog)))];
        DurYMVprey = [DurYMVprey; sum(sum(depthbinsdurDown(dd>PrcAll(1,1) & dd<PrcAll(1,2), mylog)))...
            sum(sum(depthbinsdurDown(dd>PrcAll(2,1) & dd<PrcAll(2,2), mylog)))...
            sum(sum(depthbinsdurDown(dd>PrcAll(3,1) & dd<PrcAll(3,2), mylog)))];
        
        SqYMTotV = [SqYMTotV; sum(DSER(mylog))];
        HaYMTotV = [HaYMTotV; sum(PreyDSER(mylog,1))];
        LoYMTotV = [LoYMTotV; sum(PreyDSER(mylog,2))];
        MyYMTotV = [MyYMTotV; sum(PreyDSER(mylog,3))];
        
    end
end
PreyYMTotV = [HaYMTotV LoYMTotV MyYMTotV]; % put into one matrix for plotting below

% timeseries squid: encounter rate and effort: Fig. 2?
figure
hold on
H(1) = bar(DateYMV, DurYMV, 'BarWidth',1);
set(H(1),'facecolor',[0.6 0.6 0.6],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1.5)
ylabel('ROV effort (hours) ', 'Color',[0.5 0.5 0.5], 'FontSize', 14, 'FontWeight', 'bold')
set(gca, 'XTick', DateYMV(1:12:end))
datetick('x', 'mm/yy', 'keepticks', 'keeplimits')
h1 = gca;
h2 = axes('Position',get(h1,'Position'));
H(2) = plot(DateYMV,log(SqYMTotV+1), 'k', 'LineWidth', 1.2);
set(h2,'YAxisLocation','right','Color','none','XTickLabel',[])
set(h2,'XLim',get(h1,'XLim'), 'YLim', [0 25], 'Layer','top')
xlabel('Date ', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Humboldt squid Encounter Rate (#/hour) ', 'Color','k', 'FontSize', 14, 'FontWeight', 'bold')
FigName = ['EncounterRate_Effort_Sq.pdf'];
annotate_JS(Mfilename, gcf, FigName)
% print('-dpdf', [DirFATEsem '_figs/' FigName]);


% prey timeseries figures with effort
for p = 1:length(Preyfiles)
    figure
    hold on
    H(1) = bar(DateYMV, DurYMVprey(:,p), 'BarWidth',1);
    set(H(1),'facecolor',[0.6 0.6 0.6],'EdgeColor',[0.5 0.5 0.5],'LineWidth',1.5)
    ylabel('ROV effort in prey depth swath (hours) ', 'Color',[0.5 0.5 0.5], 'FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'XTick', DateYMV(1:12:end))
    datetick('x', 'mm/yy', 'keepticks', 'keeplimits')
    h1 = gca;
    h2 = axes('Position',get(h1,'Position'));
    H(2) = plot(DateYMV,log(PreyYMTotV(:,p)+1), 'k', 'LineWidth', 1.2);
    set(h2,'YAxisLocation','right','Color','none','XTickLabel',[])
    set(h2,'XLim',get(h1,'XLim'), 'YLim', [0 25], 'Layer','top')
    xlabel('Date ', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel([Preyfiles{p} ' encounter rate (#/hour) '], 'Color','k', 'FontSize', 14, 'FontWeight', 'bold')
    FigName = ['EncounterRate_Effort_' Preyfiles{p} '.pdf'];
    annotate_JS(Mfilename, gcf, FigName)
%     print('-dpdf', [DirFATEsem '_figs/' FigName]);
end

% figure; hold on
% H(1) = bar(dd, sqplot, 'BarWidth',1);
% set(H(1),'facecolor',[0.5 0.5 0.5],'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.5)
% set(gca, 'XLim', [0 2500])
% ylabel('Raw count ', 'Color','k', 'FontSize', 14, 'FontWeight', 'bold')
% box off
% h1 = gca;
% h2 = axes('Position',get(h1,'Position'));
% H(2) = bar(dd, ADtotFull, 'BarWidth',1);
% hold off
% set(H(2),'facecolor','none','EdgeColor','k','LineWidth',1.5)
% set(h2,'YAxisLocation','right','Color','none','XTickLabel',[])
% set(h2,'XLim',get(h1,'XLim'),'Layer','top')
% xlabel('Depth (50 m bins) ', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('Adjusted count ', 'FontSize', 14, 'FontWeight', 'bold');


%%

disp('Completed sem_prey.m')
% ===== EOF [sem_prey.m] ======
