%% A code loading a DAT file and apply the dPhi/dt fucntion to the data
clear all;
clc
%% Read in required data data from OMNI hourly set:
FID=fopen('omni2_all_years.dat','r');
FormatStr = ['%d %d %d %*s %*s %*s %*s %*s %*s %f ' ...
             '%*s %*s %*s %*s %*s %f %f %*s %*s %*s ' ...
             '%*s %*s %*s %*s %f %*s %*s %*s %f %*s ' ...
             '%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s ' ...
             '%f %*s %*s %*s %*s %*s %*s %*s %*s %*s ' ...
             '%*s %*s %*s %*s %*s'];
        
Data=textscan(FID,FormatStr);
fclose(FID);
% Extract data to vs_Rawariables:
Years_All = double(Data{1});
Days_All = double(Data{2});
Hours_All = double(Data{3});
BTs_All = Data{4};
Bys_All = Data{5};
Bzs_All = Data{6};
vs_All = Data{7};
ps_All = Data{8};
Dsts_All = Data{9};
% Clear the Data vs_Rawariable to free memory:
clearvars Data;
% Generate mask to select desired date range. As it's written below,it'll only work for wholeyears, not fractions of a year:
StartYear = 1984;
EndYear = 1994;
YearRangeMask = Years_All >= StartYear & ...
                Years_All <= EndYear;
 
% Extract Years_Raw range using the mask:
Years_Raw = Years_All(YearRangeMask);
Days_Raw = Days_All(YearRangeMask);
Hours_Raw = Hours_All(YearRangeMask);
BTs_Raw = BTs_All(YearRangeMask);
Bys_Raw = Bys_All(YearRangeMask);
Bzs_Raw = Bzs_All(YearRangeMask);
vs_Raw = vs_All(YearRangeMask);
ps_Raw = ps_All(YearRangeMask);
Dsts_Raw = Dsts_All(YearRangeMask);

% Clear unneeded vs_Rawariables to free memory:
 clearvars Years_All Days_All Hours_All BTs_All Bys_All ...
            Bzs_All ps_All vs_All Dsts_All...

theta_c = atan2(abs(Bys_Raw),Bzs_Raw); %the clock angle(radians)

%% Convs_Rawert to Julian day
JD = YearDayOfYearHourToMJD2000o0(Years_Raw, Days_Raw, Hours_Raw);

%% Calculate Coupling Function
dPhi_dt =  vs_Raw.^(4/3) .* BTs_Raw.^(2/3) .* (sin(theta_c./2)).^(8/3); 
dPhi_mod_BeforeIntegration = sqrt(ps_Raw).*dPhi_dt;

%% Masking -- marking the fault elements
MaskForGoodData = BTs_Raw<999.9 & ...
                  Bys_Raw<999.9 & ...
                  Bzs_Raw<999.9 & ...
                  vs_Raw<9999 & ...
                  ps_Raw<99.99;

ConvolvedMasks=conv(double(MaskForGoodData),ones(72,1),'valid');
MaskForGoodData72HrsBeforePadding = ConvolvedMasks > 71.5;
MaskForGoodData72Hrs = [false(71,1); MaskForGoodData72HrsBeforePadding];

%% 72-hour time integration 
n = 0:71;
weight_function = (0.95.^n);

dPhi_dt_AfterIntegrationBeforePadding = conv(dPhi_mod_BeforeIntegration, weight_function,'valid')/sum(weight_function);

dPhi_dt_AfterIntegration = [false(71,1); dPhi_dt_AfterIntegrationBeforePadding];

ps_AfterIntegrationBeforePadding = conv(sqrt(ps_Raw), weight_function, 'valid')/sum(weight_function);

ps_AfterIntegration = [false(71,1); ps_AfterIntegrationBeforePadding];

Dsts_Raw_mod = Dsts_Raw - 18.9*(ps_AfterIntegration); 

%% Pick out the FALSE values
dPhi_dt_AfterIntegration(MaskForGoodData72Hrs==0,:)=[];
Dsts_Raw_mod(MaskForGoodData72Hrs==0,:)=[];
JD(MaskForGoodData72Hrs==0,:)=[];

%% Correlation Coefficient
coef = corrcoef(dPhi_dt_AfterIntegration, Dsts_Raw_mod);

%% Straight line fit to Corrected Dst versus the coupling function
[Gradient,YIntercept] = BestFitLineFromPerpendicularOffsets(dPhi_dt_AfterIntegration, Dsts_Raw_mod);

%% Create plots
figure(1)
scatter(dPhi_dt_AfterIntegration,Dsts_Raw_mod,'k.');
xlabel({'sqrt(p)dPhi/dt'},'FontSize',12);
ylabel({'Dsts - 18.9sqrt(p)'},'FontSize',12);
title(['\fontsize{12} Year(s) ' num2str(Years_Raw(1,1)) ' to ' num2str(Years_Raw(end,1)) ' Correlation Coefficient ' num2str(coef(2,1))])
hold on
x=0:1:max(dPhi_dt_AfterIntegration);
plot(x, Gradient*x+YIntercept, 'r');
hold off
%%%%%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%




