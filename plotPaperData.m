%% Plots Data recorded from the OHC

clear all
close all
clc

saveFiles = 1; % Save the matfiles if 1.
smoothNum = 150; %How many samples to smooth over (30 fps)
%% 2 Robots

%2Robot ML OHC Data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/OHCData/2Robot/ML/Speed4_12
load OHCDataExperimentSpeed4_12.mat

%2Robot Stepping OHC Data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/OHCData/2Robot/Stepping
load OHCDataExperimentStepping2Bot.mat

%2Robot Stepping and ML Reference Data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/RobotData/2Robot/ML/Speed4_12
load RobotDataExperimentSpeed4_12.mat

cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/RobotData/2Robot/Stepping
load RobotDataExperimentStepping2Bot.mat

%% ML Speed 4 12

% Average and STD of all 10 runs
absVelLoad4_12(absVelLoad4_12(:,:)==0 | absVelLoad4_12(:,:)>=6)=nan;

avgLoadVel = smooth(nanmean(absVelLoad4_12,2),20);
stdLoadVel = smooth(nanstd(absVelLoad4_12,0,2),20);

% Plot the data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/Figures/PaperFigs
endStep = 2840;
minSpeed2 = 4;

figure(1)% Plot Average 
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
shadedErrorBar(timeLoad4_12(1:endStep,2),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
hold on
ylim([0 8])
xlim([20 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4)
plot(timeLoad4_12(1:endStep,2),avgLoadVel(1:endStep),'-b','LineWidth',3.5) 
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'MLSpeed4_12AvgVelocity','fig')
    print(gcf, 'MLSpeed4_12AvgVelocity.pdf', '-dpdf','-r600');
    print(gcf, 'MLSpeed4_12AvgVelocity.eps', '-depsc2','-r600');
    saveas(gcf,'MLSpeed4_12AvgVelocity','tif')
end

colorMatrix = rand(size(absVelLoad4_12,2),3);
smoothN = 100;
figure(2)%Plot all Load Velocities with the Mean
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(timeLoad4_12(1:end-1,1),smooth(absVelLoad4_12(:,1),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(1,:))
hold on
plot(timeLoad4_12(1:end-1,2),smooth(absVelLoad4_12(:,2),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(2,:))
plot(timeLoad4_12(1:end-1,3),smooth(absVelLoad4_12(:,3),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(3,:))
plot(timeLoad4_12(1:end-1,4),smooth(absVelLoad4_12(:,4),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(4,:))
plot(timeLoad4_12(1:end-1,5),smooth(absVelLoad4_12(:,5),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(5,:))
plot(timeLoad4_12(1:end-1,6),smooth(absVelLoad4_12(:,6),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(6,:))
plot(timeLoad4_12(1:end-1,7),smooth(absVelLoad4_12(:,7),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(7,:))
plot(timeLoad4_12(1:end-1,8),smooth(absVelLoad4_12(:,8),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(8,:))
plot(timeLoad4_12(1:end-1,9),smooth(absVelLoad4_12(:,9),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(9,:))
plot(timeLoad4_12(1:end-1,10),smooth(absVelLoad4_12(:,10),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(10,:))
plot(timeLoad4_12(1:endStep,2),smooth(avgLoadVel(1:endStep),20),'LineWidth',6,'Color','r')
ylim([0 8])
xlim([20 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'MLSpeed4_12LoadVelocity','fig')
    print(gcf, 'MLSpeed4_12LoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'MLSpeed4_12LoadVelocity','tif')
end

% Plot Reference Velocities.
refVelBot304_12 = -refVelBot304_12;
refVelBot324_12 = -refVelBot324_12;

RN = 1;
saveFName1 = sprintf('MLSpeed4_12Trial%dTrackedVelocity',RN);
saveFName2 = sprintf('MLSpeed4_12Trial%dReferenceVelocity',RN);

figure(3)%Tracked Velocity of an Individual Trial
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
%set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
plot(timeBot324_12(1:endStep-smoothNum,RN),smooth(absVelBot324_12(1:endStep-smoothNum,RN),smoothNum),'LineWidth',3.5,'Color','r')
hold on
plot(timeBot304_12(1:endStep-smoothNum,RN),smooth(absVelBot304_12(1:endStep-smoothNum,RN),smoothNum),'LineWidth',3.5,'Color','m')
plot(timeLoad4_12(1:endStep,RN),smooth(absVelLoad4_12(1:endStep,RN),smoothNum),'LineWidth',4,'Color','b')
legend('Robot 1','Robot 2','Load','Location','SouthEast', 'AutoUpdate','off')
ylim([0 8])
xlim([0 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,saveFName1,'fig')
    print(gcf,saveFName1, '-dpdf','-r600');
    print(gcf, saveFName1, '-depsc2','-r600');
    saveas(gcf,saveFName1,'tif')
end

figure(6)%Reference Velocities of an Individual Trial
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(refTimeBot324_12(:,RN),refVelBot324_12(:,RN),'LineWidth',4.5,'Color','r')
hold on
plot(refTimeBot304_12(:,RN),refVelBot304_12(:,RN),'LineWidth',4.5,'Color','m')
legend('Robot 1','Robot 2','Location','SouthEast','AutoUpdate','off')
ylim([0 8])
xlim([0 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4,'LineStyle','--')
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,saveFName2,'fig')
    print(gcf,saveFName2, '-dpdf','-r600');
    print(gcf, saveFName2, '-depsc2','-r600');
    saveas(gcf,saveFName2,'tif')
end

%% Stepping 2 Bots
steppingPrediction = (8-(0.5*91.8*0.045)/(1/2*1));

% Average and STD of all 10 runs
absVelLoadStepping2Bot(absVelLoadStepping2Bot(:,:)==0)=nan;

avgLoadVel = smooth(nanmean(absVelLoadStepping2Bot,2),20);
stdLoadVel = smooth(nanstd(absVelLoadStepping2Bot,0,2),20);

avgLoadVelStepping2Bot = nanmean(avgLoadVel(1:2100));
stdLoadVelStepping2Bot = nanstd(avgLoadVel(1:2100));


%% 3 Robots

%3Robot ML OHC Data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/OHCData/3Robot/ML/Speed4_8_12 
load OHCDataExperimentSpeed4_8_12.mat

%3Robot Stepping OHC Data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/OHCData/3Robot/Stepping
load OHCDataExperimentStepping3Bot.mat

%3Robot Stepping and ML Reference Data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/RobotData/3Robot/ML/Speed4_8_12
load RobotDataExperimentSpeed4_8_12.mat

cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/RobotData/3Robot/Stepping
load RobotDataExperimentStepping3Bot.mat

%% ML 4 8 12

% Average and STD of all 10 runs
absVelLoad4_8_12(absVelLoad4_8_12(:,:)==0|absVelLoad4_8_12(:,:)>=10)=nan;

avgLoadVel = smooth(nanmean(absVelLoad4_8_12,2),20);
stdLoadVel = smooth(nanstd(absVelLoad4_8_12,0,2),20);

% Plot the data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/Figures/PaperFigs
endStep = 3130;
minSpeed2 = 4;

figure(1)% Plot Average
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
%set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
shadedErrorBar(timeLoad4_8_12(1:endStep,4),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
ylim([0 8])
xlim([20 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4)
hold on
plot(timeLoad4_8_12(1:endStep,4),avgLoadVel(1:endStep),'-b','LineWidth',3.5)
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
%legend('Robot 1','Load','Robot 2')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'MLSpeed4_8_12AvgVelocity','fig')
    print(gcf, 'MLSpeed4_8_12AvgVelocity.pdf', '-dpdf','-r600');
    print(gcf, 'MLSpeed4_8_12AvgVelocity.eps', '-depsc2','-r600');
    saveas(gcf,'MLSpeed4_8_12AvgVelocity','tif')
end

colorMatrix = rand(size(absVelLoad4_8_12,2),3);
smoothN = 50;
figure(2)%Plot all Load Velocities with the Mean
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(timeLoad4_8_12(1:end-1,1),smooth(absVelLoad4_8_12(:,1),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(1,:))
hold on
plot(timeLoad4_8_12(1:end-1,2),smooth(absVelLoad4_8_12(:,2),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(2,:))
plot(timeLoad4_8_12(1:end-1,3),smooth(absVelLoad4_8_12(:,3),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(3,:))
plot(timeLoad4_8_12(1:end-1,4),smooth(absVelLoad4_8_12(:,4),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(4,:))
plot(timeLoad4_8_12(1:end-1,5),smooth(absVelLoad4_8_12(:,5),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(5,:))
plot(timeLoad4_8_12(1:end-1,6),smooth(absVelLoad4_8_12(:,6),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(6,:))
plot(timeLoad4_8_12(1:end-1,7),smooth(absVelLoad4_8_12(:,7),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(7,:))
plot(timeLoad4_8_12(1:end-1,8),smooth(absVelLoad4_8_12(:,8),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(8,:))
plot(timeLoad4_8_12(1:end-1,9),smooth(absVelLoad4_8_12(:,9),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(9,:))
plot(timeLoad4_8_12(1:end-1,10),smooth(absVelLoad4_8_12(:,10),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(10,:))
plot(timeLoad4_8_12(1:endStep,10),avgLoadVel(1:endStep),'LineWidth',6,'Color','r')
ylim([0 8])
xlim([20 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'MLSpeed4_8_12LoadVelocity','fig')
    print(gcf, 'MLSpeed4_8_12LoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'MLSpeed4_8_12LoadVelocity','tif')
end

% Plot Reference Velocities.
refVelBot304_8_12 = -refVelBot304_8_12;
refVelBot324_8_12 = -refVelBot324_8_12;
refVelBot34_8_12 = -refVelBot34_8_12;

% Individual Trials

RN = 6;%ceil(unifrnd(1,10));
saveFName1 = sprintf('MLSpeed4_8_12Trial%dTrackedVelocity',RN);
saveFName2 = sprintf('MLSpeed4_8_12Trial%dReferenceVelocity',RN);

figure(3)% Plot of an Individual Trial
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
%set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
plot(timeBot324_8_12(50:endStep-smoothNum,RN),smooth(absVelBot324_8_12(50:endStep-smoothNum,RN),smoothNum),'LineWidth',3.5,'Color','r')
hold on
plot(timeBot304_8_12(50:endStep-smoothNum,RN),smooth(absVelBot304_8_12(50:endStep-smoothNum,RN),smoothNum),'LineWidth',3.5,'Color','m')
plot(timeBot34_8_12(10:endStep-smoothNum,RN),smooth(absVelBot34_8_12(10:endStep-smoothNum,RN),smoothNum),'LineWidth',3.5,'Color','g')
plot(timeLoad4_8_12(1:endStep,RN),smooth(absVelLoad4_8_12(1:endStep,RN),smoothNum),'LineWidth',4,'Color','b')
legend('Robot 1','Robot 2','Robot 3','Load','Location','SouthEast', 'AutoUpdate','off')
ylim([0 8])
xlim([0 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,saveFName1,'fig')
    print(gcf,saveFName1, '-dpdf','-r600');
    print(gcf, saveFName1, '-depsc2','-r600');
    saveas(gcf,saveFName1,'tif')
end

figure(6)%One Trial with Reference Velocities
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(refTimeBot324_8_12(:,RN),refVelBot324_8_12(:,RN),'LineWidth',4.5,'Color','r')
hold on
plot(refTimeBot304_8_12(:,RN),refVelBot304_8_12(:,RN),'LineWidth',4.5,'Color','m')
plot(refTimeBot34_8_12(:,RN),refVelBot34_8_12(:,RN),'LineWidth',4.5,'Color','g')
legend('Robot 1','Robot 2','Robot 3','Location','SouthEast', 'AutoUpdate','off')
ylim([0 8])
xlim([0 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4,'LineStyle','--')
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,saveFName2,'fig')
    print(gcf,saveFName2, '-dpdf','-r600');
    print(gcf, saveFName2, '-depsc2','-r600');
    saveas(gcf,saveFName2,'tif')
end

%% Stepping 3 Bots

% Average and STD of all 10 runs
absVelLoadStepping3Bot(absVelLoadStepping3Bot(:,:)==0)=nan;

avgLoadVel = smooth(nanmean(absVelLoadStepping3Bot,2),20);
stdLoadVel = smooth(nanstd(absVelLoadStepping3Bot,0,2),20);

avgLoadVelStepping3Bot = nanmean(avgLoadVel(1:2100));
stdLoadVelStepping3Bot = nanstd(avgLoadVel(1:2100));

%% 4 Robots

%4Robot ML OHC Data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/OHCData/4Robot/ML/Speed4_6_9_12
load OHCDataExperimentSpeed4_6_9_12.mat

%4Robot Stepping OHC Data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/OHCData/4Robot/Stepping
load OHCDataExperimentStepping4Bot.mat

%4Robot Stepping and ML Reference Data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/RobotData/4Robot/ML/Speed4_6_9_12
load RobotDataExperimentSpeed4_6_9_12.mat

cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/TrackingCode/DataForPaper/RobotData/4Robot/Stepping
load RobotDataExperimentStepping4Bot.mat

%% ML 4 6 9 12

% Average and STD of all 10 runs
absVelLoad4_6_9_12(absVelLoad4_6_9_12(:,:)==0 | absVelLoad4_6_9_12(:,:) > 7)=nan;

avgLoadVel = smooth(nanmean(absVelLoad4_6_9_12,2),20);
stdLoadVel = smooth(nanstd(absVelLoad4_6_9_12,0,2),20);

% Plot the data
cd /media/sean/Sean_RAID1/ASU/SeanResearch/PheenoChariot/Figures/PaperFigs
endStep = 2900;
minSpeed2 = 4;

figure(1)% Plot Average 
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
shadedErrorBar(timeLoad4_6_9_12(1:endStep,4),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
ylim([0 8])
xlim([20 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4)
hold on
plot(timeLoad4_6_9_12(1:endStep,4),avgLoadVel(1:endStep),'-b','LineWidth',3.5)
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
%legend('Robot 1','Load','Robot 2')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'MLSpeed4_6_9_12AvgVelocity','fig')
    print(gcf, 'MLSpeed4_6_9_12AvgVelocity.pdf', '-dpdf','-r600');
    print(gcf, 'MLSpeed4_6_9_12AvgVelocity.eps', '-depsc2','-r600');
    saveas(gcf,'MLSpeed4_6_9_12AvgVelocity','tif')
end

colorMatrix = rand(size(absVelLoad4_6_9_12,2),3);
smoothN = 100;
figure(2)%Plot all Load Velocities with the Mean
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(timeLoad4_6_9_12(1:end-1,1),smooth(absVelLoad4_6_9_12(:,1),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(1,:))
hold on
plot(timeLoad4_6_9_12(1:end-1,2),smooth(absVelLoad4_6_9_12(:,2),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(2,:))
plot(timeLoad4_6_9_12(1:end-1,3),smooth(absVelLoad4_6_9_12(:,3),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(3,:))
plot(timeLoad4_6_9_12(1:end-1,4),smooth(absVelLoad4_6_9_12(:,4),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(4,:))
plot(timeLoad4_6_9_12(1:end-1,5),smooth(absVelLoad4_6_9_12(:,5),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(5,:))
plot(timeLoad4_6_9_12(1:end-1,6),smooth(absVelLoad4_6_9_12(:,6),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(6,:))
plot(timeLoad4_6_9_12(1:end-1,7),smooth(absVelLoad4_6_9_12(:,7),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(7,:))
plot(timeLoad4_6_9_12(1:end-1,8),smooth(absVelLoad4_6_9_12(:,8),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(8,:))
plot(timeLoad4_6_9_12(1:end-1,9),smooth(absVelLoad4_6_9_12(:,9),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(9,:))
plot(timeLoad4_6_9_12(1:end-1,10),smooth(absVelLoad4_6_9_12(:,10),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(10,:))
plot(timeLoad4_6_9_12(1:endStep,4),smooth(avgLoadVel(1:endStep),20),'LineWidth',6,'Color','r')
ylim([0 8])
xlim([20 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'MLSpeed4_6_9_12LoadVelocity','fig')
    print(gcf, 'MLSpeed4_6_9_12LoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'MLSpeed4_6_9_12LoadVelocity','tif')
end

% Plot Reference Velocities.
refVelBot304_6_9_12 = -refVelBot304_6_9_12;
refVelBot324_6_9_12 = -refVelBot324_6_9_12;
refVelBot34_6_9_12 = -refVelBot34_6_9_12;
refVelBot464_6_9_12 = -refVelBot464_6_9_12;

% Individual Trials

RN = 9;%ceil(unifrnd(1,10));
saveFName1 = sprintf('MLSpeed4_6_9_12Trial%dTrackedVelocity',RN);
saveFName2 = sprintf('MLSpeed4_6_9_12Trial%dReferenceVelocity',RN);

figure(3)% Plot of an Individual Trial
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(timeBot324_6_9_12(50:end-1,RN),smooth(absVelBot324_6_9_12(50:end,RN),smoothNum),'LineWidth',3.5,'Color','r')
hold on
plot(timeBot304_6_9_12(100:end-1,RN),smooth(absVelBot304_6_9_12(100:end,RN),smoothNum),'LineWidth',3.5,'Color','m')
plot(timeBot34_6_9_12(50:end-1,RN),smooth(absVelBot34_6_9_12(50:end,RN),smoothNum),'LineWidth',3.5,'Color','g')
plot(timeBot464_6_9_12(100:end-1,RN),smooth(absVelBot464_6_9_12(100:end,RN),smoothNum),'LineWidth',3.5,'Color','c')
plot(timeLoad4_6_9_12(1:end-1,RN),smooth(absVelLoad4_6_9_12(:,RN),smoothNum),'LineWidth',4,'Color','b')
legend('Robot 1','Robot 2','Robot 3','Robot 4','Load','Location','SouthEast', 'AutoUpdate','off')
ylim([0 8])
xlim([0 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,saveFName1,'fig')
    print(gcf,saveFName1, '-dpdf','-r600');
    print(gcf, saveFName1, '-depsc2','-r600');
    saveas(gcf,saveFName1,'tif')
end

figure(6)%One Trial with Reference Velocities
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(refTimeBot324_6_9_12(:,RN),refVelBot324_6_9_12(:,RN),'LineWidth',4.5,'Color','r')
hold on
plot(refTimeBot304_6_9_12(:,RN),refVelBot304_6_9_12(:,RN),'LineWidth',4.5,'Color','m')
plot(refTimeBot34_6_9_12(:,RN),refVelBot34_6_9_12(:,RN),'LineWidth',4.5,'Color','g')
plot(refTimeBot464_6_9_12(:,RN),refVelBot464_6_9_12(:,RN),'LineWidth',4.5,'Color','c')
legend('Robot 1','Robot 2','Robot 3','Robot 4','Location','SouthEast', 'AutoUpdate','off')
ylim([0 8])
xlim([0 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4,'LineStyle','--')
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,saveFName2,'fig')
    print(gcf,saveFName2, '-dpdf','-r600');
    print(gcf, saveFName2, '-depsc2','-r600');
    saveas(gcf,saveFName2,'tif')
end

%% Stepping

% Average and STD of all 10 runs
absVelLoadStepping4Bot(absVelLoadStepping4Bot(:,:)==0)=nan;

avgLoadVel = smooth(nanmean(absVelLoadStepping4Bot,2),20);
stdLoadVel = smooth(nanstd(absVelLoadStepping4Bot,0,2),20);

avgLoadVelStepping4Bot = nanmean(avgLoadVel(1:2100));
stdLoadVelStepping4Bot = nanstd(avgLoadVel(1:2100));


figure(1)% Plot Stepping Average
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
h1 = errorbar([2],avgLoadVelStepping2Bot,stdLoadVelStepping2Bot,'.');
hold on
h2 = errorbar([3],avgLoadVelStepping3Bot,stdLoadVelStepping3Bot,'.');
h3 = errorbar([4],avgLoadVelStepping4Bot,stdLoadVelStepping4Bot,'.');
set(h1,'Linewidth',4,'MarkerSize',30,'Color','blue')
set(h2,'Linewidth',4,'MarkerSize',30,'Color','blue')
set(h3,'Linewidth',4,'MarkerSize',30,'Color','blue')
ylim([2 6])
xlim([1.5 4.5])
XL = xlim;
YL = ylim;
line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',4)
xlabel('Team Size','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',30,'FontWeight','Bold','XTick',[1 2 3 4])
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'SteppingAverage','fig')
    print(gcf,'SteppingAverage.pdf', '-dpdf','-r600');
    print(gcf,'SteppingAverage.eps', '-depsc2','-r600');
    saveas(gcf,'SteppingAverage','tif')
end
