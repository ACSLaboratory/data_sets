%% Plots Data recorded from the OHC

clear all
close all
clc

saveFiles = 1; % Save the matfiles if 1.
%% 2 Robots

%2Robot ML OHC Data
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\2Robot\ML\Speed6_12
load OHCDataExperimentSpeed6_12.mat

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\2Robot\ML\Speed4_12
load OHCDataExperimentSpeed4_12.mat

%2Robot Stepping OHC Data
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\2Robot\Stepping
load OHCDataExperimentStepping2Bot.mat

%2Robot Stepping and ML Reference Data
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\2Robot\ML\Speed6_12
load RobotDataExperimentSpeed6_12.mat
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\2Robot\ML\Speed4_12
load RobotDataExperimentSpeed4_12.mat
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\2Robot\Stepping
load RobotDataExperimentStepping2Bot.mat

%% ML Speed 6 12

% Average and STD of all 10 runs
absVelLoad6_12(absVelLoad6_12(:,:)==0)=nan;

avgLoadVel = smooth(nanmean(absVelLoad6_12,2),20);
stdLoadVel = smooth(nanstd(absVelLoad6_12,0,2),20);

% Plot the data
cd E:\SeanResearch\PheenoChariot\Figures\2Robot\MachineLearning
endStep = 2840;
minSpeed = 6;

figure(1)% Plot Average 
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
shadedErrorBar(timeLoad6_12(1:endStep,2),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
hold on
ylim([0 8])
xlim([20 100])
XL = xlim;
YL = ylim;
line([XL(1) XL(2)],[minSpeed,minSpeed],'Color','k','LineWidth',4)
plot(timeLoad6_12(1:endStep,2),avgLoadVel(1:endStep),'-b','LineWidth',3.5)
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed6_12AvgVelocity','fig')
    print(gcf, 'ExperimentSpeed6_12AvgVelocity.pdf', '-dpdf','-r600');
    print(gcf, 'ExperimentSpeed6_12AvgVelocity.eps', '-depsc2','-r600');
    saveas(gcf,'ExperimentSpeed6_12AvgVelocity','tif')
end

colorMatrix = rand(size(absVelLoad6_12,2),3);
smoothN = 100;
figure(2)%Plot all Load Velocities with the Mean
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(timeLoad6_12(1:end-1,1),smooth(absVelLoad6_12(:,1),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(1,:))
hold on
plot(timeLoad6_12(1:end-1,2),smooth(absVelLoad6_12(:,2),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(2,:))
plot(timeLoad6_12(1:end-1,3),smooth(absVelLoad6_12(:,3),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(3,:))
plot(timeLoad6_12(1:end-1,4),smooth(absVelLoad6_12(:,4),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(4,:))
plot(timeLoad6_12(1:end-1,5),smooth(absVelLoad6_12(:,5),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(5,:))
plot(timeLoad6_12(1:end-1,6),smooth(absVelLoad6_12(:,6),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(6,:))
plot(timeLoad6_12(1:end-1,7),smooth(absVelLoad6_12(:,7),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(7,:))
plot(timeLoad6_12(1:end-1,8),smooth(absVelLoad6_12(:,8),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(8,:))
plot(timeLoad6_12(1:end-1,9),smooth(absVelLoad6_12(:,9),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(9,:))
plot(timeLoad6_12(1:end-1,10),smooth(absVelLoad6_12(:,10),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(10,:))
plot(timeLoad6_12(1:endStep,2),smooth(avgLoadVel(1:endStep),20),'LineWidth',6,'Color','r')
ylim([0 8])
xlim([20 100])
XL = xlim;
YL = ylim;
line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed6_12LoadVelocity','fig')
    print(gcf, 'ExperimentSpeed6_12LoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeed6_12LoadVelocity','tif')
end

% Individual Trials
cd E:\SeanResearch\PheenoChariot\Figures\2Robot\MachineLearning\IndividualTrials

% Plot Reference Velocities.
refVelBot306_12 = -refVelBot306_12;
refVelBot326_12 = -refVelBot326_12;

for i = 1:10
    RN = i;%ceil(unifrnd(1,10));
    saveFName1 = sprintf('ExperimentSpeed6_12Trial%dVelocity',RN);
    saveFName2 = sprintf('ExperimentSpeed6_12Trial%dReferenceVelocity',RN);

    figure(3)% Plot of an Individual Trial
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    %set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
    plot(timeBot326_12(1:endStep,RN),smooth(absVelBot326_12(1:endStep,RN),50),'LineWidth',2,'Color','r')
    hold on
    plot(timeLoad6_12(1:endStep,RN),smooth(absVelLoad6_12(1:endStep,RN),50),'LineWidth',2,'Color','b')
    plot(timeBot306_12(1:endStep,RN),smooth(absVelBot306_12(1:endStep,RN),50),'LineWidth',2,'Color','m')
    legend('Robot 1 Tracked','Load Tracked','Robot 2 Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 100])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName1,'fig')
        print(gcf,saveFName1, '-dpdf','-r600');
        saveas(gcf,saveFName1,'tif')
    end

    figure(6)%One Trial with Reference Velocities
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    plot(refTimeBot326_12(:,RN),refVelBot326_12(:,RN),'LineWidth',3,'Color','r')
    hold on
    plot(refTimeBot306_12(:,RN),refVelBot306_12(:,RN),'LineWidth',3,'Color','m')
    plot(timeLoad6_12(1:end-1,RN),smooth(absVelLoad6_12(:,RN),smoothN),'LineWidth',2,'Color','b')
    plot(timeBot306_12(1:end-1,RN),smooth(absVelBot306_12(:,RN),smoothN),'LineWidth',2,'Color','m','LineStyle','--')    
    plot(timeBot326_12(1:end-1,RN),smooth(absVelBot326_12(:,RN),smoothN),'LineWidth',2,'Color','r','LineStyle','--')
    legend('Robot One','Robot Two','Load','Location','SouthEast')
    ylim([0 8])
    xlim([0 100])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',4,'LineStyle','--')
    plot(timeLoad6_12(1:end-1,RN),smooth(absVelLoad6_12(:,RN),smoothN),'LineWidth',3.5,'Color','b')
    xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName2,'fig')
        print(gcf,saveFName2, '-dpdf','-r600');
        print(gcf, saveFName2, '-depsc2','-r600');
        saveas(gcf,saveFName2,'tif')
    end
end

%% ML Speed 4 12

% Average and STD of all 10 runs
absVelLoad4_12(absVelLoad4_12(:,:)==0 | absVelLoad4_12(:,:)>=6)=nan;

avgLoadVel = smooth(nanmean(absVelLoad4_12,2),20);
stdLoadVel = smooth(nanstd(absVelLoad4_12,0,2),20);

% Plot the data
cd E:\SeanResearch\PheenoChariot\Figures\2Robot\MachineLearning
endStep = 2840;
minSpeed2 = 4;

figure(1)% Plot Average 
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
shadedErrorBar(timeLoad4_12(1:endStep,2),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
hold on
ylim([0 8])
xlim([20 100])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4)
plot(timeLoad4_12(1:endStep,2),avgLoadVel(1:endStep),'-b','LineWidth',3.5) 
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed4_12AvgVelocity','fig')
    print(gcf, 'ExperimentSpeed4_12AvgVelocity.pdf', '-dpdf','-r600');
    print(gcf, 'ExperimentSpeed4_12AvgVelocity.eps', '-depsc2','-r600');
    saveas(gcf,'ExperimentSpeed4_12AvgVelocity','tif')
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
xlim([20 100])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed4_12LoadVelocity','fig')
    print(gcf, 'ExperimentSpeed4_12LoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeed4_12LoadVelocity','tif')
end

% Individual Trials
cd E:\SeanResearch\PheenoChariot\Figures\2Robot\MachineLearning\IndividualTrials

% Plot Reference Velocities.
refVelBot304_12 = -refVelBot304_12;
refVelBot324_12 = -refVelBot324_12;

for i = 1:10
    RN = i;%ceil(unifrnd(1,10));
    saveFName1 = sprintf('ExperimentSpeed4_12Trial%dVelocity',RN);
    saveFName2 = sprintf('ExperimentSpeed4_12Trial%dReferenceVelocity',RN);

    figure(3)% Plot of an Individual Trial
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    %set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
    plot(timeBot324_12(1:endStep,RN),smooth(absVelBot324_12(1:endStep,RN),50),'LineWidth',2,'Color','r')
    hold on
    plot(timeLoad4_12(1:endStep,RN),smooth(absVelLoad4_12(1:endStep,RN),50),'LineWidth',2,'Color','b')
    plot(timeBot304_12(1:endStep,RN),smooth(absVelBot304_12(1:endStep,RN),50),'LineWidth',2,'Color','m')
    legend('Robot 1 Tracked','Load Tracked','Robot 2 Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 100])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName1,'fig')
        print(gcf,saveFName1, '-dpdf','-r600');
        saveas(gcf,saveFName1,'tif')
    end

    figure(6)%One Trial with Reference Velocities
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    plot(refTimeBot324_12(:,RN),refVelBot324_12(:,RN),'LineWidth',3,'Color','r')
    hold on
    plot(refTimeBot304_12(:,RN),refVelBot304_12(:,RN),'LineWidth',3,'Color','m')
    plot(timeLoad4_12(1:end-1,RN),smooth(absVelLoad4_12(:,RN),smoothN),'LineWidth',2,'Color','b')
    plot(timeBot304_12(1:end-1,RN),smooth(absVelBot304_12(:,RN),smoothN),'LineWidth',2,'Color','m','LineStyle','--')
    plot(timeBot324_12(1:end-1,RN),smooth(absVelBot324_12(:,RN),smoothN),'LineWidth',2,'Color','r','LineStyle','--')
    legend('Robot One','Robot Two','Load','Location','SouthEast')
    ylim([0 8])
    xlim([0 100])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4,'LineStyle','--')
    plot(timeLoad4_12(1:end-1,RN),smooth(absVelLoad4_12(:,RN),smoothN),'LineWidth',3.5,'Color','b')
    xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName2,'fig')
        print(gcf,saveFName2, '-dpdf','-r600');
        print(gcf, saveFName2, '-depsc2','-r600');
        saveas(gcf,saveFName2,'tif')
    end
end

%% Stepping
steppingPrediction = (8-(0.5*91.8*0.045)/(1/2*1));

% Average and STD of all 10 runs
absVelLoadStepping2Bot(absVelLoadStepping2Bot(:,:)==0)=nan;

avgLoadVel = smooth(nanmean(absVelLoadStepping2Bot,2),20);
stdLoadVel = smooth(nanstd(absVelLoadStepping2Bot,0,2),20);

avgLoadVelStepping2Bot = nanmean(avgLoadVel(1:2100));
stdLoadVelStepping2Bot = nanstd(avgLoadVel(1:2100));

% Plot the data
cd E:\SeanResearch\PheenoChariot\Figures\2Robot\Stepping
endStep = 4000;

figure(1)% Plot Average 
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
%set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
shadedErrorBar(timeLoadStepping2Bot(1:endStep,7),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
ylim([0 8])
xlim([15 90])
XL = xlim;
YL = ylim;
line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
%legend('Robot 1','Load','Robot 2')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeedStepping2BotAvgVelocity','fig')
    print(gcf, 'ExperimentSpeedStepping2BotAvgVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeedStepping2BotAvgVelocity','tif')
end

colorMatrix = rand(size(absVelLoadStepping2Bot,2),3);
smoothN = 50;
figure(2)%Plot all Load Velocities with the Mean
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(timeLoadStepping2Bot(1:end-1,1),smooth(absVelLoadStepping2Bot(:,1),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(1,:))
hold on
plot(timeLoadStepping2Bot(1:end-1,2),smooth(absVelLoadStepping2Bot(:,2),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(2,:))
plot(timeLoadStepping2Bot(1:end-1,3),smooth(absVelLoadStepping2Bot(:,3),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(3,:))
plot(timeLoadStepping2Bot(1:end-1,4),smooth(absVelLoadStepping2Bot(:,4),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(4,:))
plot(timeLoadStepping2Bot(1:end-1,5),smooth(absVelLoadStepping2Bot(:,5),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(5,:))
plot(timeLoadStepping2Bot(1:end-1,6),smooth(absVelLoadStepping2Bot(:,6),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(6,:))
plot(timeLoadStepping2Bot(1:end-1,7),smooth(absVelLoadStepping2Bot(:,7),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(7,:))
plot(timeLoadStepping2Bot(1:end-1,8),smooth(absVelLoadStepping2Bot(:,8),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(8,:))
plot(timeLoadStepping2Bot(1:end-1,9),smooth(absVelLoadStepping2Bot(:,9),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(9,:))
plot(timeLoadStepping2Bot(1:end-1,10),smooth(absVelLoadStepping2Bot(:,10),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(10,:))
plot(timeLoadStepping2Bot(1:endStep,7),smooth(avgLoadVel(1:endStep),smoothN),'LineWidth',6,'Color','r')
ylim([0 8])
xlim([15 90])
XL = xlim;
YL = ylim;
line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeedStepping2BotLoadVelocity','fig')
    print(gcf, 'ExperimentSpeedStepping2BotLoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeedStepping2BotLoadVelocity','tif')
end

% Plot Reference Velocities.
refVelBot30Stepping2Bot = -refVelBot30Stepping2Bot;
refVelBot32Stepping2Bot = -refVelBot32Stepping2Bot;

% figure(4)%Plot Force of One Robot
% set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
% plot(refTimeBot30Stepping2Bot(:,1),forceBot30Stepping2Bot(:,1),'LineWidth',2,'Color',colorMatrix(1,:))
% hold on
% plot(refTimeBot30Stepping2Bot(:,2),forceBot30Stepping2Bot(:,2),'LineWidth',2,'Color',colorMatrix(2,:))
% plot(refTimeBot30Stepping2Bot(:,3),forceBot30Stepping2Bot(:,3),'LineWidth',2,'Color',colorMatrix(3,:))
% plot(refTimeBot30Stepping2Bot(:,4),forceBot30Stepping2Bot(:,4),'LineWidth',2,'Color',colorMatrix(4,:))
% plot(refTimeBot30Stepping2Bot(:,5),forceBot30Stepping2Bot(:,5),'LineWidth',2,'Color',colorMatrix(5,:))
% plot(refTimeBot30Stepping2Bot(:,6),forceBot30Stepping2Bot(:,6),'LineWidth',2,'Color',colorMatrix(6,:))
% plot(refTimeBot30Stepping2Bot(:,7),forceBot30Stepping2Bot(:,7),'LineWidth',2,'Color',colorMatrix(7,:))
% plot(refTimeBot30Stepping2Bot(:,8),forceBot30Stepping2Bot(:,8),'LineWidth',2,'Color',colorMatrix(8,:))
% plot(refTimeBot30Stepping2Bot(:,9),forceBot30Stepping2Bot(:,9),'LineWidth',2,'Color',colorMatrix(9,:))
% plot(refTimeBot30Stepping2Bot(:,10),forceBot30Stepping2Bot(:,10),'LineWidth',2,'Color',colorMatrix(10,:))
% ylim([0 700])
% xlim([0 90])
% XL = xlim;
% YL = ylim;
% xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
% ylabel('Force Applied(10-Bit)','fontsize',24,'FontWeight','Bold')
% set(gca,'fontsize',24,'FontWeight','Bold')
% grid on
% box on
% hold off
% if (saveFiles == 1)    
%     saveas(gcf,'ExperimentSpeedStepping2BotForceBot30','fig')
%     print(gcf, 'ExperimentSpeedStepping2BotForceBot30.pdf', '-dpdf','-r600');
%     saveas(gcf,'ExperimentSpeedStepping2BotForceBot30','tif')
% end
% 
% figure(5)%Plot Force of Other Robot
% set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
% %set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
% plot(refTimeBot32Stepping2Bot(:,1),forceBot32Stepping2Bot(:,1),'LineWidth',2,'Color',colorMatrix(1,:))
% hold on
% plot(refTimeBot32Stepping2Bot(:,2),forceBot32Stepping2Bot(:,2),'LineWidth',2,'Color',colorMatrix(2,:))
% plot(refTimeBot32Stepping2Bot(:,3),forceBot32Stepping2Bot(:,3),'LineWidth',2,'Color',colorMatrix(3,:))
% plot(refTimeBot32Stepping2Bot(:,4),forceBot32Stepping2Bot(:,4),'LineWidth',2,'Color',colorMatrix(4,:))
% plot(refTimeBot32Stepping2Bot(:,5),forceBot32Stepping2Bot(:,5),'LineWidth',2,'Color',colorMatrix(5,:))
% plot(refTimeBot32Stepping2Bot(:,6),forceBot32Stepping2Bot(:,6),'LineWidth',2,'Color',colorMatrix(6,:))
% plot(refTimeBot32Stepping2Bot(:,7),forceBot32Stepping2Bot(:,7),'LineWidth',2,'Color',colorMatrix(7,:))
% plot(refTimeBot32Stepping2Bot(:,8),forceBot32Stepping2Bot(:,8),'LineWidth',2,'Color',colorMatrix(8,:))
% plot(refTimeBot32Stepping2Bot(:,9),forceBot32Stepping2Bot(:,9),'LineWidth',2,'Color',colorMatrix(9,:))
% plot(refTimeBot32Stepping2Bot(:,10),forceBot32Stepping2Bot(:,10),'LineWidth',2,'Color',colorMatrix(10,:))
% ylim([0 700])
% xlim([0 90])
% XL = xlim;
% YL = ylim;
% xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
% ylabel('Force Applied(10-Bit)','fontsize',24,'FontWeight','Bold')
% set(gca,'fontsize',24,'FontWeight','Bold')
% grid on
% box on
% hold off
% if (saveFiles == 1)    
%     saveas(gcf,'ExperimentSpeedStepping2BotForceBot32','fig')
%     print(gcf, 'ExperimentSpeedStepping2BotForceBot32.pdf', '-dpdf','-r600');
%     saveas(gcf,'ExperimentSpeedStepping2BotForceBot32','tif')
% end

% Individual Trials
cd E:\SeanResearch\PheenoChariot\Figures\2Robot\Stepping\IndividualTrials

for i = 1:10
    RN = i;%ceil(unifrnd(1,10));
    saveFName1 = sprintf('ExperimentSpeedStepping2BotTrial%dVelocity',RN);
    saveFName2 = sprintf('ExperimentSpeedStepping2BotTrial%dReferenceVelocity',RN);

    figure(3)% Plot of an Individual Trial
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    %set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
    plot(timeBot32Stepping2Bot(1:endStep,RN),smooth(absVelBot32Stepping2Bot(1:endStep,RN),50),'LineWidth',2,'Color','r')
    hold on
    plot(timeLoadStepping2Bot(1:endStep,RN),smooth(absVelLoadStepping2Bot(1:endStep,RN),50),'LineWidth',2,'Color','b')
    plot(timeBot30Stepping2Bot(1:endStep,RN),smooth(absVelBot30Stepping2Bot(1:endStep,RN),50),'LineWidth',2,'Color','m')
    legend('Robot 1 Tracked','Load Tracked','Robot 2 Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 80])
    XL = xlim;
    YL = ylim;
    line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName1,'fig')
        print(gcf,saveFName1, '-dpdf','-r600');
        saveas(gcf,saveFName1,'tif')
    end

    figure(6)%One Trial with Reference Velocities
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    plot(refTimeBot32Stepping2Bot(:,RN),refVelBot32Stepping2Bot(:,RN),'LineWidth',3,'Color','r')
    hold on
    plot(timeBot32Stepping2Bot(1:end-1,RN),smooth(absVelBot32Stepping2Bot(:,RN),smoothN),'LineWidth',2,'Color','r','LineStyle','--')
    plot(timeLoadStepping2Bot(1:end-1,RN),smooth(absVelLoadStepping2Bot(:,RN),smoothN),'LineWidth',2,'Color','b')
    plot(refTimeBot30Stepping2Bot(:,RN),refVelBot30Stepping2Bot(:,RN),'LineWidth',3,'Color','m')
    plot(timeBot30Stepping2Bot(1:end-1,RN),smooth(absVelBot30Stepping2Bot(:,RN),smoothN),'LineWidth',2,'Color','m','LineStyle','--')
    legend('Robot One Reference','Robot One Tracked','Load Tracked','Robot Two Reference','Robot Two Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 80])
    XL = xlim;
    YL = ylim;
    line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName2,'fig')
        print(gcf,saveFName2, '-dpdf','-r600');
        saveas(gcf,saveFName2,'tif')
    end
end

%% 3 Robots

%3Robot ML OHC Data
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\3Robot\ML\Speed6_9_12
load OHCDataExperimentSpeed6_9_12.mat

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\3Robot\ML\Speed4_8_12
load OHCDataExperimentSpeed4_8_12.mat

%3Robot Stepping OHC Data
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\3Robot\Stepping
load OHCDataExperimentStepping3Bot.mat

%3Robot Stepping and ML Reference Data
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\3Robot\ML\Speed6_9_12
load RobotDataExperimentSpeed6_9_12.mat

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\3Robot\ML\Speed4_8_12
load RobotDataExperimentSpeed4_8_12.mat

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\3Robot\Stepping
load RobotDataExperimentStepping3Bot.mat

%% ML 6 9 12

% Average and STD of all 10 runs
absVelLoad6_9_12(absVelLoad6_9_12(:,:)==0|absVelLoad6_9_12(:,:)>=10)=nan;

avgLoadVel = smooth(nanmean(absVelLoad6_9_12,2),20);
stdLoadVel = smooth(nanstd(absVelLoad6_9_12,0,2),20);

% Plot the data
cd E:\SeanResearch\PheenoChariot\Figures\3Robot\MachineLearning
endStep = 2650;
minSpeed = 6;

figure(1)% Plot Average
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
%set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
shadedErrorBar(timeLoad6_9_12(1:endStep,5),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
ylim([0 8])
xlim([20 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',4)
hold on
plot(timeLoad6_9_12(1:endStep,5),avgLoadVel(1:endStep),'-b','LineWidth',3.5);
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
%legend('Robot 1','Load','Robot 2')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed6_9_12AvgVelocity','fig')
    print(gcf, 'ExperimentSpeed6_9_12AvgVelocity.pdf', '-dpdf','-r600');
    print(gcf, 'ExperimentSpeed6_9_12AvgVelocity.eps', '-depsc2','-r600');
    saveas(gcf,'ExperimentSpeed6_9_12AvgVelocity','tif')
end

colorMatrix = rand(size(absVelLoad6_9_12,2),3);
smoothN = 50;
figure(2)%Plot all Load Velocities with the Mean
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(timeLoad6_9_12(1:end-1,1),smooth(absVelLoad6_9_12(:,1),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(1,:))
hold on
plot(timeLoad6_9_12(1:end-1,2),smooth(absVelLoad6_9_12(:,2),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(2,:))
plot(timeLoad6_9_12(1:end-1,3),smooth(absVelLoad6_9_12(:,3),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(3,:))
plot(timeLoad6_9_12(1:end-1,4),smooth(absVelLoad6_9_12(:,4),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(4,:))
plot(timeLoad6_9_12(1:end-1,5),smooth(absVelLoad6_9_12(:,5),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(5,:))
plot(timeLoad6_9_12(1:end-1,6),smooth(absVelLoad6_9_12(:,6),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(6,:))
plot(timeLoad6_9_12(1:end-1,7),smooth(absVelLoad6_9_12(:,7),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(7,:))
plot(timeLoad6_9_12(1:end-1,8),smooth(absVelLoad6_9_12(:,8),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(8,:))
plot(timeLoad6_9_12(1:end-1,9),smooth(absVelLoad6_9_12(:,9),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(9,:))
plot(timeLoad6_9_12(1:end-1,10),smooth(absVelLoad6_9_12(:,10),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(10,:))
plot(timeLoad6_9_12(1:endStep,10),avgLoadVel(1:endStep),'LineWidth',6,'Color','r')
ylim([0 8])
xlim([15 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed6_9_12LoadVelocity','fig')
    print(gcf, 'ExperimentSpeed6_9_12LoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeed6_9_12LoadVelocity','tif')
end

% Plot Reference Velocities.
refVelBot306_9_12 = -refVelBot306_9_12;
refVelBot326_9_12 = -refVelBot326_9_12;
refVelBot36_9_12 = -refVelBot36_9_12;

% Individual Trials
cd E:\SeanResearch\PheenoChariot\Figures\3Robot\MachineLearning\IndividualTrials

for i = 1:10
    RN = i;%ceil(unifrnd(1,10));
    saveFName1 = sprintf('ExperimentSpeed6_9_12Trial%dVelocity',RN);
    saveFName2 = sprintf('ExperimentSpeed6_9_12Trial%dReferenceVelocity',RN);

    figure(3)% Plot of an Individual Trial
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    %set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
    plot(timeBot326_9_12(1:endStep,RN),smooth(absVelBot326_9_12(1:endStep,RN),50),'LineWidth',2,'Color','r')
    hold on
    plot(timeLoad6_9_12(1:endStep,RN),smooth(absVelLoad6_9_12(1:endStep,RN),50),'LineWidth',2,'Color','b')
    plot(timeBot306_9_12(1:endStep,RN),smooth(absVelBot306_9_12(1:endStep,RN),50),'LineWidth',2,'Color','m')
    plot(timeBot36_9_12(1:endStep,RN),smooth(absVelBot36_9_12(1:endStep,RN),50),'LineWidth',2,'Color','g')
    legend('Robot 1 Tracked','Load Tracked','Robot 2 Tracked','Robot 3 Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 90])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName1,'fig')
        print(gcf,saveFName1, '-dpdf','-r600');
        saveas(gcf,saveFName1,'tif')
    end

    figure(6)%One Trial with Reference Velocities
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    plot(refTimeBot326_9_12(:,RN),refVelBot326_9_12(:,RN),'LineWidth',3,'Color','r')
    hold on
    plot(refTimeBot306_9_12(:,RN),refVelBot306_9_12(:,RN),'LineWidth',3,'Color','m')
    plot(refTimeBot36_9_12(:,RN),refVelBot36_9_12(:,RN),'LineWidth',3,'Color','g')
    plot(timeLoad6_9_12(1:end-1,RN),smooth(absVelLoad6_9_12(:,RN),smoothN),'LineWidth',2,'Color','b')
    plot(timeBot326_9_12(1:end-1,RN),smooth(absVelBot326_9_12(:,RN),smoothN),'LineWidth',2,'Color','r','LineStyle','--')
    plot(timeBot306_9_12(1:end-1,RN),smooth(absVelBot306_9_12(:,RN),smoothN),'LineWidth',2,'Color','m','LineStyle','--')   
    plot(timeBot36_9_12(1:end-1,RN),smooth(absVelBot36_9_12(:,RN),smoothN),'LineWidth',2,'Color','g','LineStyle','--')
    legend('Robot One','Robot Two','Robot Three','Load','Location','SouthEast')
    ylim([0 8])
    xlim([0 90])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',4,'LineStyle','--')
    plot(timeLoad6_9_12(1:end-1,RN),smooth(absVelLoad6_9_12(:,RN),smoothN),'LineWidth',3.5,'Color','b')
    xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName2,'fig')
        print(gcf,saveFName2, '-dpdf','-r600');
        print(gcf, saveFName2, '-depsc2','-r600');
        saveas(gcf,saveFName2,'tif')
    end
end

%% ML 4 8 12

% Average and STD of all 10 runs
absVelLoad4_8_12(absVelLoad4_8_12(:,:)==0|absVelLoad4_8_12(:,:)>=10)=nan;

avgLoadVel = smooth(nanmean(absVelLoad4_8_12,2),20);
stdLoadVel = smooth(nanstd(absVelLoad4_8_12,0,2),20);

% Plot the data
cd E:\SeanResearch\PheenoChariot\Figures\3Robot\MachineLearning
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
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed4_8_12AvgVelocity','fig')
    print(gcf, 'ExperimentSpeed4_8_12AvgVelocity.pdf', '-dpdf','-r600');
    print(gcf, 'ExperimentSpeed4_8_12AvgVelocity.eps', '-depsc2','-r600');
    saveas(gcf,'ExperimentSpeed4_8_12AvgVelocity','tif')
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
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed4_8_12LoadVelocity','fig')
    print(gcf, 'ExperimentSpeed4_8_12LoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeed4_8_12LoadVelocity','tif')
end

% Plot Reference Velocities.
refVelBot304_8_12 = -refVelBot304_8_12;
refVelBot324_8_12 = -refVelBot324_8_12;
refVelBot34_8_12 = -refVelBot34_8_12;

% Individual Trials
cd E:\SeanResearch\PheenoChariot\Figures\3Robot\MachineLearning\IndividualTrials

for i = 1:10
    RN = i;%ceil(unifrnd(1,10));
    saveFName1 = sprintf('ExperimentSpeed4_8_12Trial%dVelocity',RN);
    saveFName2 = sprintf('ExperimentSpeed4_8_12Trial%dReferenceVelocity',RN);

    figure(3)% Plot of an Individual Trial
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    %set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
    plot(timeBot324_8_12(1:endStep,RN),smooth(absVelBot324_8_12(1:endStep,RN),50),'LineWidth',2,'Color','r')
    hold on
    plot(timeLoad4_8_12(1:endStep,RN),smooth(absVelLoad4_8_12(1:endStep,RN),50),'LineWidth',2,'Color','b')
    plot(timeBot304_8_12(1:endStep,RN),smooth(absVelBot304_8_12(1:endStep,RN),50),'LineWidth',2,'Color','m')
    plot(timeBot34_8_12(1:endStep,RN),smooth(absVelBot34_8_12(1:endStep,RN),50),'LineWidth',2,'Color','g')
    legend('Robot 1 Tracked','Load Tracked','Robot 2 Tracked','Robot 3 Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 90])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName1,'fig')
        print(gcf,saveFName1, '-dpdf','-r600');
        saveas(gcf,saveFName1,'tif')
    end

    figure(6)%One Trial with Reference Velocities
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    plot(refTimeBot324_8_12(:,RN),refVelBot324_8_12(:,RN),'LineWidth',3,'Color','r')
    hold on
    plot(refTimeBot304_8_12(:,RN),refVelBot304_8_12(:,RN),'LineWidth',3,'Color','m')
    plot(refTimeBot34_8_12(:,RN),refVelBot34_8_12(:,RN),'LineWidth',3,'Color','g')
    plot(timeLoad4_8_12(1:end-1,RN),smooth(absVelLoad4_8_12(:,RN),smoothN),'LineWidth',2,'Color','b')
    plot(timeBot324_8_12(1:end-1,RN),smooth(absVelBot324_8_12(:,RN),smoothN),'LineWidth',2,'Color','r','LineStyle','--')
    plot(timeBot304_8_12(1:end-1,RN),smooth(absVelBot304_8_12(:,RN),smoothN),'LineWidth',2,'Color','m','LineStyle','--')
    plot(timeBot34_8_12(1:end-1,RN),smooth(absVelBot34_8_12(:,RN),smoothN),'LineWidth',2,'Color','g','LineStyle','--')
    legend('Robot One','Robot Two','Robot Three','Load','Location','SouthEast')
    ylim([0 8])
    xlim([0 90])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4,'LineStyle','--')
    plot(timeLoad4_8_12(1:end-1,RN),smooth(absVelLoad4_8_12(:,RN),smoothN),'LineWidth',3.5,'Color','b')
    xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName2,'fig')
        print(gcf,saveFName2, '-dpdf','-r600');
        print(gcf, saveFName2, '-depsc2','-r600');
        saveas(gcf,saveFName2,'tif')
    end
end

%% Stepping

% Average and STD of all 10 runs
absVelLoadStepping3Bot(absVelLoadStepping3Bot(:,:)==0)=nan;

avgLoadVel = smooth(nanmean(absVelLoadStepping3Bot,2),20);
stdLoadVel = smooth(nanstd(absVelLoadStepping3Bot,0,2),20);

avgLoadVelStepping3Bot = nanmean(avgLoadVel(1:2100));
stdLoadVelStepping3Bot = nanstd(avgLoadVel(1:2100));

% Plot the data
cd E:\SeanResearch\PheenoChariot\Figures\3Robot\Stepping
endStep = 2500;

figure(1)% Plot Average 
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
%set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
shadedErrorBar(timeLoadStepping3Bot(1:endStep,1),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
ylim([0 8])
xlim([10 90])
XL = xlim;
YL = ylim;
line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
%legend('Robot 1','Load','Robot 2')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeedStepping3BotAvgVelocity','fig')
    print(gcf, 'ExperimentSpeedStepping3BotAvgVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeedStepping3BotAvgVelocity','tif')
end

colorMatrix = rand(size(absVelLoadStepping3Bot,2),3);
smoothN = 50;
figure(2)%Plot all Load Velocities with the Mean
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(timeLoadStepping3Bot(1:end-1,1),smooth(absVelLoadStepping3Bot(:,1),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(1,:))
hold on
plot(timeLoadStepping3Bot(1:end-1,2),smooth(absVelLoadStepping3Bot(:,2),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(2,:))
plot(timeLoadStepping3Bot(1:end-1,3),smooth(absVelLoadStepping3Bot(:,3),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(3,:))
plot(timeLoadStepping3Bot(1:end-1,4),smooth(absVelLoadStepping3Bot(:,4),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(4,:))
plot(timeLoadStepping3Bot(1:end-1,5),smooth(absVelLoadStepping3Bot(:,5),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(5,:))
plot(timeLoadStepping3Bot(1:end-1,6),smooth(absVelLoadStepping3Bot(:,6),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(6,:))
plot(timeLoadStepping3Bot(1:end-1,7),smooth(absVelLoadStepping3Bot(:,7),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(7,:))
plot(timeLoadStepping3Bot(1:end-1,8),smooth(absVelLoadStepping3Bot(:,8),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(8,:))
plot(timeLoadStepping3Bot(1:end-1,9),smooth(absVelLoadStepping3Bot(:,9),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(9,:))
plot(timeLoadStepping3Bot(1:end-1,10),smooth(absVelLoadStepping3Bot(:,10),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(10,:))
plot(timeLoadStepping3Bot(1:endStep,1),smooth(avgLoadVel(1:endStep),20),'LineWidth',6,'Color','r')
ylim([0 8])
xlim([10 90])
XL = xlim;
YL = ylim;
line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeedStepping3BotLoadVelocity','fig')
    print(gcf, 'ExperimentSpeedStepping3BotLoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeedStepping3BotLoadVelocity','tif')
end

% Plot Reference Velocities.
refVelBot3Stepping3Bot = -refVelBot3Stepping3Bot;
refVelBot30Stepping3Bot = -refVelBot30Stepping3Bot;
refVelBot32Stepping3Bot = -refVelBot32Stepping3Bot;

% Individual Trials
cd E:\SeanResearch\PheenoChariot\Figures\3Robot\Stepping\IndividualTrials

for i = 1:10
    RN = i;%ceil(unifrnd(1,10));
    saveFName1 = sprintf('ExperimentSpeedStepping3BotTrial%dVelocity',RN);
    saveFName2 = sprintf('ExperimentSpeedStepping3BotTrial%dReferenceVelocity',RN);

    figure(3)% Plot of an Individual Trial
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    %set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
    plot(timeBot32Stepping3Bot(1:endStep,RN),smooth(absVelBot32Stepping3Bot(1:endStep,RN),50),'LineWidth',2,'Color','r')
    hold on
    plot(timeLoadStepping3Bot(1:endStep,RN),smooth(absVelLoadStepping3Bot(1:endStep,RN),50),'LineWidth',2,'Color','b')
    plot(timeBot30Stepping3Bot(1:endStep,RN),smooth(absVelBot30Stepping3Bot(1:endStep,RN),50),'LineWidth',2,'Color','m')
    plot(timeBot3Stepping3Bot(1:endStep,RN),smooth(absVelBot3Stepping3Bot(1:endStep,RN),50),'LineWidth',2,'Color','g')
    legend('Robot 1 Tracked','Load Tracked','Robot 2 Tracked','Robot 3 Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 80])
    XL = xlim;
    YL = ylim;
    line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName1,'fig')
        print(gcf,saveFName1, '-dpdf','-r600');
        saveas(gcf,saveFName1,'tif')
    end

    figure(6)%One Trial with Reference Velocities
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    plot(refTimeBot32Stepping3Bot(:,RN),refVelBot32Stepping3Bot(:,RN),'LineWidth',3,'Color','r')
    hold on
    plot(timeBot32Stepping3Bot(1:end-1,RN),smooth(absVelBot32Stepping3Bot(:,RN),smoothN),'LineWidth',2,'Color','r','LineStyle','--')
    plot(timeLoadStepping3Bot(1:end-1,RN),smooth(absVelLoadStepping3Bot(:,RN),smoothN),'LineWidth',2,'Color','b')
    plot(refTimeBot30Stepping3Bot(:,RN),refVelBot30Stepping3Bot(:,RN),'LineWidth',3,'Color','m')
    plot(timeBot30Stepping3Bot(1:end-1,RN),smooth(absVelBot30Stepping3Bot(:,RN),smoothN),'LineWidth',2,'Color','m','LineStyle','--')
    plot(refTimeBot3Stepping3Bot(:,RN),refVelBot3Stepping3Bot(:,RN),'LineWidth',3,'Color','g')
    plot(timeBot3Stepping3Bot(1:end-1,RN),smooth(absVelBot3Stepping3Bot(:,RN),smoothN),'LineWidth',2,'Color','g','LineStyle','--')
    legend('Robot One Reference','Robot One Tracked','Load Tracked','Robot Two Reference','Robot Two Tracked','Robot Three Reference','Robot Three Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 80])
    XL = xlim;
    YL = ylim;
    line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName2,'fig')
        print(gcf,saveFName2, '-dpdf','-r600');
        saveas(gcf,saveFName2,'tif')
    end
end

%% 4 Robots

%4Robot ML OHC Data
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\4Robot\ML\Speed6_8_10_12
load OHCDataExperimentSpeed6_8_10_12.mat

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\4Robot\ML\Speed4_6_9_12
load OHCDataExperimentSpeed4_6_9_12.mat

%4Robot Stepping OHC Data
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\4Robot\Stepping
load OHCDataExperimentStepping4Bot.mat

%4Robot Stepping and ML Reference Data
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\4Robot\ML\Speed6_8_10_12
load RobotDataExperimentSpeed6_8_10_12.mat
cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\4Robot\ML\Speed4_6_9_12
load RobotDataExperimentSpeed4_6_9_12.mat

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\4Robot\Stepping
load RobotDataExperimentStepping4Bot.mat

%% ML

% Average and STD of all 10 runs
absVelLoad6_8_10_12(absVelLoad6_8_10_12(:,:)==0 | absVelLoad6_8_10_12(:,:) > 12)=nan;

avgLoadVel = smooth(nanmean(absVelLoad6_8_10_12,2),20);
stdLoadVel = smooth(nanstd(absVelLoad6_8_10_12,0,2),20);

% Plot the data
cd E:\SeanResearch\PheenoChariot\Figures\4Robot\MachineLearning
endStep = 2400;
minSpeed = 6;

figure(1)% Plot Average
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
%set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
shadedErrorBar(timeLoad6_8_10_12(1:endStep,3),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
ylim([0 8])
xlim([20 80])
XL = xlim;
YL = ylim;
line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',4)
hold on
plot(timeLoad6_8_10_12(1:endStep,3),avgLoadVel(1:endStep),'-b','LineWidth',3.5)
xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
%legend('Robot 1','Load','Robot 2')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed6_8_10_12AvgVelocity','fig')
    print(gcf, 'ExperimentSpeed6_8_10_12AvgVelocity.pdf', '-dpdf','-r600');
    print(gcf, 'ExperimentSpeed6_8_10_12AvgVelocity.eps', '-depsc2','-r600');
    saveas(gcf,'ExperimentSpeed6_8_10_12AvgVelocity','tif')
end

colorMatrix = rand(size(absVelLoad6_8_10_12,2),3);
smoothN = 100;
figure(2)%Plot all Load Velocities with the Mean
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(timeLoad6_8_10_12(1:end-1,1),smooth(absVelLoad6_8_10_12(:,1),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(1,:))
hold on
plot(timeLoad6_8_10_12(1:end-1,2),smooth(absVelLoad6_8_10_12(:,2),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(2,:))
plot(timeLoad6_8_10_12(1:end-1,3),smooth(absVelLoad6_8_10_12(:,3),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(3,:))
plot(timeLoad6_8_10_12(1:end-1,4),smooth(absVelLoad6_8_10_12(:,4),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(4,:))
plot(timeLoad6_8_10_12(1:end-1,5),smooth(absVelLoad6_8_10_12(:,5),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(5,:))
plot(timeLoad6_8_10_12(1:end-1,6),smooth(absVelLoad6_8_10_12(:,6),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(6,:))
plot(timeLoad6_8_10_12(1:end-1,7),smooth(absVelLoad6_8_10_12(:,7),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(7,:))
plot(timeLoad6_8_10_12(1:end-1,8),smooth(absVelLoad6_8_10_12(:,8),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(8,:))
plot(timeLoad6_8_10_12(1:end-1,9),smooth(absVelLoad6_8_10_12(:,9),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(9,:))
plot(timeLoad6_8_10_12(1:end-1,10),smooth(absVelLoad6_8_10_12(:,10),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(10,:))
plot(timeLoad6_8_10_12(1:endStep,8),smooth(avgLoadVel(1:endStep),20),'LineWidth',6,'Color','r')
ylim([0 8])
xlim([20 80])
XL = xlim;
YL = ylim;
line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed6_8_10_12LoadVelocity','fig')
    print(gcf, 'ExperimentSpeed6_8_10_12LoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeed6_8_10_12LoadVelocity','tif')
end

% Plot Reference Velocities.
refVelBot306_8_10_12 = -refVelBot306_8_10_12;
refVelBot326_8_10_12 = -refVelBot326_8_10_12;
refVelBot36_8_10_12 = -refVelBot36_8_10_12;
refVelBot466_8_10_12 = -refVelBot466_8_10_12;

% Individual Trials
cd E:\SeanResearch\PheenoChariot\Figures\4Robot\MachineLearning\IndividualTrials

for i = 1:10
    RN = i;%ceil(unifrnd(1,10));
    saveFName1 = sprintf('ExperimentSpeed6_8_10_12Trial%dVelocity',RN);
    saveFName2 = sprintf('ExperimentSpeed6_8_10_12Trial%dReferenceVelocity',RN);

    figure(3)% Plot of an Individual Trial
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    %set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
    plot(timeBot326_8_10_12(1:end-1,RN),smooth(absVelBot326_8_10_12(:,RN),50),'LineWidth',2,'Color','r')
    hold on
    plot(timeLoad6_8_10_12(1:end-1,RN),smooth(absVelLoad6_8_10_12(:,RN),50),'LineWidth',2,'Color','b')
    plot(timeBot306_8_10_12(1:end-1,RN),smooth(absVelBot306_8_10_12(:,RN),50),'LineWidth',2,'Color','m')
    plot(timeBot36_8_10_12(1:end-1,RN),smooth(absVelBot36_8_10_12(:,RN),50),'LineWidth',2,'Color','g')
    plot(timeBot466_8_10_12(1:end-1,RN),smooth(absVelBot466_8_10_12(:,RN),50),'LineWidth',2,'Color','c')
    legend('Robot 1 Tracked','Load Tracked','Robot 2 Tracked','Robot 3 Tracked','Robot 4 Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 80])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName1,'fig')
        print(gcf,saveFName1, '-dpdf','-r600');
        saveas(gcf,saveFName1,'tif')
    end

    figure(6)%One Trial with Reference Velocities
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    plot(refTimeBot326_8_10_12(:,RN),refVelBot326_8_10_12(:,RN),'LineWidth',3,'Color','r')
    hold on
    plot(refTimeBot306_8_10_12(:,RN),refVelBot306_8_10_12(:,RN),'LineWidth',3,'Color','m')
    plot(refTimeBot36_8_10_12(:,RN),refVelBot36_8_10_12(:,RN),'LineWidth',3,'Color','g')
    plot(refTimeBot466_8_10_12(:,RN),refVelBot466_8_10_12(:,RN),'LineWidth',3,'Color','c')
    plot(timeLoad6_8_10_12(1:end-1,RN),smooth(absVelLoad6_8_10_12(:,RN),smoothN),'LineWidth',2,'Color','b')
    plot(timeBot326_8_10_12(1:end-1,RN),smooth(absVelBot326_8_10_12(:,RN),smoothN),'LineWidth',2,'Color','r','LineStyle','--')
    plot(timeBot306_8_10_12(1:end-1,RN),smooth(absVelBot306_8_10_12(:,RN),smoothN),'LineWidth',2,'Color','m','LineStyle','--')
    plot(timeBot36_8_10_12(1:end-1,RN),smooth(absVelBot36_8_10_12(:,RN),smoothN),'LineWidth',2,'Color','g','LineStyle','--')
    plot(timeBot466_8_10_12(1:end-1,RN),smooth(absVelBot466_8_10_12(:,RN),smoothN),'LineWidth',2,'Color','c','LineStyle','--')
    legend('Robot One','Robot Two','Robot Three','Robot Four','Load Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 80])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed,minSpeed],'Color','k','LineWidth',4,'LineStyle','--')
    plot(timeLoad6_8_10_12(1:end-1,RN),smooth(absVelLoad6_8_10_12(:,RN),smoothN),'LineWidth',3.5,'Color','b')
    xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName2,'fig')
        print(gcf,saveFName2, '-dpdf','-r600');
        print(gcf, saveFName2, '-depsc2','-r600');
        saveas(gcf,saveFName2,'tif')
    end
end

%% ML

% Average and STD of all 10 runs
absVelLoad4_6_9_12(absVelLoad4_6_9_12(:,:)==0 | absVelLoad4_6_9_12(:,:) > 7)=nan;

avgLoadVel = smooth(nanmean(absVelLoad4_6_9_12,2),20);
stdLoadVel = smooth(nanstd(absVelLoad4_6_9_12,0,2),20);

% Plot the data
cd E:\SeanResearch\PheenoChariot\Figures\4Robot\MachineLearning
endStep = 2900;
minSpeed2 = 4;

figure(1)% Plot Average 
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
%set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
shadedErrorBar(timeLoad4_6_9_12(1:endStep,4),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
ylim([0 8])
xlim([20 90])
XL = xlim;
YL = ylim;
line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4)
hold on
plot(timeLoad4_6_9_12(1:endStep,4),avgLoadVel(1:endStep),'-b','LineWidth',3.5)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
%legend('Robot 1','Load','Robot 2')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed4_6_9_12AvgVelocity','fig')
    print(gcf, 'ExperimentSpeed4_6_9_12AvgVelocity.pdf', '-dpdf','-r600');
    print(gcf, 'ExperimentSpeed4_6_9_12AvgVelocity.eps', '-depsc2','-r600');
    saveas(gcf,'ExperimentSpeed4_6_9_12AvgVelocity','tif')
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
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeed4_6_9_12LoadVelocity','fig')
    print(gcf, 'ExperimentSpeed4_6_9_12LoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeed4_6_9_12LoadVelocity','tif')
end

% Plot Reference Velocities.
refVelBot304_6_9_12 = -refVelBot304_6_9_12;
refVelBot324_6_9_12 = -refVelBot324_6_9_12;
refVelBot34_6_9_12 = -refVelBot34_6_9_12;
refVelBot464_6_9_12 = -refVelBot464_6_9_12;

% Individual Trials
cd E:\SeanResearch\PheenoChariot\Figures\4Robot\MachineLearning\IndividualTrials

for i = 1:10
    RN = i;%ceil(unifrnd(1,10));
    saveFName1 = sprintf('ExperimentSpeed4_6_9_12Trial%dVelocity',RN);
    saveFName2 = sprintf('ExperimentSpeed4_6_9_12Trial%dReferenceVelocity',RN);

    figure(3)% Plot of an Individual Trial
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    %set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
    plot(timeBot324_6_9_12(1:end-1,RN),smooth(absVelBot324_6_9_12(:,RN),50),'LineWidth',2,'Color','r')
    hold on
    plot(timeLoad4_6_9_12(1:end-1,RN),smooth(absVelLoad4_6_9_12(:,RN),50),'LineWidth',2,'Color','b')
    plot(timeBot304_6_9_12(1:end-1,RN),smooth(absVelBot304_6_9_12(:,RN),50),'LineWidth',2,'Color','m')
    plot(timeBot34_6_9_12(1:end-1,RN),smooth(absVelBot34_6_9_12(:,RN),50),'LineWidth',2,'Color','g')
    plot(timeBot464_6_9_12(1:end-1,RN),smooth(absVelBot464_6_9_12(:,RN),50),'LineWidth',2,'Color','c')
    legend('Robot 1 Tracked','Load Tracked','Robot 2 Tracked','Robot 3 Tracked','Robot 4 Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 90])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName1,'fig')
        print(gcf,saveFName1, '-dpdf','-r600');
        saveas(gcf,saveFName1,'tif')
    end

    figure(6)%One Trial with Reference Velocities
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    plot(refTimeBot324_6_9_12(:,RN),refVelBot324_6_9_12(:,RN),'LineWidth',3,'Color','r')
    hold on
    plot(refTimeBot304_6_9_12(:,RN),refVelBot304_6_9_12(:,RN),'LineWidth',3,'Color','m')
    plot(refTimeBot34_6_9_12(:,RN),refVelBot34_6_9_12(:,RN),'LineWidth',3,'Color','g')
    plot(refTimeBot464_6_9_12(:,RN),refVelBot464_6_9_12(:,RN),'LineWidth',3,'Color','c')
    plot(timeLoad4_6_9_12(1:end-1,RN),smooth(absVelLoad4_6_9_12(:,RN),smoothN),'LineWidth',2,'Color','b')
    plot(timeBot324_6_9_12(1:end-1,RN),smooth(absVelBot324_6_9_12(:,RN),smoothN),'LineWidth',2,'Color','r','LineStyle','--')
    plot(timeBot304_6_9_12(1:end-1,RN),smooth(absVelBot304_6_9_12(:,RN),smoothN),'LineWidth',2,'Color','m','LineStyle','--')
    plot(timeBot34_6_9_12(1:end-1,RN),smooth(absVelBot34_6_9_12(:,RN),smoothN),'LineWidth',2,'Color','g','LineStyle','--')
    plot(timeBot464_6_9_12(1:end-1,RN),smooth(absVelBot464_6_9_12(:,RN),smoothN),'LineWidth',2,'Color','c','LineStyle','--')
    legend('Robot One','Robot Two','Robot Three','Robot Four','Load','Location','SouthEast')
    ylim([0 8])
    xlim([0 90])
    XL = xlim;
    YL = ylim;
    line(XL,[minSpeed2,minSpeed2],'Color','k','LineWidth',4,'LineStyle','--')
    plot(timeLoad4_6_9_12(1:end-1,RN),smooth(absVelLoad4_6_9_12(:,RN),smoothN),'LineWidth',3.5,'Color','b')
    xlabel('Time(s)','fontsize',30,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName2,'fig')
        print(gcf,saveFName2, '-dpdf','-r600');
        print(gcf, saveFName2, '-depsc2','-r600');
        saveas(gcf,saveFName2,'tif')
    end
end

%% Stepping

% Average and STD of all 10 runs
absVelLoadStepping4Bot(absVelLoadStepping4Bot(:,:)==0)=nan;

avgLoadVel = smooth(nanmean(absVelLoadStepping4Bot,2),20);
stdLoadVel = smooth(nanstd(absVelLoadStepping4Bot,0,2),20);

avgLoadVelStepping4Bot = nanmean(avgLoadVel(1:2100));
stdLoadVelStepping4Bot = nanstd(avgLoadVel(1:2100));

% Plot the data
cd E:\SeanResearch\PheenoChariot\Figures\4Robot\Stepping
endStep = 2800;

figure(1)% Plot Average 
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
%set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
shadedErrorBar(timeLoadStepping4Bot(1:endStep,6),avgLoadVel(1:endStep),stdLoadVel(1:endStep),{'-b','LineWidth',2});
ylim([0 8])
xlim([10 90])
XL = xlim;
YL = ylim;
line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
%legend('Robot 1','Load','Robot 2')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeedStepping4BotAvgVelocity','fig')
    print(gcf, 'ExperimentSpeedStepping4BotAvgVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeedStepping4BotAvgVelocity','tif')
end

colorMatrix = rand(size(absVelLoadStepping4Bot,2),3);
smoothN = 50;
figure(2)%Plot all Load Velocities with the Mean
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
plot(timeLoadStepping4Bot(1:end-1,1),smooth(absVelLoadStepping4Bot(:,1),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(1,:))
hold on
plot(timeLoadStepping4Bot(1:end-1,2),smooth(absVelLoadStepping4Bot(:,2),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(2,:))
plot(timeLoadStepping4Bot(1:end-1,3),smooth(absVelLoadStepping4Bot(:,3),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(3,:))
plot(timeLoadStepping4Bot(1:end-1,4),smooth(absVelLoadStepping4Bot(:,4),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(4,:))
plot(timeLoadStepping4Bot(1:end-1,5),smooth(absVelLoadStepping4Bot(:,5),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(5,:))
plot(timeLoadStepping4Bot(1:end-1,6),smooth(absVelLoadStepping4Bot(:,6),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(6,:))
plot(timeLoadStepping4Bot(1:end-1,7),smooth(absVelLoadStepping4Bot(:,7),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(7,:))
plot(timeLoadStepping4Bot(1:end-1,8),smooth(absVelLoadStepping4Bot(:,8),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(8,:))
plot(timeLoadStepping4Bot(1:end-1,9),smooth(absVelLoadStepping4Bot(:,9),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(9,:))
plot(timeLoadStepping4Bot(1:end-1,10),smooth(absVelLoadStepping4Bot(:,10),smoothN),'LineStyle','--','LineWidth',2,'Color',colorMatrix(10,:))
plot(timeLoadStepping4Bot(1:endStep,6),smooth(avgLoadVel(1:endStep),smoothN),'LineWidth',6,'Color','r')
ylim([0 8])
xlim([10 90])
XL = xlim;
YL = ylim;
line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeedStepping4BotLoadVelocity','fig')
    print(gcf, 'ExperimentSpeedStepping4BotLoadVelocity.pdf', '-dpdf','-r600');
    saveas(gcf,'ExperimentSpeedStepping4BotLoadVelocity','tif')
end

figure(4)% Plot Stepping Average
set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
%set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
h1 = errorbar([2],avgLoadVelStepping2Bot,stdLoadVelStepping2Bot,'.');
hold on
h2 = errorbar([3],avgLoadVelStepping3Bot,stdLoadVelStepping3Bot,'.');
h3 = errorbar([4],avgLoadVelStepping4Bot,stdLoadVelStepping4Bot,'.');
set(h1,'Linewidth',4,'MarkerSize',30)
set(h2,'Linewidth',4,'MarkerSize',30)
set(h3,'Linewidth',4,'MarkerSize',30)
ylim([2 6])
xlim([1.5 4.5])
XL = xlim;
YL = ylim;
line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',4)
xlabel('Team Size','fontsize',30,'FontWeight','Bold')
ylabel('Velocity(cm/s)','fontsize',30,'FontWeight','Bold')
set(gca,'fontsize',24,'FontWeight','Bold')
grid on
box on
hold off
if (saveFiles == 1)    
    saveas(gcf,'ExperimentSpeedStepping','fig')
    print(gcf, 'ExperimentSpeedStepping.pdf', '-dpdf','-r600');
    print(gcf, 'ExperimentSpeedStepping.eps', '-depsc2','-r600');
    saveas(gcf,'ExperimentSpeedStepping','tif')
end

% Plot Reference Velocities.
refVelBot3Stepping4Bot = -refVelBot3Stepping4Bot;
refVelBot30Stepping4Bot = -refVelBot30Stepping4Bot;
refVelBot32Stepping4Bot = -refVelBot32Stepping4Bot;
refVelBot46Stepping4Bot = -refVelBot46Stepping4Bot;

% Individual Trials
cd E:\SeanResearch\PheenoChariot\Figures\4Robot\Stepping\IndividualTrials

for i = 1:10
    RN = i;%ceil(unifrnd(1,10));
    saveFName1 = sprintf('ExperimentSpeedStepping4BotTrial%dVelocity',RN);
    saveFName2 = sprintf('ExperimentSpeedStepping4BotTrial%dReferenceVelocity',RN);

    figure(3)% Plot of an Individual Trial
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    %set(gcf,'Position', [0, 0, 1920, 1080],'PaperOrientation','landscape')
    plot(timeBot32Stepping4Bot(1:end-1,RN),smooth(absVelBot32Stepping4Bot(:,RN),50),'LineWidth',2,'Color','r')
    hold on
    plot(timeLoadStepping4Bot(1:end-1,RN),smooth(absVelLoadStepping4Bot(:,RN),50),'LineWidth',2,'Color','b')
    plot(timeBot30Stepping4Bot(1:end-1,RN),smooth(absVelBot30Stepping4Bot(:,RN),50),'LineWidth',2,'Color','m')
    plot(timeBot3Stepping4Bot(1:end-1,RN),smooth(absVelBot3Stepping4Bot(:,RN),50),'LineWidth',2,'Color','g')
    plot(timeBot46Stepping4Bot(1:end-1,RN),smooth(absVelBot46Stepping4Bot(:,RN),50),'LineWidth',2,'Color','c')
    legend('Robot 1 Tracked','Load Tracked','Robot 2 Tracked','Robot 3 Tracked','Robot 4 Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 80])
    XL = xlim;
    YL = ylim;
    line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName1,'fig')
        print(gcf,saveFName1, '-dpdf','-r600');
        saveas(gcf,saveFName1,'tif')
    end

    figure(6)%One Trial with Reference Velocities
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1])
    plot(refTimeBot32Stepping4Bot(:,RN),refVelBot32Stepping4Bot(:,RN),'LineWidth',3,'Color','r')
    hold on
    plot(timeBot32Stepping4Bot(1:end-1,RN),smooth(absVelBot32Stepping4Bot(:,RN),smoothN),'LineWidth',2,'Color','r','LineStyle','--')
    plot(timeLoadStepping4Bot(1:end-1,RN),smooth(absVelLoadStepping4Bot(:,RN),smoothN),'LineWidth',2,'Color','b')
    plot(refTimeBot30Stepping4Bot(:,RN),refVelBot30Stepping4Bot(:,RN),'LineWidth',3,'Color','m')
    plot(timeBot30Stepping4Bot(1:end-1,RN),smooth(absVelBot30Stepping4Bot(:,RN),smoothN),'LineWidth',2,'Color','m','LineStyle','--')
    plot(refTimeBot3Stepping4Bot(:,RN),refVelBot3Stepping4Bot(:,RN),'LineWidth',3,'Color','g')
    plot(timeBot3Stepping4Bot(1:end-1,RN),smooth(absVelBot3Stepping4Bot(:,RN),smoothN),'LineWidth',2,'Color','g','LineStyle','--')
    plot(refTimeBot46Stepping4Bot(:,RN),refVelBot46Stepping4Bot(:,RN),'LineWidth',3,'Color','c')
    plot(timeBot46Stepping4Bot(1:end-1,RN),smooth(absVelBot46Stepping4Bot(:,RN),smoothN),'LineWidth',2,'Color','c','LineStyle','--')
    legend('Robot One Reference','Robot One Tracked','Load Tracked','Robot Two Reference','Robot Two Tracked',...
        'Robot Three Reference','Robot Three Tracked','Robot Four Reference','Robot Four Tracked','Location','SouthEast')
    ylim([0 8])
    xlim([0 80])
    XL = xlim;
    YL = ylim;
    line(XL,[steppingPrediction,steppingPrediction],'Color','k','LineWidth',2)
    xlabel('Time(s)','fontsize',24,'FontWeight','Bold')
    ylabel('Velocity(cm/s)','fontsize',24,'FontWeight','Bold')
    set(gca,'fontsize',24,'FontWeight','Bold')
    grid on
    box on
    hold off
    if (saveFiles == 1)    
        saveas(gcf,saveFName2,'fig')
        print(gcf,saveFName2, '-dpdf','-r600');
        saveas(gcf,saveFName2,'tif')
    end
end