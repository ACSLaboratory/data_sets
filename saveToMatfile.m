%% Take OHC data and save it as mat files
%Sean Wilson
%ASU 2016

clear all
close all
clc

%% OHC Data for 2 Robots with ML Speed 6 12
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\2Robot\ML\Speed6_12

% Determine the largest file
maxSizeLoad = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:numTrials
    fname = sprintf('experiment%dSpeed6_12.csv',jj);
    data = csvread(fname,1,0);
    
    tempLoad = data(data(:,1)==40,:);
    maxSizeLoad = max(maxSizeLoad,length(tempLoad));
    
    tempBot30 = data(data(:,1)==30,:);
    maxSizeBot30 = max(maxSizeBot30,length(tempBot30));
    
    tempBot32 = data(data(:,1)==32,:);
    maxSizeBot32 = max(maxSizeBot32,length(tempBot32)); 
end

    xPosLoad6_12 = nan(maxSizeLoad,numTrials);yPosLoad6_12 = nan(maxSizeLoad,numTrials);...
        thetaLoad6_12 = nan(maxSizeLoad,numTrials); timeLoad6_12 = nan(maxSizeLoad,numTrials);
    xPosBot306_12 = nan(maxSizeBot30,numTrials);yPosBot306_12 = nan(maxSizeBot30,numTrials);...
        thetaBot306_12 =nan(maxSizeBot30,numTrials); timeBot306_12 = nan(maxSizeBot30,numTrials);
    xPosBot326_12 = nan(maxSizeBot32,numTrials);yPosBot326_12 = nan(maxSizeBot32,numTrials);...
        thetaBot326_12 = nan(maxSizeBot32,numTrials); timeBot326_12 = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials
    fname = sprintf('experiment%dSpeed6_12.csv',ii);
    data = csvread(fname,1,0);

    tempLoad = data(data(:,1)==40,:);
    tempBot30 = data(data(:,1)==30,:);
    tempBot32 = data(data(:,1)==32,:);

    %Separate Load Things
    xPosLoad6_12(1:length(tempLoad(:,2)),ii) = tempLoad(:,2);
    yPosLoad6_12(1:length(tempLoad(:,3)),ii) = tempLoad(:,3);
    thetaLoad6_12(1:length(tempLoad(:,4)),ii) = tempLoad(:,4);
    timeLoad6_12(1:length(tempLoad(:,5)),ii) = tempLoad(:,5);
    
    %Separate Bot30 Things
    xPosBot306_12(1:length(tempBot30(:,2)),ii) = tempBot30(:,2);
    yPosBot306_12(1:length(tempBot30(:,3)),ii) = tempBot30(:,3);
    thetaBot306_12(1:length(tempBot30(:,4)),ii) = tempBot30(:,4);
    timeBot306_12(1:length(tempBot30(:,5)),ii) = tempBot30(:,5);
    
    %Separate Bot32 Things
    xPosBot326_12(1:length(tempBot32(:,2)),ii) = tempBot32(:,2);
    yPosBot326_12(1:length(tempBot32(:,3)),ii) = tempBot32(:,3);
    thetaBot326_12(1:length(tempBot32(:,4)),ii) = tempBot32(:,4);
    timeBot326_12(1:length(tempBot32(:,5)),ii) = tempBot32(:,5);
end
    
% Determine Velocities

    velxLoad6_12 = (xPosLoad6_12(2:end,:)-xPosLoad6_12(1:end-1,:))./(timeLoad6_12(2:end,:)-timeLoad6_12(1:end-1,:));
    velyLoad6_12 = (yPosLoad6_12(2:end,:)-yPosLoad6_12(1:end-1,:))./(timeLoad6_12(2:end,:)-timeLoad6_12(1:end-1,:));
    velxLoad6_12(isnan(velxLoad6_12)) = 0;
    velyLoad6_12(isnan(velyLoad6_12)) = 0;
    
    absVelLoad6_12 = sqrt(velxLoad6_12.^2 + velyLoad6_12.^2);
    
    velxBot306_12 = (xPosBot306_12(2:end,:)-xPosBot306_12(1:end-1,:))./(timeBot306_12(2:end,:)-timeBot306_12(1:end-1,:));
    velyBot306_12 = (yPosBot306_12(2:end,:)-yPosBot306_12(1:end-1,:))./(timeBot306_12(2:end,:)-timeBot306_12(1:end-1,:));
    velxBot306_12(isnan(velxBot306_12)) = 0;
    velyBot306_12(isnan(velyBot306_12)) = 0;
    
    absVelBot306_12 = sqrt(velxBot306_12.^2 + velyBot306_12.^2);
    
    velxBot326_12 = (xPosBot326_12(2:end,:)-xPosBot326_12(1:end-1,:))./(timeBot326_12(2:end,:)-timeBot326_12(1:end-1,:));
    velyBot326_12 = (yPosBot326_12(2:end,:)-yPosBot326_12(1:end-1,:))./(timeBot326_12(2:end,:)-timeBot326_12(1:end-1,:));
    velxBot326_12(isnan(velxBot326_12)) = 0;
    velyBot326_12(isnan(velyBot326_12)) = 0;
    
    absVelBot326_12 = sqrt(velxBot326_12.^2 + velyBot326_12.^2);
    
if saveFiles == 1
    save('OHCDataExperimentSpeed6_12.mat','xPosLoad6_12','yPosLoad6_12',...
        'thetaLoad6_12','timeLoad6_12','xPosBot306_12','yPosBot306_12',...
        'thetaBot306_12','timeBot306_12','xPosBot326_12','yPosBot326_12',...
        'thetaBot326_12','timeBot326_12','absVelLoad6_12','absVelBot306_12',...
        'absVelBot326_12');
end

clear all
close all
clc
% Robot Data
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\2Robot\ML\Speed6_12
% Determine the largest file
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:10
    fname1 = sprintf('experiment%dBot30Speed6_12.csv',jj);
    fname2 = sprintf('experiment%dBot32Speed6_12.csv',jj);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    
    maxSizeBot30 = max(maxSizeBot30,length(data30));
    maxSizeBot32 = max(maxSizeBot32,length(data32));
end
    
    refVelBot306_12 = nan(maxSizeBot30,numTrials);
    refVelBot326_12 = nan(maxSizeBot32,numTrials);
    
    forceBot306_12 = nan(maxSizeBot30,numTrials);
    forceBot326_12 = nan(maxSizeBot32,numTrials);
    
    refTimeBot306_12 = nan(maxSizeBot30,numTrials);
    refTimeBot326_12 = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials

    fname1 = sprintf('experiment%dBot30Speed6_12.csv',ii);
    fname2 = sprintf('experiment%dBot32Speed6_12.csv',ii);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);

    %Separate Vel,Force,Time Things
    refVelBot306_12(1:length(data30),ii) = data30(:,1);
    refVelBot326_12(1:length(data32),ii) = data32(:,1);
    
    forceBot306_12(1:length(data30),ii) = data30(:,3);
    forceBot326_12(1:length(data32),ii) = data32(:,3);
    
    refTimeBot306_12(1:length(data30),ii) = data30(:,2)/1000;
    refTimeBot326_12(1:length(data32),ii) = data32(:,2)/1000;

end

if saveFiles == 1
    save('RobotDataExperimentSpeed6_12.mat','refVelBot306_12','refVelBot326_12',...
        'forceBot306_12','forceBot326_12','refTimeBot306_12','refTimeBot326_12');
end

clear all
close all
clc

%% OHC Data for 2 Robots with ML Speed 4 12
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\2Robot\ML\Speed4_12

% Determine the largest file
maxSizeLoad = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:numTrials
    fname = sprintf('experiment%dSpeed4_12.csv',jj);
    data = csvread(fname,1,0);
    
    tempLoad = data(data(:,1)==40,:);
    maxSizeLoad = max(maxSizeLoad,length(tempLoad));
    
    tempBot30 = data(data(:,1)==30,:);
    maxSizeBot30 = max(maxSizeBot30,length(tempBot30));
    
    tempBot32 = data(data(:,1)==32,:);
    maxSizeBot32 = max(maxSizeBot32,length(tempBot32)); 
end

    xPosLoad4_12 = nan(maxSizeLoad,numTrials);yPosLoad4_12 = nan(maxSizeLoad,numTrials);...
        thetaLoad4_12 = nan(maxSizeLoad,numTrials); timeLoad4_12 = nan(maxSizeLoad,numTrials);
    xPosBot304_12 = nan(maxSizeBot30,numTrials);yPosBot304_12 = nan(maxSizeBot30,numTrials);...
        thetaBot304_12 =nan(maxSizeBot30,numTrials); timeBot304_12 = nan(maxSizeBot30,numTrials);
    xPosBot324_12 = nan(maxSizeBot32,numTrials);yPosBot324_12 = nan(maxSizeBot32,numTrials);...
        thetaBot324_12 = nan(maxSizeBot32,numTrials); timeBot324_12 = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials
    fname = sprintf('experiment%dSpeed4_12.csv',ii);
    data = csvread(fname,1,0);

    tempLoad = data(data(:,1)==40,:);
    tempBot30 = data(data(:,1)==30,:);
    tempBot32 = data(data(:,1)==32,:);

    %Separate Load Things
    xPosLoad4_12(1:length(tempLoad(:,2)),ii) = tempLoad(:,2);
    yPosLoad4_12(1:length(tempLoad(:,3)),ii) = tempLoad(:,3);
    thetaLoad4_12(1:length(tempLoad(:,4)),ii) = tempLoad(:,4);
    timeLoad4_12(1:length(tempLoad(:,5)),ii) = tempLoad(:,5);
    
    %Separate Bot30 Things
    xPosBot304_12(1:length(tempBot30(:,2)),ii) = tempBot30(:,2);
    yPosBot304_12(1:length(tempBot30(:,3)),ii) = tempBot30(:,3);
    thetaBot304_12(1:length(tempBot30(:,4)),ii) = tempBot30(:,4);
    timeBot304_12(1:length(tempBot30(:,5)),ii) = tempBot30(:,5);
    
    %Separate Bot32 Things
    xPosBot324_12(1:length(tempBot32(:,2)),ii) = tempBot32(:,2);
    yPosBot324_12(1:length(tempBot32(:,3)),ii) = tempBot32(:,3);
    thetaBot324_12(1:length(tempBot32(:,4)),ii) = tempBot32(:,4);
    timeBot324_12(1:length(tempBot32(:,5)),ii) = tempBot32(:,5);
end
    
% Determine Velocities

    velxLoad4_12 = (xPosLoad4_12(2:end,:)-xPosLoad4_12(1:end-1,:))./(timeLoad4_12(2:end,:)-timeLoad4_12(1:end-1,:));
    velyLoad4_12 = (yPosLoad4_12(2:end,:)-yPosLoad4_12(1:end-1,:))./(timeLoad4_12(2:end,:)-timeLoad4_12(1:end-1,:));
    velxLoad4_12(isnan(velxLoad4_12)) = 0;
    velyLoad4_12(isnan(velyLoad4_12)) = 0;
    
    absVelLoad4_12 = sqrt(velxLoad4_12.^2 + velyLoad4_12.^2);
    
    velxBot304_12 = (xPosBot304_12(2:end,:)-xPosBot304_12(1:end-1,:))./(timeBot304_12(2:end,:)-timeBot304_12(1:end-1,:));
    velyBot304_12 = (yPosBot304_12(2:end,:)-yPosBot304_12(1:end-1,:))./(timeBot304_12(2:end,:)-timeBot304_12(1:end-1,:));
    velxBot304_12(isnan(velxBot304_12)) = 0;
    velyBot304_12(isnan(velyBot304_12)) = 0;
    
    absVelBot304_12 = sqrt(velxBot304_12.^2 + velyBot304_12.^2);
    
    velxBot324_12 = (xPosBot324_12(2:end,:)-xPosBot324_12(1:end-1,:))./(timeBot324_12(2:end,:)-timeBot324_12(1:end-1,:));
    velyBot324_12 = (yPosBot324_12(2:end,:)-yPosBot324_12(1:end-1,:))./(timeBot324_12(2:end,:)-timeBot324_12(1:end-1,:));
    velxBot324_12(isnan(velxBot324_12)) = 0;
    velyBot324_12(isnan(velyBot324_12)) = 0;
    
    absVelBot324_12 = sqrt(velxBot324_12.^2 + velyBot324_12.^2);
    
if saveFiles == 1
    save('OHCDataExperimentSpeed4_12.mat','xPosLoad4_12','yPosLoad4_12',...
        'thetaLoad4_12','timeLoad4_12','xPosBot304_12','yPosBot304_12',...
        'thetaBot304_12','timeBot304_12','xPosBot324_12','yPosBot324_12',...
        'thetaBot324_12','timeBot324_12','absVelLoad4_12','absVelBot304_12',...
        'absVelBot324_12');
end

clear all
close all
clc
% Robot Data
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\2Robot\ML\Speed4_12
% Determine the largest file
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:10
    fname1 = sprintf('experiment%dBot30Speed4_12.csv',jj);
    fname2 = sprintf('experiment%dBot32Speed4_12.csv',jj);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    
    maxSizeBot30 = max(maxSizeBot30,length(data30));
    maxSizeBot32 = max(maxSizeBot32,length(data32));
end
    
    refVelBot304_12 = nan(maxSizeBot30,numTrials);
    refVelBot324_12 = nan(maxSizeBot32,numTrials);
    
    forceBot304_12 = nan(maxSizeBot30,numTrials);
    forceBot324_12 = nan(maxSizeBot32,numTrials);
    
    refTimeBot304_12 = nan(maxSizeBot30,numTrials);
    refTimeBot324_12 = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials

    fname1 = sprintf('experiment%dBot30Speed4_12.csv',ii);
    fname2 = sprintf('experiment%dBot32Speed4_12.csv',ii);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);

    %Separate Vel,Force,Time Things
    refVelBot304_12(1:length(data30),ii) = data30(:,1);
    refVelBot324_12(1:length(data32),ii) = data32(:,1);
    
    forceBot304_12(1:length(data30),ii) = data30(:,3);
    forceBot324_12(1:length(data32),ii) = data32(:,3);
    
    refTimeBot304_12(1:length(data30),ii) = data30(:,2)/1000;
    refTimeBot324_12(1:length(data32),ii) = data32(:,2)/1000;

end

if saveFiles == 1
    save('RobotDataExperimentSpeed4_12.mat','refVelBot304_12','refVelBot324_12',...
        'forceBot304_12','forceBot324_12','refTimeBot304_12','refTimeBot324_12');
end

clear all
close all
clc


%% OHC Data for 2 Robots Stepping
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\2Robot\Stepping

% Determine the largest file
maxSizeLoad = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:10
    fname = sprintf('experiment%dStepping8.csv',jj);
    data = csvread(fname,1,0);
    
    tempLoad = data(data(:,1)==40,:);
    maxSizeLoad = max(maxSizeLoad,length(tempLoad));
    
    tempBot30 = data(data(:,1)==30,:);
    maxSizeBot30 = max(maxSizeBot30,length(tempBot30));
    
    tempBot32 = data(data(:,1)==32,:);
    maxSizeBot32 = max(maxSizeBot32,length(tempBot32)); 
end

    xPosLoadStepping2Bot = nan(maxSizeLoad,numTrials);yPosLoadStepping2Bot = nan(maxSizeLoad,numTrials);...
        thetaLoadStepping2Bot = nan(maxSizeLoad,numTrials); timeLoadStepping2Bot = nan(maxSizeLoad,numTrials);
    xPosBot30Stepping2Bot = nan(maxSizeBot30,numTrials);yPosBot30Stepping2Bot = nan(maxSizeBot30,numTrials);...
        thetaBot30Stepping2Bot =nan(maxSizeBot30,numTrials); timeBot30Stepping2Bot = nan(maxSizeBot30,numTrials);
    xPosBot32Stepping2Bot = nan(maxSizeBot32,numTrials);yPosBot32Stepping2Bot = nan(maxSizeBot32,numTrials);...
        thetaBot32Stepping2Bot = nan(maxSizeBot32,numTrials); timeBot32Stepping2Bot = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials
    fname = sprintf('experiment%dStepping8.csv',ii);
    data = csvread(fname,1,0);

    tempLoad = data(data(:,1)==40,:);
    tempBot30 = data(data(:,1)==30,:);
    tempBot32 = data(data(:,1)==32,:);

    %Separate Load Things
    xPosLoadStepping2Bot(1:length(tempLoad(:,2)),ii) = tempLoad(:,2);
    yPosLoadStepping2Bot(1:length(tempLoad(:,3)),ii) = tempLoad(:,3);
    thetaLoadStepping2Bot(1:length(tempLoad(:,4)),ii) = tempLoad(:,4);
    timeLoadStepping2Bot(1:length(tempLoad(:,5)),ii) = tempLoad(:,5);
    
    %Separate Bot30 Things
    xPosBot30Stepping2Bot(1:length(tempBot30(:,2)),ii) = tempBot30(:,2);
    yPosBot30Stepping2Bot(1:length(tempBot30(:,3)),ii) = tempBot30(:,3);
    thetaBot30Stepping2Bot(1:length(tempBot30(:,4)),ii) = tempBot30(:,4);
    timeBot30Stepping2Bot(1:length(tempBot30(:,5)),ii) = tempBot30(:,5);
    
    %Separate Bot32 Things
    xPosBot32Stepping2Bot(1:length(tempBot32(:,2)),ii) = tempBot32(:,2);
    yPosBot32Stepping2Bot(1:length(tempBot32(:,3)),ii) = tempBot32(:,3);
    thetaBot32Stepping2Bot(1:length(tempBot32(:,4)),ii) = tempBot32(:,4);
    timeBot32Stepping2Bot(1:length(tempBot32(:,5)),ii) = tempBot32(:,5);
end
    
% Determine Velocities

    velxLoadStepping2Bot = (xPosLoadStepping2Bot(2:end,:)-xPosLoadStepping2Bot(1:end-1,:))./(timeLoadStepping2Bot(2:end,:)-timeLoadStepping2Bot(1:end-1,:));
    velyLoadStepping2Bot = (yPosLoadStepping2Bot(2:end,:)-yPosLoadStepping2Bot(1:end-1,:))./(timeLoadStepping2Bot(2:end,:)-timeLoadStepping2Bot(1:end-1,:));
    velxLoadStepping2Bot(isnan(velxLoadStepping2Bot)) = 0;
    velyLoadStepping2Bot(isnan(velyLoadStepping2Bot)) = 0;
    
    absVelLoadStepping2Bot = sqrt(velxLoadStepping2Bot.^2 + velyLoadStepping2Bot.^2);
    
    velxBot30Stepping2Bot = (xPosBot30Stepping2Bot(2:end,:)-xPosBot30Stepping2Bot(1:end-1,:))./(timeBot30Stepping2Bot(2:end,:)-timeBot30Stepping2Bot(1:end-1,:));
    velyBot30Stepping2Bot = (yPosBot30Stepping2Bot(2:end,:)-yPosBot30Stepping2Bot(1:end-1,:))./(timeBot30Stepping2Bot(2:end,:)-timeBot30Stepping2Bot(1:end-1,:));
    velxBot30Stepping2Bot(isnan(velxBot30Stepping2Bot)) = 0;
    velyBot30Stepping2Bot(isnan(velyBot30Stepping2Bot)) = 0;
    
    absVelBot30Stepping2Bot = sqrt(velxBot30Stepping2Bot.^2 + velyBot30Stepping2Bot.^2);
    
    velxBot32Stepping2Bot = (xPosBot32Stepping2Bot(2:end,:)-xPosBot32Stepping2Bot(1:end-1,:))./(timeBot32Stepping2Bot(2:end,:)-timeBot32Stepping2Bot(1:end-1,:));
    velyBot32Stepping2Bot = (yPosBot32Stepping2Bot(2:end,:)-yPosBot32Stepping2Bot(1:end-1,:))./(timeBot32Stepping2Bot(2:end,:)-timeBot32Stepping2Bot(1:end-1,:));
    velxBot32Stepping2Bot(isnan(velxBot32Stepping2Bot)) = 0;
    velyBot32Stepping2Bot(isnan(velyBot32Stepping2Bot)) = 0;
    
    absVelBot32Stepping2Bot = sqrt(velxBot32Stepping2Bot.^2 + velyBot32Stepping2Bot.^2);
    
if saveFiles == 1
    save('OHCDataExperimentStepping2Bot.mat','xPosLoadStepping2Bot','yPosLoadStepping2Bot',...
        'thetaLoadStepping2Bot','timeLoadStepping2Bot','xPosBot30Stepping2Bot','yPosBot30Stepping2Bot',...
        'thetaBot30Stepping2Bot','timeBot30Stepping2Bot','xPosBot32Stepping2Bot','yPosBot32Stepping2Bot',...
        'thetaBot32Stepping2Bot','timeBot32Stepping2Bot','absVelLoadStepping2Bot','absVelBot30Stepping2Bot',...
        'absVelBot32Stepping2Bot');
end

clear all
close all
clc
% Robot Data
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\2Robot\Stepping
% Determine the largest file
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:10
    fname1 = sprintf('experiment%dBot30Stepping8.csv',jj);
    fname2 = sprintf('experiment%dBot32Stepping8.csv',jj);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    
    maxSizeBot30 = max(maxSizeBot30,length(data30));
    maxSizeBot32 = max(maxSizeBot32,length(data32));
end
    
    refVelBot30Stepping2Bot = nan(maxSizeBot30,numTrials);
    refVelBot32Stepping2Bot = nan(maxSizeBot32,numTrials);
    
    forceBot30Stepping2Bot = nan(maxSizeBot30,numTrials);
    forceBot32Stepping2Bot = nan(maxSizeBot32,numTrials);
    
    refTimeBot30Stepping2Bot = nan(maxSizeBot30,numTrials);
    refTimeBot32Stepping2Bot = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials

    fname1 = sprintf('experiment%dBot30Stepping8.csv',ii);
    fname2 = sprintf('experiment%dBot32Stepping8.csv',ii);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);

    %Separate Vel,Force,Time Things
    refVelBot30Stepping2Bot(1:length(data30),ii) = data30(:,1);
    refVelBot32Stepping2Bot(1:length(data32),ii) = data32(:,1);
    
    forceBot30Stepping2Bot(1:length(data30),ii) = data30(:,3);
    forceBot32Stepping2Bot(1:length(data32),ii) = data32(:,3);
    
    refTimeBot30Stepping2Bot(1:length(data30),ii) = data30(:,2)/1000;
    refTimeBot32Stepping2Bot(1:length(data32),ii) = data32(:,2)/1000;

end

if saveFiles == 1
    save('RobotDataExperimentStepping2Bot.mat','refVelBot30Stepping2Bot','refVelBot32Stepping2Bot',...
        'forceBot30Stepping2Bot','forceBot32Stepping2Bot','refTimeBot30Stepping2Bot','refTimeBot32Stepping2Bot');
end

%% OHC Data for 3 Robots with ML 6 9 12
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\3Robot\ML\Speed6_9_12

% Determine the largest file
maxSizeLoad = 0;
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:10
    fname = sprintf('experiment%dSpeed6_9_12.csv',jj);
    data = csvread(fname,1,0);
    
    tempLoad = data(data(:,1)==40,:);
    maxSizeLoad = max(maxSizeLoad,length(tempLoad));
    
    tempBot3 = data(data(:,1)==3,:);
    maxSizeBot3 = max(maxSizeBot3,length(tempBot3));
    
    tempBot30 = data(data(:,1)==30,:);
    maxSizeBot30 = max(maxSizeBot30,length(tempBot30));
    
    tempBot32 = data(data(:,1)==32,:);
    maxSizeBot32 = max(maxSizeBot32,length(tempBot32)); 
end

    xPosLoad6_9_12 = nan(maxSizeLoad,numTrials);yPosLoad6_9_12 = nan(maxSizeLoad,numTrials);...
        thetaLoad6_9_12 = nan(maxSizeLoad,numTrials); timeLoad6_9_12 = nan(maxSizeLoad,numTrials);
    xPosBot36_9_12 = nan(maxSizeBot3,numTrials);yPosBot36_9_12 = nan(maxSizeBot3,numTrials);...
        thetaBot36_9_12 = nan(maxSizeBot3,numTrials); timeBot36_9_12 = nan(maxSizeBot3,numTrials);
    xPosBot306_9_12 = nan(maxSizeBot30,numTrials);yPosBot306_9_12 = nan(maxSizeBot30,numTrials);...
        thetaBot306_9_12 =nan(maxSizeBot30,numTrials); timeBot306_9_12 = nan(maxSizeBot30,numTrials);
    xPosBot326_9_12 = nan(maxSizeBot32,numTrials);yPosBot326_9_12 = nan(maxSizeBot32,numTrials);...
        thetaBot326_9_12 = nan(maxSizeBot32,numTrials); timeBot326_9_12 = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials
    fname = sprintf('experiment%dSpeed6_9_12.csv',ii);
    data = csvread(fname,1,0);

    tempLoad = data(data(:,1)==40,:);
    tempBot3 = data(data(:,1)==3,:);
    tempBot30 = data(data(:,1)==30,:);
    tempBot32 = data(data(:,1)==32,:);

    %Separate Load Things
    xPosLoad6_9_12(1:length(tempLoad(:,2)),ii) = tempLoad(:,2);
    yPosLoad6_9_12(1:length(tempLoad(:,3)),ii) = tempLoad(:,3);
    thetaLoad6_9_12(1:length(tempLoad(:,4)),ii) = tempLoad(:,4);
    timeLoad6_9_12(1:length(tempLoad(:,5)),ii) = tempLoad(:,5);
    
    %Separate Bot3 Things
    xPosBot36_9_12(1:length(tempBot3(:,2)),ii) = tempBot3(:,2);
    yPosBot36_9_12(1:length(tempBot3(:,3)),ii) = tempBot3(:,3);
    thetaBot36_9_12(1:length(tempBot3(:,4)),ii) = tempBot3(:,4);
    timeBot36_9_12(1:length(tempBot3(:,5)),ii) = tempBot3(:,5);
    
    %Separate Bot30 Things
    xPosBot306_9_12(1:length(tempBot30(:,2)),ii) = tempBot30(:,2);
    yPosBot306_9_12(1:length(tempBot30(:,3)),ii) = tempBot30(:,3);
    thetaBot306_9_12(1:length(tempBot30(:,4)),ii) = tempBot30(:,4);
    timeBot306_9_12(1:length(tempBot30(:,5)),ii) = tempBot30(:,5);
    
    %Separate Bot32 Things
    xPosBot326_9_12(1:length(tempBot32(:,2)),ii) = tempBot32(:,2);
    yPosBot326_9_12(1:length(tempBot32(:,3)),ii) = tempBot32(:,3);
    thetaBot326_9_12(1:length(tempBot32(:,4)),ii) = tempBot32(:,4);
    timeBot326_9_12(1:length(tempBot32(:,5)),ii) = tempBot32(:,5);
end
    
% Determine Velocities

    velxLoad6_9_12 = (xPosLoad6_9_12(2:end,:)-xPosLoad6_9_12(1:end-1,:))./(timeLoad6_9_12(2:end,:)-timeLoad6_9_12(1:end-1,:));
    velyLoad6_9_12 = (yPosLoad6_9_12(2:end,:)-yPosLoad6_9_12(1:end-1,:))./(timeLoad6_9_12(2:end,:)-timeLoad6_9_12(1:end-1,:));
    velxLoad6_9_12(isnan(velxLoad6_9_12)) = 0;
    velyLoad6_9_12(isnan(velyLoad6_9_12)) = 0;
    
    absVelLoad6_9_12 = sqrt(velxLoad6_9_12.^2 + velyLoad6_9_12.^2);
    
    velxBot36_9_12 = (xPosBot36_9_12(2:end,:)-xPosBot36_9_12(1:end-1,:))./(timeBot36_9_12(2:end,:)-timeBot36_9_12(1:end-1,:));
    velyBot36_9_12 = (yPosBot36_9_12(2:end,:)-yPosBot36_9_12(1:end-1,:))./(timeBot36_9_12(2:end,:)-timeBot36_9_12(1:end-1,:));
    velxBot36_9_12(isnan(velxBot36_9_12)) = 0;
    velyBot36_9_12(isnan(velyBot36_9_12)) = 0;
    
    absVelBot36_9_12 = sqrt(velxBot36_9_12.^2 + velyBot36_9_12.^2);
    
    velxBot306_9_12 = (xPosBot306_9_12(2:end,:)-xPosBot306_9_12(1:end-1,:))./(timeBot306_9_12(2:end,:)-timeBot306_9_12(1:end-1,:));
    velyBot306_9_12 = (yPosBot306_9_12(2:end,:)-yPosBot306_9_12(1:end-1,:))./(timeBot306_9_12(2:end,:)-timeBot306_9_12(1:end-1,:));
    velxBot306_9_12(isnan(velxBot306_9_12)) = 0;
    velyBot306_9_12(isnan(velyBot306_9_12)) = 0;
    
    absVelBot306_9_12 = sqrt(velxBot306_9_12.^2 + velyBot306_9_12.^2);
    
    velxBot326_9_12 = (xPosBot326_9_12(2:end,:)-xPosBot326_9_12(1:end-1,:))./(timeBot326_9_12(2:end,:)-timeBot326_9_12(1:end-1,:));
    velyBot326_9_12 = (yPosBot326_9_12(2:end,:)-yPosBot326_9_12(1:end-1,:))./(timeBot326_9_12(2:end,:)-timeBot326_9_12(1:end-1,:));
    velxBot326_9_12(isnan(velxBot326_9_12)) = 0;
    velyBot326_9_12(isnan(velyBot326_9_12)) = 0;
    
    absVelBot326_9_12 = sqrt(velxBot326_9_12.^2 + velyBot326_9_12.^2);
    
if saveFiles == 1
    save('OHCDataExperimentSpeed6_9_12.mat','xPosLoad6_9_12','yPosLoad6_9_12',...
        'thetaLoad6_9_12','timeLoad6_9_12','xPosBot36_9_12','yPosBot36_9_12',...
        'thetaBot36_9_12','timeBot36_9_12','xPosBot306_9_12','yPosBot306_9_12',...
        'thetaBot306_9_12','timeBot306_9_12','xPosBot326_9_12','yPosBot326_9_12',...
        'thetaBot326_9_12','timeBot326_9_12','absVelLoad6_9_12','absVelBot306_9_12',...
        'absVelBot326_9_12','absVelBot36_9_12');
end

clear all
close all
clc
% Robot Data
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\3Robot\ML\Speed6_9_12
% Determine the largest file
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:10
    fname1 = sprintf('experiment%dBot30Speed6_9_12.csv',jj);
    fname2 = sprintf('experiment%dBot32Speed6_9_12.csv',jj);
    fname3 = sprintf('experiment%dBot03Speed6_9_12.csv',jj);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);
    
    maxSizeBot3 = max(maxSizeBot3,length(data3));
    maxSizeBot30 = max(maxSizeBot30,length(data30));
    maxSizeBot32 = max(maxSizeBot32,length(data32));
end
    
    refVelBot36_9_12 = nan(maxSizeBot3,numTrials);
    refVelBot306_9_12 = nan(maxSizeBot30,numTrials);
    refVelBot326_9_12 = nan(maxSizeBot32,numTrials);
    
    forceBot36_9_12 = nan(maxSizeBot3,numTrials);
    forceBot306_9_12 = nan(maxSizeBot30,numTrials);
    forceBot326_9_12 = nan(maxSizeBot32,numTrials);
    
    refTimeBot36_9_12 = nan(maxSizeBot3,numTrials);
    refTimeBot306_9_12 = nan(maxSizeBot30,numTrials);
    refTimeBot326_9_12 = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials

    fname1 = sprintf('experiment%dBot30Speed6_9_12.csv',ii);
    fname2 = sprintf('experiment%dBot32Speed6_9_12.csv',ii);
    fname3 = sprintf('experiment%dBot03Speed6_9_12.csv',ii);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);

    %Separate Vel,Force,Time Things
    refVelBot36_9_12(1:length(data3),ii) = data3(:,1) - 0.5;
    refVelBot306_9_12(1:length(data30),ii) = data30(:,1);
    refVelBot326_9_12(1:length(data32),ii) = data32(:,1);
    
    forceBot36_9_12(1:length(data3),ii) = data3(:,3);
    forceBot306_9_12(1:length(data30),ii) = data30(:,3);
    forceBot326_9_12(1:length(data32),ii) = data32(:,3);
    
    refTimeBot36_9_12(1:length(data3),ii) = data3(:,2)/1000;
    refTimeBot306_9_12(1:length(data30),ii) = data30(:,2)/1000;
    refTimeBot326_9_12(1:length(data32),ii) = data32(:,2)/1000;

end

if saveFiles == 1
    save('RobotDataExperimentSpeed6_9_12.mat','refVelBot306_9_12','refVelBot326_9_12',...
        'refVelBot36_9_12','forceBot306_9_12','forceBot326_9_12','forceBot36_9_12',...
        'refTimeBot306_9_12','refTimeBot326_9_12','refTimeBot36_9_12');
end

clear all
close all
clc

%% OHC Data for 3 Robots with ML 4 8 12
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\3Robot\ML\Speed4_8_12

% Determine the largest file
maxSizeLoad = 0;
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:10
    fname = sprintf('experiment%dSpeed4_8_12.csv',jj);
    data = csvread(fname,1,0);
    
    tempLoad = data(data(:,1)==40,:);
    maxSizeLoad = max(maxSizeLoad,length(tempLoad));
    
    tempBot3 = data(data(:,1)==3,:);
    maxSizeBot3 = max(maxSizeBot3,length(tempBot3));
    
    tempBot30 = data(data(:,1)==30,:);
    maxSizeBot30 = max(maxSizeBot30,length(tempBot30));
    
    tempBot32 = data(data(:,1)==32,:);
    maxSizeBot32 = max(maxSizeBot32,length(tempBot32)); 
end

    xPosLoad4_8_12 = nan(maxSizeLoad,numTrials);yPosLoad4_8_12 = nan(maxSizeLoad,numTrials);...
        thetaLoad4_8_12 = nan(maxSizeLoad,numTrials); timeLoad4_8_12 = nan(maxSizeLoad,numTrials);
    xPosBot34_8_12 = nan(maxSizeBot3,numTrials);yPosBot34_8_12 = nan(maxSizeBot3,numTrials);...
        thetaBot34_8_12 = nan(maxSizeBot3,numTrials); timeBot34_8_12 = nan(maxSizeBot3,numTrials);
    xPosBot304_8_12 = nan(maxSizeBot30,numTrials);yPosBot304_8_12 = nan(maxSizeBot30,numTrials);...
        thetaBot304_8_12 =nan(maxSizeBot30,numTrials); timeBot304_8_12 = nan(maxSizeBot30,numTrials);
    xPosBot324_8_12 = nan(maxSizeBot32,numTrials);yPosBot324_8_12 = nan(maxSizeBot32,numTrials);...
        thetaBot324_8_12 = nan(maxSizeBot32,numTrials); timeBot324_8_12 = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials
    fname = sprintf('experiment%dSpeed4_8_12.csv',ii);
    data = csvread(fname,1,0);

    tempLoad = data(data(:,1)==40,:);
    tempBot3 = data(data(:,1)==3,:);
    tempBot30 = data(data(:,1)==30,:);
    tempBot32 = data(data(:,1)==32,:);

    %Separate Load Things
    xPosLoad4_8_12(1:length(tempLoad(:,2)),ii) = tempLoad(:,2);
    yPosLoad4_8_12(1:length(tempLoad(:,3)),ii) = tempLoad(:,3);
    thetaLoad4_8_12(1:length(tempLoad(:,4)),ii) = tempLoad(:,4);
    timeLoad4_8_12(1:length(tempLoad(:,5)),ii) = tempLoad(:,5);
    
    %Separate Bot3 Things
    xPosBot34_8_12(1:length(tempBot3(:,2)),ii) = tempBot3(:,2);
    yPosBot34_8_12(1:length(tempBot3(:,3)),ii) = tempBot3(:,3);
    thetaBot34_8_12(1:length(tempBot3(:,4)),ii) = tempBot3(:,4);
    timeBot34_8_12(1:length(tempBot3(:,5)),ii) = tempBot3(:,5);
    
    %Separate Bot30 Things
    xPosBot304_8_12(1:length(tempBot30(:,2)),ii) = tempBot30(:,2);
    yPosBot304_8_12(1:length(tempBot30(:,3)),ii) = tempBot30(:,3);
    thetaBot304_8_12(1:length(tempBot30(:,4)),ii) = tempBot30(:,4);
    timeBot304_8_12(1:length(tempBot30(:,5)),ii) = tempBot30(:,5);
    
    %Separate Bot32 Things
    xPosBot324_8_12(1:length(tempBot32(:,2)),ii) = tempBot32(:,2);
    yPosBot324_8_12(1:length(tempBot32(:,3)),ii) = tempBot32(:,3);
    thetaBot324_8_12(1:length(tempBot32(:,4)),ii) = tempBot32(:,4);
    timeBot324_8_12(1:length(tempBot32(:,5)),ii) = tempBot32(:,5);
end
    
% Determine Velocities

    velxLoad4_8_12 = (xPosLoad4_8_12(2:end,:)-xPosLoad4_8_12(1:end-1,:))./(timeLoad4_8_12(2:end,:)-timeLoad4_8_12(1:end-1,:));
    velyLoad4_8_12 = (yPosLoad4_8_12(2:end,:)-yPosLoad4_8_12(1:end-1,:))./(timeLoad4_8_12(2:end,:)-timeLoad4_8_12(1:end-1,:));
    velxLoad4_8_12(isnan(velxLoad4_8_12)) = 0;
    velyLoad4_8_12(isnan(velyLoad4_8_12)) = 0;
    
    absVelLoad4_8_12 = sqrt(velxLoad4_8_12.^2 + velyLoad4_8_12.^2);
    
    velxBot34_8_12 = (xPosBot34_8_12(2:end,:)-xPosBot34_8_12(1:end-1,:))./(timeBot34_8_12(2:end,:)-timeBot34_8_12(1:end-1,:));
    velyBot34_8_12 = (yPosBot34_8_12(2:end,:)-yPosBot34_8_12(1:end-1,:))./(timeBot34_8_12(2:end,:)-timeBot34_8_12(1:end-1,:));
    velxBot34_8_12(isnan(velxBot34_8_12)) = 0;
    velyBot34_8_12(isnan(velyBot34_8_12)) = 0;
    
    absVelBot34_8_12 = sqrt(velxBot34_8_12.^2 + velyBot34_8_12.^2);
    
    velxBot304_8_12 = (xPosBot304_8_12(2:end,:)-xPosBot304_8_12(1:end-1,:))./(timeBot304_8_12(2:end,:)-timeBot304_8_12(1:end-1,:));
    velyBot304_8_12 = (yPosBot304_8_12(2:end,:)-yPosBot304_8_12(1:end-1,:))./(timeBot304_8_12(2:end,:)-timeBot304_8_12(1:end-1,:));
    velxBot304_8_12(isnan(velxBot304_8_12)) = 0;
    velyBot304_8_12(isnan(velyBot304_8_12)) = 0;
    
    absVelBot304_8_12 = sqrt(velxBot304_8_12.^2 + velyBot304_8_12.^2);
    
    velxBot324_8_12 = (xPosBot324_8_12(2:end,:)-xPosBot324_8_12(1:end-1,:))./(timeBot324_8_12(2:end,:)-timeBot324_8_12(1:end-1,:));
    velyBot324_8_12 = (yPosBot324_8_12(2:end,:)-yPosBot324_8_12(1:end-1,:))./(timeBot324_8_12(2:end,:)-timeBot324_8_12(1:end-1,:));
    velxBot324_8_12(isnan(velxBot324_8_12)) = 0;
    velyBot324_8_12(isnan(velyBot324_8_12)) = 0;
    
    absVelBot324_8_12 = sqrt(velxBot324_8_12.^2 + velyBot324_8_12.^2);
    
if saveFiles == 1
    save('OHCDataExperimentSpeed4_8_12.mat','xPosLoad4_8_12','yPosLoad4_8_12',...
        'thetaLoad4_8_12','timeLoad4_8_12','xPosBot34_8_12','yPosBot34_8_12',...
        'thetaBot34_8_12','timeBot34_8_12','xPosBot304_8_12','yPosBot304_8_12',...
        'thetaBot304_8_12','timeBot304_8_12','xPosBot324_8_12','yPosBot324_8_12',...
        'thetaBot324_8_12','timeBot324_8_12','absVelLoad4_8_12','absVelBot304_8_12',...
        'absVelBot324_8_12','absVelBot34_8_12');
end

clear all
close all
clc
% Robot Data
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\3Robot\ML\Speed4_8_12
% Determine the largest file
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:10
    fname1 = sprintf('experiment%dBot30Speed4_8_12.csv',jj);
    fname2 = sprintf('experiment%dBot32Speed4_8_12.csv',jj);
    fname3 = sprintf('experiment%dBot03Speed4_8_12.csv',jj);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);
    
    maxSizeBot3 = max(maxSizeBot3,length(data3));
    maxSizeBot30 = max(maxSizeBot30,length(data30));
    maxSizeBot32 = max(maxSizeBot32,length(data32));
end
    
    refVelBot34_8_12 = nan(maxSizeBot3,numTrials);
    refVelBot304_8_12 = nan(maxSizeBot30,numTrials);
    refVelBot324_8_12 = nan(maxSizeBot32,numTrials);
    
    forceBot34_8_12 = nan(maxSizeBot3,numTrials);
    forceBot304_8_12 = nan(maxSizeBot30,numTrials);
    forceBot324_8_12 = nan(maxSizeBot32,numTrials);
    
    refTimeBot34_8_12 = nan(maxSizeBot3,numTrials);
    refTimeBot304_8_12 = nan(maxSizeBot30,numTrials);
    refTimeBot324_8_12 = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials

    fname1 = sprintf('experiment%dBot30Speed4_8_12.csv',ii);
    fname2 = sprintf('experiment%dBot32Speed4_8_12.csv',ii);
    fname3 = sprintf('experiment%dBot03Speed4_8_12.csv',ii);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);

    %Separate Vel,Force,Time Things
    refVelBot34_8_12(1:length(data3),ii) = data3(:,1);
    refVelBot304_8_12(1:length(data30),ii) = data30(:,1);
    refVelBot324_8_12(1:length(data32),ii) = data32(:,1);
    
    forceBot34_8_12(1:length(data3),ii) = data3(:,3);
    forceBot304_8_12(1:length(data30),ii) = data30(:,3);
    forceBot324_8_12(1:length(data32),ii) = data32(:,3);
    
    refTimeBot34_8_12(1:length(data3),ii) = data3(:,2)/1000;
    refTimeBot304_8_12(1:length(data30),ii) = data30(:,2)/1000;
    refTimeBot324_8_12(1:length(data32),ii) = data32(:,2)/1000;

end

if saveFiles == 1
    save('RobotDataExperimentSpeed4_8_12.mat','refVelBot304_8_12','refVelBot324_8_12',...
        'refVelBot34_8_12','forceBot304_8_12','forceBot324_8_12','forceBot34_8_12',...
        'refTimeBot304_8_12','refTimeBot324_8_12','refTimeBot34_8_12');
end


clear all
close all
clc

%% OHC Data for 3 Robots Stepping
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\3Robot\Stepping

% Determine the largest file
maxSizeLoad = 0;
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:10
    fname = sprintf('experiment%dStepping8.csv',jj);
    data = csvread(fname,1,0);
    
    tempLoad = data(data(:,1)==40,:);
    maxSizeLoad = max(maxSizeLoad,length(tempLoad));
    
    tempBot3 = data(data(:,1)==3,:);
    maxSizeBot3 = max(maxSizeBot3,length(tempBot3));
    
    tempBot30 = data(data(:,1)==30,:);
    maxSizeBot30 = max(maxSizeBot30,length(tempBot30));
    
    tempBot32 = data(data(:,1)==32,:);
    maxSizeBot32 = max(maxSizeBot32,length(tempBot32)); 
end

    xPosLoadStepping3Bot = nan(maxSizeLoad,numTrials);yPosLoadStepping3Bot = nan(maxSizeLoad,numTrials);...
        thetaLoadStepping3Bot = nan(maxSizeLoad,numTrials); timeLoadStepping3Bot = nan(maxSizeLoad,numTrials);
    xPosBot3Stepping3Bot = nan(maxSizeBot3,numTrials);yPosBot3Stepping3Bot = nan(maxSizeBot3,numTrials);...
        thetaBot3Stepping3Bot =nan(maxSizeBot3,numTrials); timeBot3Stepping3Bot = nan(maxSizeBot3,numTrials);
    xPosBot30Stepping3Bot = nan(maxSizeBot30,numTrials);yPosBot30Stepping3Bot = nan(maxSizeBot30,numTrials);...
        thetaBot30Stepping3Bot =nan(maxSizeBot30,numTrials); timeBot30Stepping3Bot = nan(maxSizeBot30,numTrials);
    xPosBot32Stepping3Bot = nan(maxSizeBot32,numTrials);yPosBot32Stepping3Bot = nan(maxSizeBot32,numTrials);...
        thetaBot32Stepping3Bot = nan(maxSizeBot32,numTrials); timeBot32Stepping3Bot = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials
    fname = sprintf('experiment%dStepping8.csv',ii);
    data = csvread(fname,1,0);

    tempLoad = data(data(:,1)==40,:);
    tempBot3 = data(data(:,1)==3,:);
    tempBot30 = data(data(:,1)==30,:);
    tempBot32 = data(data(:,1)==32,:);

    %Separate Load Things
    xPosLoadStepping3Bot(1:length(tempLoad(:,2)),ii) = tempLoad(:,2);
    yPosLoadStepping3Bot(1:length(tempLoad(:,3)),ii) = tempLoad(:,3);
    thetaLoadStepping3Bot(1:length(tempLoad(:,4)),ii) = tempLoad(:,4);
    timeLoadStepping3Bot(1:length(tempLoad(:,5)),ii) = tempLoad(:,5);
    
    %Separate Bot3 Things
    xPosBot3Stepping3Bot(1:length(tempBot3(:,2)),ii) = tempBot3(:,2);
    yPosBot3Stepping3Bot(1:length(tempBot3(:,3)),ii) = tempBot3(:,3);
    thetaBot3Stepping3Bot(1:length(tempBot3(:,4)),ii) = tempBot3(:,4);
    timeBot3Stepping3Bot(1:length(tempBot3(:,5)),ii) = tempBot3(:,5);
    
    %Separate Bot30 Things
    xPosBot30Stepping3Bot(1:length(tempBot30(:,2)),ii) = tempBot30(:,2);
    yPosBot30Stepping3Bot(1:length(tempBot30(:,3)),ii) = tempBot30(:,3);
    thetaBot30Stepping3Bot(1:length(tempBot30(:,4)),ii) = tempBot30(:,4);
    timeBot30Stepping3Bot(1:length(tempBot30(:,5)),ii) = tempBot30(:,5);
    
    %Separate Bot32 Things
    xPosBot32Stepping3Bot(1:length(tempBot32(:,2)),ii) = tempBot32(:,2);
    yPosBot32Stepping3Bot(1:length(tempBot32(:,3)),ii) = tempBot32(:,3);
    thetaBot32Stepping3Bot(1:length(tempBot32(:,4)),ii) = tempBot32(:,4);
    timeBot32Stepping3Bot(1:length(tempBot32(:,5)),ii) = tempBot32(:,5);
end
    
% Determine Velocities

    velxLoadStepping3Bot = (xPosLoadStepping3Bot(2:end,:)-xPosLoadStepping3Bot(1:end-1,:))./(timeLoadStepping3Bot(2:end,:)-timeLoadStepping3Bot(1:end-1,:));
    velyLoadStepping3Bot = (yPosLoadStepping3Bot(2:end,:)-yPosLoadStepping3Bot(1:end-1,:))./(timeLoadStepping3Bot(2:end,:)-timeLoadStepping3Bot(1:end-1,:));
    velxLoadStepping3Bot(isnan(velxLoadStepping3Bot)) = 0;
    velyLoadStepping3Bot(isnan(velyLoadStepping3Bot)) = 0;
    
    absVelLoadStepping3Bot = sqrt(velxLoadStepping3Bot.^2 + velyLoadStepping3Bot.^2)-0.15;
    
    velxBot3Stepping3Bot = (xPosBot3Stepping3Bot(2:end,:)-xPosBot3Stepping3Bot(1:end-1,:))./(timeBot3Stepping3Bot(2:end,:)-timeBot3Stepping3Bot(1:end-1,:));
    velyBot3Stepping3Bot = (yPosBot3Stepping3Bot(2:end,:)-yPosBot3Stepping3Bot(1:end-1,:))./(timeBot3Stepping3Bot(2:end,:)-timeBot3Stepping3Bot(1:end-1,:));
    velxBot3Stepping3Bot(isnan(velxBot3Stepping3Bot)) = 0;
    velyBot3Stepping3Bot(isnan(velyBot3Stepping3Bot)) = 0;
    
    absVelBot3Stepping3Bot = sqrt(velxBot3Stepping3Bot.^2 + velyBot3Stepping3Bot.^2)-0.15;
    
    velxBot30Stepping3Bot = (xPosBot30Stepping3Bot(2:end,:)-xPosBot30Stepping3Bot(1:end-1,:))./(timeBot30Stepping3Bot(2:end,:)-timeBot30Stepping3Bot(1:end-1,:));
    velyBot30Stepping3Bot = (yPosBot30Stepping3Bot(2:end,:)-yPosBot30Stepping3Bot(1:end-1,:))./(timeBot30Stepping3Bot(2:end,:)-timeBot30Stepping3Bot(1:end-1,:));
    velxBot30Stepping3Bot(isnan(velxBot30Stepping3Bot)) = 0;
    velyBot30Stepping3Bot(isnan(velyBot30Stepping3Bot)) = 0;
    
    absVelBot30Stepping3Bot = sqrt(velxBot30Stepping3Bot.^2 + velyBot30Stepping3Bot.^2)-0.15;
    
    velxBot32Stepping3Bot = (xPosBot32Stepping3Bot(2:end,:)-xPosBot32Stepping3Bot(1:end-1,:))./(timeBot32Stepping3Bot(2:end,:)-timeBot32Stepping3Bot(1:end-1,:));
    velyBot32Stepping3Bot = (yPosBot32Stepping3Bot(2:end,:)-yPosBot32Stepping3Bot(1:end-1,:))./(timeBot32Stepping3Bot(2:end,:)-timeBot32Stepping3Bot(1:end-1,:));
    velxBot32Stepping3Bot(isnan(velxBot32Stepping3Bot)) = 0;
    velyBot32Stepping3Bot(isnan(velyBot32Stepping3Bot)) = 0;
    
    absVelBot32Stepping3Bot = sqrt(velxBot32Stepping3Bot.^2 + velyBot32Stepping3Bot.^2)-0.15;
    
if saveFiles == 1
    save('OHCDataExperimentStepping3Bot.mat','xPosLoadStepping3Bot','yPosLoadStepping3Bot',...
        'thetaLoadStepping3Bot','timeLoadStepping3Bot','xPosBot3Stepping3Bot','yPosBot3Stepping3Bot',...
        'thetaBot3Stepping3Bot','timeBot3Stepping3Bot','xPosBot30Stepping3Bot','yPosBot30Stepping3Bot',...
        'thetaBot30Stepping3Bot','timeBot30Stepping3Bot','xPosBot32Stepping3Bot','yPosBot32Stepping3Bot',...
        'thetaBot32Stepping3Bot','timeBot32Stepping3Bot','absVelLoadStepping3Bot','absVelBot30Stepping3Bot',...
        'absVelBot32Stepping3Bot','absVelBot3Stepping3Bot');
end

clear all
close all
clc
% Robot Data
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\3Robot\Stepping
% Determine the largest file
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;

numTrials = 10;

for jj = 1:10
    fname1 = sprintf('experiment%dBot30Stepping8.csv',jj);
    fname2 = sprintf('experiment%dBot32Stepping8.csv',jj);
    fname3 = sprintf('experiment%dBot03Stepping8.csv',jj);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);
    
    maxSizeBot3 = max(maxSizeBot3,length(data3));
    maxSizeBot30 = max(maxSizeBot30,length(data30));
    maxSizeBot32 = max(maxSizeBot32,length(data32));
end
    
    refVelBot3Stepping3Bot = nan(maxSizeBot3,numTrials);
    refVelBot30Stepping3Bot = nan(maxSizeBot30,numTrials);
    refVelBot32Stepping3Bot = nan(maxSizeBot32,numTrials);
    
    forceBot3Stepping3Bot = nan(maxSizeBot3,numTrials);
    forceBot30Stepping3Bot = nan(maxSizeBot30,numTrials);
    forceBot32Stepping3Bot = nan(maxSizeBot32,numTrials);
    
    refTimeBot3Stepping3Bot = nan(maxSizeBot3,numTrials);
    refTimeBot30Stepping3Bot = nan(maxSizeBot30,numTrials);
    refTimeBot32Stepping3Bot = nan(maxSizeBot32,numTrials);
    
for ii = 1:numTrials

    fname1 = sprintf('experiment%dBot30Stepping8.csv',ii);
    fname2 = sprintf('experiment%dBot32Stepping8.csv',ii);
    fname3 = sprintf('experiment%dBot03Stepping8.csv',ii);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);

    %Separate Vel,Force,Time Things
    refVelBot3Stepping3Bot(1:length(data3),ii) = data3(:,1);
    refVelBot30Stepping3Bot(1:length(data30),ii) = data30(:,1);
    refVelBot32Stepping3Bot(1:length(data32),ii) = data32(:,1);
    
    forceBot3Stepping3Bot(1:length(data3),ii) = data3(:,3);
    forceBot30Stepping3Bot(1:length(data30),ii) = data30(:,3);
    forceBot32Stepping3Bot(1:length(data32),ii) = data32(:,3);
    
    refTimeBot3Stepping3Bot(1:length(data3),ii) = data3(:,2)/1000;
    refTimeBot30Stepping3Bot(1:length(data30),ii) = data30(:,2)/1000;
    refTimeBot32Stepping3Bot(1:length(data32),ii) = data32(:,2)/1000;

end

if saveFiles == 1
    save('RobotDataExperimentStepping3Bot.mat','refVelBot30Stepping3Bot','refVelBot32Stepping3Bot','refVelBot3Stepping3Bot',...
        'forceBot30Stepping3Bot','forceBot32Stepping3Bot','forceBot3Stepping3Bot','refTimeBot30Stepping3Bot','refTimeBot32Stepping3Bot','refTimeBot3Stepping3Bot');
end

%% OHC Data for 4 Robots with ML 6 8 10 12
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\4Robot\ML\Speed6_8_10_12

% Determine the largest file
maxSizeLoad = 0;
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;
maxSizeBot46 = 0;

numTrials = 10;

for jj = 1:numTrials
    fname = sprintf('experiment%dSpeed6_8_10_12.csv',jj);
    data = csvread(fname,1,0);
    
    tempLoad = data(data(:,1)==40,:);
    maxSizeLoad = max(maxSizeLoad,length(tempLoad));
    
    tempBot3 = data(data(:,1)==3,:);
    maxSizeBot3 = max(maxSizeBot3,length(tempBot3));
    
    tempBot30 = data(data(:,1)==30,:);
    maxSizeBot30 = max(maxSizeBot30,length(tempBot30));
    
    tempBot32 = data(data(:,1)==32,:);
    maxSizeBot32 = max(maxSizeBot32,length(tempBot32));
    
    tempBot46 = data((data(:,1)==46 | data(:,1)==66),:);
    maxSizeBot46 = max(maxSizeBot46,length(tempBot46)); 
end

    xPosLoad6_8_10_12 = nan(maxSizeLoad,numTrials);yPosLoad6_8_10_12 = nan(maxSizeLoad,numTrials);...
        thetaLoad6_8_10_12 = nan(maxSizeLoad,numTrials); timeLoad6_8_10_12 = nan(maxSizeLoad,numTrials);
    xPosBot36_8_10_12 = nan(maxSizeBot3,numTrials);yPosBot36_8_10_12 = nan(maxSizeBot3,numTrials);...
        thetaBot36_8_10_12 = nan(maxSizeBot3,numTrials); timeBot36_8_10_12 = nan(maxSizeBot3,numTrials);
    xPosBot306_8_10_12 = nan(maxSizeBot30,numTrials);yPosBot306_8_10_12 = nan(maxSizeBot30,numTrials);...
        thetaBot306_8_10_12 =nan(maxSizeBot30,numTrials); timeBot306_8_10_12 = nan(maxSizeBot30,numTrials);
    xPosBot326_8_10_12 = nan(maxSizeBot32,numTrials);yPosBot326_8_10_12 = nan(maxSizeBot32,numTrials);...
        thetaBot326_8_10_12 = nan(maxSizeBot32,numTrials); timeBot326_8_10_12 = nan(maxSizeBot32,numTrials);
    xPosBot466_8_10_12 = nan(maxSizeBot46,numTrials);yPosBot466_8_10_12 = nan(maxSizeBot46,numTrials);...
        thetaBot466_8_10_12 = nan(maxSizeBot46,numTrials); timeBot466_8_10_12 = nan(maxSizeBot46,numTrials);
    
for ii = 1:numTrials
    fname = sprintf('experiment%dSpeed6_8_10_12.csv',ii);
    data = csvread(fname,1,0);

    tempLoad = data(data(:,1)==40,:);
    tempBot30 = data(data(:,1)==30,:);
    tempBot32 = data(data(:,1)==32,:);
    tempBot46 = data((data(:,1)==46 | data(:,1)==66),:);

    %Separate Load Things
    xPosLoad6_8_10_12(1:length(tempLoad(:,2)),ii) = tempLoad(:,2);
    yPosLoad6_8_10_12(1:length(tempLoad(:,3)),ii) = tempLoad(:,3);
    thetaLoad6_8_10_12(1:length(tempLoad(:,4)),ii) = tempLoad(:,4);
    timeLoad6_8_10_12(1:length(tempLoad(:,5)),ii) = tempLoad(:,5);
    
    %Separate Bot3 Things
    xPosBot36_8_10_12(1:length(tempBot3(:,2)),ii) = tempBot3(:,2);
    yPosBot36_8_10_12(1:length(tempBot3(:,3)),ii) = tempBot3(:,3);
    thetaBot36_8_10_12(1:length(tempBot3(:,4)),ii) = tempBot3(:,4);
    timeBot36_8_10_12(1:length(tempBot3(:,5)),ii) = tempBot3(:,5);
    
    %Separate Bot30 Things
    xPosBot306_8_10_12(1:length(tempBot30(:,2)),ii) = tempBot30(:,2);
    yPosBot306_8_10_12(1:length(tempBot30(:,3)),ii) = tempBot30(:,3);
    thetaBot306_8_10_12(1:length(tempBot30(:,4)),ii) = tempBot30(:,4);
    timeBot306_8_10_12(1:length(tempBot30(:,5)),ii) = tempBot30(:,5);
    
    %Separate Bot32 Things
    xPosBot326_8_10_12(1:length(tempBot32(:,2)),ii) = tempBot32(:,2);
    yPosBot326_8_10_12(1:length(tempBot32(:,3)),ii) = tempBot32(:,3);
    thetaBot326_8_10_12(1:length(tempBot32(:,4)),ii) = tempBot32(:,4);
    timeBot326_8_10_12(1:length(tempBot32(:,5)),ii) = tempBot32(:,5);
    
    %Separate Bot46 Things
    xPosBot466_8_10_12(1:length(tempBot46(:,2)),ii) = tempBot46(:,2);
    yPosBot466_8_10_12(1:length(tempBot46(:,3)),ii) = tempBot46(:,3);
    thetaBot466_8_10_12(1:length(tempBot46(:,4)),ii) = tempBot46(:,4);
    timeBot466_8_10_12(1:length(tempBot46(:,5)),ii) = tempBot46(:,5);
end
    
% Determine Velocities

    velxLoad6_8_10_12 = (xPosLoad6_8_10_12(2:end,:)-xPosLoad6_8_10_12(1:end-1,:))./(timeLoad6_8_10_12(2:end,:)-timeLoad6_8_10_12(1:end-1,:));
    velyLoad6_8_10_12 = (yPosLoad6_8_10_12(2:end,:)-yPosLoad6_8_10_12(1:end-1,:))./(timeLoad6_8_10_12(2:end,:)-timeLoad6_8_10_12(1:end-1,:));
    velxLoad6_8_10_12(isnan(velxLoad6_8_10_12)) = 0;
    velyLoad6_8_10_12(isnan(velyLoad6_8_10_12)) = 0;
    
    absVelLoad6_8_10_12 = sqrt(velxLoad6_8_10_12.^2 + velyLoad6_8_10_12.^2);
    
    velxBot36_8_10_12 = (xPosBot36_8_10_12(2:end,:)-xPosBot36_8_10_12(1:end-1,:))./(timeBot36_8_10_12(2:end,:)-timeBot36_8_10_12(1:end-1,:));
    velyBot36_8_10_12 = (yPosBot36_8_10_12(2:end,:)-yPosBot36_8_10_12(1:end-1,:))./(timeBot36_8_10_12(2:end,:)-timeBot36_8_10_12(1:end-1,:));
    velxBot36_8_10_12(isnan(velxBot36_8_10_12)) = 0;
    velyBot36_8_10_12(isnan(velyBot36_8_10_12)) = 0;
    
    absVelBot36_8_10_12 = sqrt(velxBot36_8_10_12.^2 + velyBot36_8_10_12.^2);
    
    velxBot306_8_10_12 = (xPosBot306_8_10_12(2:end,:)-xPosBot306_8_10_12(1:end-1,:))./(timeBot306_8_10_12(2:end,:)-timeBot306_8_10_12(1:end-1,:));
    velyBot306_8_10_12 = (yPosBot306_8_10_12(2:end,:)-yPosBot306_8_10_12(1:end-1,:))./(timeBot306_8_10_12(2:end,:)-timeBot306_8_10_12(1:end-1,:));
    velxBot306_8_10_12(isnan(velxBot306_8_10_12)) = 0;
    velyBot306_8_10_12(isnan(velyBot306_8_10_12)) = 0;
    
    absVelBot306_8_10_12 = sqrt(velxBot306_8_10_12.^2 + velyBot306_8_10_12.^2);
    
    velxBot326_8_10_12 = (xPosBot326_8_10_12(2:end,:)-xPosBot326_8_10_12(1:end-1,:))./(timeBot326_8_10_12(2:end,:)-timeBot326_8_10_12(1:end-1,:));
    velyBot326_8_10_12 = (yPosBot326_8_10_12(2:end,:)-yPosBot326_8_10_12(1:end-1,:))./(timeBot326_8_10_12(2:end,:)-timeBot326_8_10_12(1:end-1,:));
    velxBot326_8_10_12(isnan(velxBot326_8_10_12)) = 0;
    velyBot326_8_10_12(isnan(velyBot326_8_10_12)) = 0;
    
    absVelBot326_8_10_12 = sqrt(velxBot326_8_10_12.^2 + velyBot326_8_10_12.^2);
    
    velxBot466_8_10_12 = (xPosBot466_8_10_12(2:end,:)-xPosBot466_8_10_12(1:end-1,:))./(timeBot466_8_10_12(2:end,:)-timeBot466_8_10_12(1:end-1,:));
    velyBot466_8_10_12 = (yPosBot466_8_10_12(2:end,:)-yPosBot466_8_10_12(1:end-1,:))./(timeBot466_8_10_12(2:end,:)-timeBot466_8_10_12(1:end-1,:));
    velxBot466_8_10_12(isnan(velxBot466_8_10_12)) = 0;
    velyBot466_8_10_12(isnan(velyBot466_8_10_12)) = 0;
    
    absVelBot466_8_10_12 = sqrt(velxBot466_8_10_12.^2 + velyBot466_8_10_12.^2);
    
if saveFiles == 1
    save('OHCDataExperimentSpeed6_8_10_12.mat','xPosLoad6_8_10_12','yPosLoad6_8_10_12',...
        'thetaLoad6_8_10_12','timeLoad6_8_10_12','xPosBot36_8_10_12','yPosBot36_8_10_12',...
        'thetaBot36_8_10_12','timeBot36_8_10_12','xPosBot306_8_10_12','yPosBot306_8_10_12',...
        'thetaBot306_8_10_12','timeBot306_8_10_12','xPosBot326_8_10_12','yPosBot326_8_10_12',...
        'thetaBot326_8_10_12','timeBot326_8_10_12','xPosBot466_8_10_12','yPosBot466_8_10_12',...
        'thetaBot466_8_10_12','timeBot466_8_10_12','absVelLoad6_8_10_12','absVelBot306_8_10_12',...
        'absVelBot326_8_10_12','absVelBot36_8_10_12','absVelBot466_8_10_12');
end

clear all
close all
clc
% Robot Data
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\4Robot\ML\Speed6_8_10_12
% Determine the largest file
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;
maxSizeBot46 = 0;

numTrials = 10;

for jj = 1:10
    fname1 = sprintf('experiment%dBot30Speed6_8_10_12.csv',jj);
    fname2 = sprintf('experiment%dBot32Speed6_8_10_12.csv',jj);
    fname3 = sprintf('experiment%dBot03Speed6_8_10_12.csv',jj);
    fname4 = sprintf('experiment%dBot46Speed6_8_10_12.csv',jj);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);
    data46 = csvread(fname4,1,0);
    
    maxSizeBot3 = max(maxSizeBot3,length(data3));
    maxSizeBot30 = max(maxSizeBot30,length(data30));
    maxSizeBot32 = max(maxSizeBot32,length(data32));
    maxSizeBot46 = max(maxSizeBot46,length(data46));
end
    
    refVelBot36_8_10_12 = nan(maxSizeBot3,numTrials);
    refVelBot306_8_10_12 = nan(maxSizeBot30,numTrials);
    refVelBot326_8_10_12 = nan(maxSizeBot32,numTrials);
    refVelBot466_8_10_12 = nan(maxSizeBot46,numTrials);
    
    forceBot36_8_10_12 = nan(maxSizeBot3,numTrials);
    forceBot306_8_10_12 = nan(maxSizeBot30,numTrials);
    forceBot326_8_10_12 = nan(maxSizeBot32,numTrials);
    forceBot466_8_10_12 = nan(maxSizeBot46,numTrials);
    
    refTimeBot36_8_10_12 = nan(maxSizeBot3,numTrials);
    refTimeBot306_8_10_12 = nan(maxSizeBot30,numTrials);
    refTimeBot326_8_10_12 = nan(maxSizeBot32,numTrials);
    refTimeBot466_8_10_12 = nan(maxSizeBot46,numTrials);
    
for ii = 1:numTrials

    fname1 = sprintf('experiment%dBot30Speed6_8_10_12.csv',ii);
    fname2 = sprintf('experiment%dBot32Speed6_8_10_12.csv',ii);
    fname3 = sprintf('experiment%dBot03Speed6_8_10_12.csv',ii);
    fname4 = sprintf('experiment%dBot46Speed6_8_10_12.csv',ii);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);
    data46 = csvread(fname4,1,0);

    %Separate Vel,Force,Time Things
    refVelBot36_8_10_12(1:length(data3),ii) = data3(:,1);
    refVelBot306_8_10_12(1:length(data30),ii) = data30(:,1);
    refVelBot326_8_10_12(1:length(data32),ii) = data32(:,1);
    refVelBot466_8_10_12(1:length(data46),ii) = data46(:,1) + 1.1;
    
    forceBot36_8_10_12(1:length(data3),ii) = data3(:,3);
    forceBot306_8_10_12(1:length(data30),ii) = data30(:,3);
    forceBot326_8_10_12(1:length(data32),ii) = data32(:,3);
    forceBot466_8_10_12(1:length(data46),ii) = data46(:,3);
    
    refTimeBot36_8_10_12(1:length(data3),ii) = data3(:,2)/1000;
    refTimeBot306_8_10_12(1:length(data30),ii) = data30(:,2)/1000;
    refTimeBot326_8_10_12(1:length(data32),ii) = data32(:,2)/1000;
    refTimeBot466_8_10_12(1:length(data46),ii) = data46(:,2)/1000;
end

if saveFiles == 1
    save('RobotDataExperimentSpeed6_8_10_12.mat','refVelBot306_8_10_12','refVelBot326_8_10_12',...
        'refVelBot36_8_10_12','refVelBot466_8_10_12','forceBot306_8_10_12','forceBot326_8_10_12',...
        'forceBot36_8_10_12','forceBot466_8_10_12',...
        'refTimeBot306_8_10_12','refTimeBot326_8_10_12','refTimeBot36_8_10_12','refTimeBot466_8_10_12');
end

%% OHC Data for 4 Robots with ML 4 6 9 12
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\4Robot\ML\Speed4_6_9_12

% Determine the largest file
maxSizeLoad = 0;
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;
maxSizeBot46 = 0;

numTrials = 10;

for jj = 1:numTrials
    fname = sprintf('experiment%dSpeed4_6_9_12.csv',jj);
    data = csvread(fname,1,0);
    
    tempLoad = data(data(:,1)==40,:);
    maxSizeLoad = max(maxSizeLoad,length(tempLoad));
    %loadIndex = [find(timeLoad4_6_9_12(:,1)>=20,1), find(timeLoad4_6_9_12(:,1)>=100,1)];
    
    tempBot3 = data(data(:,1)==3,:);
    maxSizeBot3 = max(maxSizeBot3,length(tempBot3));
    
    tempBot30 = data(data(:,1)==30,:);
    maxSizeBot30 = max(maxSizeBot30,length(tempBot30));
    
    tempBot32 = data(data(:,1)==32,:);
    maxSizeBot32 = max(maxSizeBot32,length(tempBot32));
    
    tempBot46 = data((data(:,1)==46 | data(:,1)==66),:);
    maxSizeBot46 = max(maxSizeBot46,length(tempBot46)); 
end

    xPosLoad4_6_9_12 = nan(maxSizeLoad,numTrials);yPosLoad4_6_9_12 = nan(maxSizeLoad,numTrials);...
        thetaLoad4_6_9_12 = nan(maxSizeLoad,numTrials); timeLoad4_6_9_12 = nan(maxSizeLoad,numTrials);
    xPosBot34_6_9_12 = nan(maxSizeBot3,numTrials);yPosBot34_6_9_12 = nan(maxSizeBot3,numTrials);...
        thetaBot34_6_9_12 = nan(maxSizeBot3,numTrials); timeBot34_6_9_12 = nan(maxSizeBot3,numTrials);
    xPosBot304_6_9_12 = nan(maxSizeBot30,numTrials);yPosBot304_6_9_12 = nan(maxSizeBot30,numTrials);...
        thetaBot304_6_9_12 =nan(maxSizeBot30,numTrials); timeBot304_6_9_12 = nan(maxSizeBot30,numTrials);
    xPosBot324_6_9_12 = nan(maxSizeBot32,numTrials);yPosBot324_6_9_12 = nan(maxSizeBot32,numTrials);...
        thetaBot324_6_9_12 = nan(maxSizeBot32,numTrials); timeBot324_6_9_12 = nan(maxSizeBot32,numTrials);
    xPosBot464_6_9_12 = nan(maxSizeBot46,numTrials);yPosBot464_6_9_12 = nan(maxSizeBot46,numTrials);...
        thetaBot464_6_9_12 = nan(maxSizeBot46,numTrials); timeBot464_6_9_12 = nan(maxSizeBot46,numTrials);
    
for ii = 1:numTrials
    fname = sprintf('experiment%dSpeed4_6_9_12.csv',ii);
    data = csvread(fname,1,0);

    tempLoad = data(data(:,1)==40,:);
    tempBot30 = data(data(:,1)==30,:);
    tempBot32 = data(data(:,1)==32,:);
    tempBot46 = data((data(:,1)==46 | data(:,1)==66),:);

    %Separate Load Things
    xPosLoad4_6_9_12(1:length(tempLoad(:,2)),ii) = tempLoad(:,2);
    yPosLoad4_6_9_12(1:length(tempLoad(:,3)),ii) = tempLoad(:,3);
    thetaLoad4_6_9_12(1:length(tempLoad(:,4)),ii) = tempLoad(:,4);
    timeLoad4_6_9_12(1:length(tempLoad(:,5)),ii) = tempLoad(:,5);
    
    %Separate Bot3 Things
    xPosBot34_6_9_12(1:length(tempBot3(:,2)),ii) = tempBot3(:,2);
    yPosBot34_6_9_12(1:length(tempBot3(:,3)),ii) = tempBot3(:,3);
    thetaBot34_6_9_12(1:length(tempBot3(:,4)),ii) = tempBot3(:,4);
    timeBot34_6_9_12(1:length(tempBot3(:,5)),ii) = tempBot3(:,5);
    
    %Separate Bot30 Things
    xPosBot304_6_9_12(1:length(tempBot30(:,2)),ii) = tempBot30(:,2);
    yPosBot304_6_9_12(1:length(tempBot30(:,3)),ii) = tempBot30(:,3);
    thetaBot304_6_9_12(1:length(tempBot30(:,4)),ii) = tempBot30(:,4);
    timeBot304_6_9_12(1:length(tempBot30(:,5)),ii) = tempBot30(:,5);
    
    %Separate Bot32 Things
    xPosBot324_6_9_12(1:length(tempBot32(:,2)),ii) = tempBot32(:,2);
    yPosBot324_6_9_12(1:length(tempBot32(:,3)),ii) = tempBot32(:,3);
    thetaBot324_6_9_12(1:length(tempBot32(:,4)),ii) = tempBot32(:,4);
    timeBot324_6_9_12(1:length(tempBot32(:,5)),ii) = tempBot32(:,5);
    
    %Separate Bot46 Things
    xPosBot464_6_9_12(1:length(tempBot46(:,2)),ii) = tempBot46(:,2);
    yPosBot464_6_9_12(1:length(tempBot46(:,3)),ii) = tempBot46(:,3);
    thetaBot464_6_9_12(1:length(tempBot46(:,4)),ii) = tempBot46(:,4);
    timeBot464_6_9_12(1:length(tempBot46(:,5)),ii) = tempBot46(:,5);
end
    
% Determine Velocities

    velxLoad4_6_9_12 = (xPosLoad4_6_9_12(2:end,:)-xPosLoad4_6_9_12(1:end-1,:))./(timeLoad4_6_9_12(2:end,:)-timeLoad4_6_9_12(1:end-1,:));
    velyLoad4_6_9_12 = (yPosLoad4_6_9_12(2:end,:)-yPosLoad4_6_9_12(1:end-1,:))./(timeLoad4_6_9_12(2:end,:)-timeLoad4_6_9_12(1:end-1,:));
    velxLoad4_6_9_12(isnan(velxLoad4_6_9_12)) = 0;
    velyLoad4_6_9_12(isnan(velyLoad4_6_9_12)) = 0;
    
    absVelLoad4_6_9_12 = sqrt(velxLoad4_6_9_12.^2 + velyLoad4_6_9_12.^2);
    
    velxBot34_6_9_12 = (xPosBot34_6_9_12(2:end,:)-xPosBot34_6_9_12(1:end-1,:))./(timeBot34_6_9_12(2:end,:)-timeBot34_6_9_12(1:end-1,:));
    velyBot34_6_9_12 = (yPosBot34_6_9_12(2:end,:)-yPosBot34_6_9_12(1:end-1,:))./(timeBot34_6_9_12(2:end,:)-timeBot34_6_9_12(1:end-1,:));
    velxBot34_6_9_12(isnan(velxBot34_6_9_12)) = 0;
    velyBot34_6_9_12(isnan(velyBot34_6_9_12)) = 0;
    
    absVelBot34_6_9_12 = sqrt(velxBot34_6_9_12.^2 + velyBot34_6_9_12.^2);
    
    velxBot304_6_9_12 = (xPosBot304_6_9_12(2:end,:)-xPosBot304_6_9_12(1:end-1,:))./(timeBot304_6_9_12(2:end,:)-timeBot304_6_9_12(1:end-1,:));
    velyBot304_6_9_12 = (yPosBot304_6_9_12(2:end,:)-yPosBot304_6_9_12(1:end-1,:))./(timeBot304_6_9_12(2:end,:)-timeBot304_6_9_12(1:end-1,:));
    velxBot304_6_9_12(isnan(velxBot304_6_9_12)) = 0;
    velyBot304_6_9_12(isnan(velyBot304_6_9_12)) = 0;
    
    absVelBot304_6_9_12 = sqrt(velxBot304_6_9_12.^2 + velyBot304_6_9_12.^2);
    
    velxBot324_6_9_12 = (xPosBot324_6_9_12(2:end,:)-xPosBot324_6_9_12(1:end-1,:))./(timeBot324_6_9_12(2:end,:)-timeBot324_6_9_12(1:end-1,:));
    velyBot324_6_9_12 = (yPosBot324_6_9_12(2:end,:)-yPosBot324_6_9_12(1:end-1,:))./(timeBot324_6_9_12(2:end,:)-timeBot324_6_9_12(1:end-1,:));
    velxBot324_6_9_12(isnan(velxBot324_6_9_12)) = 0;
    velyBot324_6_9_12(isnan(velyBot324_6_9_12)) = 0;
    
    absVelBot324_6_9_12 = sqrt(velxBot324_6_9_12.^2 + velyBot324_6_9_12.^2);
    
    velxBot464_6_9_12 = (xPosBot464_6_9_12(2:end,:)-xPosBot464_6_9_12(1:end-1,:))./(timeBot464_6_9_12(2:end,:)-timeBot464_6_9_12(1:end-1,:));
    velyBot464_6_9_12 = (yPosBot464_6_9_12(2:end,:)-yPosBot464_6_9_12(1:end-1,:))./(timeBot464_6_9_12(2:end,:)-timeBot464_6_9_12(1:end-1,:));
    velxBot464_6_9_12(isnan(velxBot464_6_9_12)) = 0;
    velyBot464_6_9_12(isnan(velyBot464_6_9_12)) = 0;
    
    absVelBot464_6_9_12 = sqrt(velxBot464_6_9_12.^2 + velyBot464_6_9_12.^2);
    
if saveFiles == 1
    save('OHCDataExperimentSpeed4_6_9_12.mat','xPosLoad4_6_9_12','yPosLoad4_6_9_12',...
        'thetaLoad4_6_9_12','timeLoad4_6_9_12','xPosBot34_6_9_12','yPosBot34_6_9_12',...
        'thetaBot34_6_9_12','timeBot34_6_9_12','xPosBot304_6_9_12','yPosBot304_6_9_12',...
        'thetaBot304_6_9_12','timeBot304_6_9_12','xPosBot324_6_9_12','yPosBot324_6_9_12',...
        'thetaBot324_6_9_12','timeBot324_6_9_12','xPosBot464_6_9_12','yPosBot464_6_9_12',...
        'thetaBot464_6_9_12','timeBot464_6_9_12','absVelLoad4_6_9_12','absVelBot304_6_9_12',...
        'absVelBot324_6_9_12','absVelBot34_6_9_12','absVelBot464_6_9_12');
end

clear all
close all
clc
% Robot Data
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\4Robot\ML\Speed4_6_9_12
% Determine the largest file
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;
maxSizeBot46 = 0;

numTrials = 10;

for jj = 1:10
    fname1 = sprintf('experiment%dBot30Speed4_6_9_12.csv',jj);
    fname2 = sprintf('experiment%dBot32Speed4_6_9_12.csv',jj);
    fname3 = sprintf('experiment%dBot03Speed4_6_9_12.csv',jj);
    fname4 = sprintf('experiment%dBot46Speed4_6_9_12.csv',jj);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);
    data46 = csvread(fname4,1,0);
    
    maxSizeBot3 = max(maxSizeBot3,length(data3));
    maxSizeBot30 = max(maxSizeBot30,length(data30));
    maxSizeBot32 = max(maxSizeBot32,length(data32));
    maxSizeBot46 = max(maxSizeBot46,length(data46));
end
    
    refVelBot34_6_9_12 = nan(maxSizeBot3,numTrials);
    refVelBot304_6_9_12 = nan(maxSizeBot30,numTrials);
    refVelBot324_6_9_12 = nan(maxSizeBot32,numTrials);
    refVelBot464_6_9_12 = nan(maxSizeBot46,numTrials);
    
    forceBot34_6_9_12 = nan(maxSizeBot3,numTrials);
    forceBot304_6_9_12 = nan(maxSizeBot30,numTrials);
    forceBot324_6_9_12 = nan(maxSizeBot32,numTrials);
    forceBot464_6_9_12 = nan(maxSizeBot46,numTrials);
    
    refTimeBot34_6_9_12 = nan(maxSizeBot3,numTrials);
    refTimeBot304_6_9_12 = nan(maxSizeBot30,numTrials);
    refTimeBot324_6_9_12 = nan(maxSizeBot32,numTrials);
    refTimeBot464_6_9_12 = nan(maxSizeBot46,numTrials);
    
for ii = 1:numTrials

    fname1 = sprintf('experiment%dBot30Speed4_6_9_12.csv',ii);
    fname2 = sprintf('experiment%dBot32Speed4_6_9_12.csv',ii);
    fname3 = sprintf('experiment%dBot03Speed4_6_9_12.csv',ii);
    fname4 = sprintf('experiment%dBot46Speed4_6_9_12.csv',ii);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);
    data46 = csvread(fname4,1,0);

    %Separate Vel,Force,Time Things
    refVelBot34_6_9_12(1:length(data3),ii) = data3(:,1);
    refVelBot304_6_9_12(1:length(data30),ii) = data30(:,1);
    refVelBot324_6_9_12(1:length(data32),ii) = data32(:,1);
    refVelBot464_6_9_12(1:length(data46),ii) = data46(:,1) + 1.1;
    
    forceBot34_6_9_12(1:length(data3),ii) = data3(:,3);
    forceBot304_6_9_12(1:length(data30),ii) = data30(:,3);
    forceBot324_6_9_12(1:length(data32),ii) = data32(:,3);
    forceBot464_6_9_12(1:length(data46),ii) = data46(:,3);
    
    refTimeBot34_6_9_12(1:length(data3),ii) = data3(:,2)/1000;
    refTimeBot304_6_9_12(1:length(data30),ii) = data30(:,2)/1000;
    refTimeBot324_6_9_12(1:length(data32),ii) = data32(:,2)/1000;
    refTimeBot464_6_9_12(1:length(data46),ii) = data46(:,2)/1000;
end

if saveFiles == 1
    save('RobotDataExperimentSpeed4_6_9_12.mat','refVelBot304_6_9_12','refVelBot324_6_9_12',...
        'refVelBot34_6_9_12','refVelBot464_6_9_12','forceBot304_6_9_12','forceBot324_6_9_12',...
        'forceBot34_6_9_12','forceBot464_6_9_12',...
        'refTimeBot304_6_9_12','refTimeBot324_6_9_12','refTimeBot34_6_9_12','refTimeBot464_6_9_12');
end


%% OHC Data for 4 Robots Stepping
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\OHCData\4Robot\Stepping

% Determine the largest file
maxSizeLoad = 0;
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;
maxSizeBot46 = 0;

numTrials = 10;

for jj = 1:10
    fname = sprintf('experiment%dStepping8.csv',jj);
    data = csvread(fname,1,0);
    
    tempLoad = data(data(:,1)==40,:);
    maxSizeLoad = max(maxSizeLoad,length(tempLoad));
    
    tempBot3 = data(data(:,1)==3,:);
    maxSizeBot3 = max(maxSizeBot3,length(tempBot3));
    
    tempBot30 = data(data(:,1)==30,:);
    maxSizeBot30 = max(maxSizeBot30,length(tempBot30));
    
    tempBot32 = data(data(:,1)==32,:);
    maxSizeBot32 = max(maxSizeBot32,length(tempBot32));
    
    tempBot46 = data(data(:,1)==(46 | 66),:);
    maxSizeBot46 = max(maxSizeBot46,length(tempBot46)); 
end

    xPosLoadStepping4Bot = nan(maxSizeLoad,numTrials);yPosLoadStepping4Bot = nan(maxSizeLoad,numTrials);...
        thetaLoadStepping4Bot = nan(maxSizeLoad,numTrials); timeLoadStepping4Bot = nan(maxSizeLoad,numTrials);
    xPosBot3Stepping4Bot = nan(maxSizeBot3,numTrials);yPosBot3Stepping4Bot = nan(maxSizeBot3,numTrials);...
        thetaBot3Stepping4Bot =nan(maxSizeBot3,numTrials); timeBot3Stepping4Bot = nan(maxSizeBot3,numTrials);
    xPosBot30Stepping4Bot = nan(maxSizeBot30,numTrials);yPosBot30Stepping4Bot = nan(maxSizeBot30,numTrials);...
        thetaBot30Stepping4Bot =nan(maxSizeBot30,numTrials); timeBot30Stepping4Bot = nan(maxSizeBot30,numTrials);
    xPosBot32Stepping4Bot = nan(maxSizeBot32,numTrials);yPosBot32Stepping4Bot = nan(maxSizeBot32,numTrials);...
        thetaBot32Stepping4Bot = nan(maxSizeBot32,numTrials); timeBot32Stepping4Bot = nan(maxSizeBot32,numTrials);
    xPosBot46Stepping4Bot = nan(maxSizeBot46,numTrials);yPosBot46Stepping4Bot = nan(maxSizeBot46,numTrials);...
        thetaBot46Stepping4Bot = nan(maxSizeBot46,numTrials); timeBot46Stepping4Bot = nan(maxSizeBot46,numTrials);
    
for ii = 1:numTrials
    fname = sprintf('experiment%dStepping8.csv',ii);
    data = csvread(fname,1,0);

    tempLoad = data(data(:,1)==40,:);
    tempBot3 = data(data(:,1)==3,:);
    tempBot30 = data(data(:,1)==30,:);
    tempBot32 = data(data(:,1)==32,:);
    tempBot46 = data(data(:,1)==(46 | 66),:);

    %Separate Load Things
    xPosLoadStepping4Bot(1:length(tempLoad(:,2)),ii) = tempLoad(:,2);
    yPosLoadStepping4Bot(1:length(tempLoad(:,3)),ii) = tempLoad(:,3);
    thetaLoadStepping4Bot(1:length(tempLoad(:,4)),ii) = tempLoad(:,4);
    timeLoadStepping4Bot(1:length(tempLoad(:,5)),ii) = tempLoad(:,5);
    
    %Separate Bot3 Things
    xPosBot3Stepping4Bot(1:length(tempBot3(:,2)),ii) = tempBot3(:,2);
    yPosBot3Stepping4Bot(1:length(tempBot3(:,3)),ii) = tempBot3(:,3);
    thetaBot3Stepping4Bot(1:length(tempBot3(:,4)),ii) = tempBot3(:,4);
    timeBot3Stepping4Bot(1:length(tempBot3(:,5)),ii) = tempBot3(:,5);
    
    %Separate Bot30 Things
    xPosBot30Stepping4Bot(1:length(tempBot30(:,2)),ii) = tempBot30(:,2);
    yPosBot30Stepping4Bot(1:length(tempBot30(:,3)),ii) = tempBot30(:,3);
    thetaBot30Stepping4Bot(1:length(tempBot30(:,4)),ii) = tempBot30(:,4);
    timeBot30Stepping4Bot(1:length(tempBot30(:,5)),ii) = tempBot30(:,5);
    
    %Separate Bot32 Things
    xPosBot32Stepping4Bot(1:length(tempBot32(:,2)),ii) = tempBot32(:,2);
    yPosBot32Stepping4Bot(1:length(tempBot32(:,3)),ii) = tempBot32(:,3);
    thetaBot32Stepping4Bot(1:length(tempBot32(:,4)),ii) = tempBot32(:,4);
    timeBot32Stepping4Bot(1:length(tempBot32(:,5)),ii) = tempBot32(:,5);
    
    %Separate Bot46 Things
    xPosBot46Stepping4Bot(1:length(tempBot46(:,2)),ii) = tempBot46(:,2);
    yPosBot46Stepping4Bot(1:length(tempBot46(:,3)),ii) = tempBot46(:,3);
    thetaBot46Stepping4Bot(1:length(tempBot46(:,4)),ii) = tempBot46(:,4);
    timeBot46Stepping4Bot(1:length(tempBot46(:,5)),ii) = tempBot46(:,5);
end
    
% Determine Velocities

    velxLoadStepping4Bot = (xPosLoadStepping4Bot(2:end,:)-xPosLoadStepping4Bot(1:end-1,:))./(timeLoadStepping4Bot(2:end,:)-timeLoadStepping4Bot(1:end-1,:));
    velyLoadStepping4Bot = (yPosLoadStepping4Bot(2:end,:)-yPosLoadStepping4Bot(1:end-1,:))./(timeLoadStepping4Bot(2:end,:)-timeLoadStepping4Bot(1:end-1,:));
    velxLoadStepping4Bot(isnan(velxLoadStepping4Bot)) = 0;
    velyLoadStepping4Bot(isnan(velyLoadStepping4Bot)) = 0;
    
    absVelLoadStepping4Bot = sqrt(velxLoadStepping4Bot.^2 + velyLoadStepping4Bot.^2);
    
    velxBot3Stepping4Bot = (xPosBot3Stepping4Bot(2:end,:)-xPosBot3Stepping4Bot(1:end-1,:))./(timeBot3Stepping4Bot(2:end,:)-timeBot3Stepping4Bot(1:end-1,:));
    velyBot3Stepping4Bot = (yPosBot3Stepping4Bot(2:end,:)-yPosBot3Stepping4Bot(1:end-1,:))./(timeBot3Stepping4Bot(2:end,:)-timeBot3Stepping4Bot(1:end-1,:));
    velxBot3Stepping4Bot(isnan(velxBot3Stepping4Bot)) = 0;
    velyBot3Stepping4Bot(isnan(velyBot3Stepping4Bot)) = 0;
    
    absVelBot3Stepping4Bot = sqrt(velxBot3Stepping4Bot.^2 + velyBot3Stepping4Bot.^2);
    
    velxBot30Stepping4Bot = (xPosBot30Stepping4Bot(2:end,:)-xPosBot30Stepping4Bot(1:end-1,:))./(timeBot30Stepping4Bot(2:end,:)-timeBot30Stepping4Bot(1:end-1,:));
    velyBot30Stepping4Bot = (yPosBot30Stepping4Bot(2:end,:)-yPosBot30Stepping4Bot(1:end-1,:))./(timeBot30Stepping4Bot(2:end,:)-timeBot30Stepping4Bot(1:end-1,:));
    velxBot30Stepping4Bot(isnan(velxBot30Stepping4Bot)) = 0;
    velyBot30Stepping4Bot(isnan(velyBot30Stepping4Bot)) = 0;
    
    absVelBot30Stepping4Bot = sqrt(velxBot30Stepping4Bot.^2 + velyBot30Stepping4Bot.^2);
    
    velxBot32Stepping4Bot = (xPosBot32Stepping4Bot(2:end,:)-xPosBot32Stepping4Bot(1:end-1,:))./(timeBot32Stepping4Bot(2:end,:)-timeBot32Stepping4Bot(1:end-1,:));
    velyBot32Stepping4Bot = (yPosBot32Stepping4Bot(2:end,:)-yPosBot32Stepping4Bot(1:end-1,:))./(timeBot32Stepping4Bot(2:end,:)-timeBot32Stepping4Bot(1:end-1,:));
    velxBot32Stepping4Bot(isnan(velxBot32Stepping4Bot)) = 0;
    velyBot32Stepping4Bot(isnan(velyBot32Stepping4Bot)) = 0;
    
    absVelBot32Stepping4Bot = sqrt(velxBot32Stepping4Bot.^2 + velyBot32Stepping4Bot.^2);
    
    velxBot46Stepping4Bot = (xPosBot46Stepping4Bot(2:end,:)-xPosBot46Stepping4Bot(1:end-1,:))./(timeBot46Stepping4Bot(2:end,:)-timeBot46Stepping4Bot(1:end-1,:));
    velyBot46Stepping4Bot = (yPosBot46Stepping4Bot(2:end,:)-yPosBot46Stepping4Bot(1:end-1,:))./(timeBot46Stepping4Bot(2:end,:)-timeBot46Stepping4Bot(1:end-1,:));
    velxBot46Stepping4Bot(isnan(velxBot46Stepping4Bot)) = 0;
    velyBot46Stepping4Bot(isnan(velyBot46Stepping4Bot)) = 0;
    
    absVelBot46Stepping4Bot = sqrt(velxBot46Stepping4Bot.^2 + velyBot46Stepping4Bot.^2);
    
if saveFiles == 1
    save('OHCDataExperimentStepping4Bot.mat','xPosLoadStepping4Bot','yPosLoadStepping4Bot',...
        'thetaLoadStepping4Bot','timeLoadStepping4Bot','xPosBot3Stepping4Bot','yPosBot3Stepping4Bot',...
        'thetaBot3Stepping4Bot','timeBot3Stepping4Bot','xPosBot30Stepping4Bot','yPosBot30Stepping4Bot',...
        'thetaBot30Stepping4Bot','timeBot30Stepping4Bot','xPosBot32Stepping4Bot','yPosBot32Stepping4Bot',...
        'thetaBot32Stepping4Bot','timeBot32Stepping4Bot','xPosBot46Stepping4Bot','yPosBot46Stepping4Bot',...
        'thetaBot46Stepping4Bot','timeBot46Stepping4Bot','absVelLoadStepping4Bot','absVelBot30Stepping4Bot',...
        'absVelBot32Stepping4Bot','absVelBot3Stepping4Bot','absVelBot46Stepping4Bot');
end

clear all
close all
clc
% Robot Data
saveFiles = 1; % Save the matfiles if 1.

cd E:\SeanResearch\PheenoChariot\TrackingCode\DataForPaper\RobotData\4Robot\Stepping
% Determine the largest file
maxSizeBot3 = 0;
maxSizeBot30 = 0;
maxSizeBot32 = 0;
maxSizeBot46 = 0;

numTrials = 10;

for jj = 1:10
    fname1 = sprintf('experiment%dBot30Stepping8.csv',jj);
    fname2 = sprintf('experiment%dBot32Stepping8.csv',jj);
    fname3 = sprintf('experiment%dBot03Stepping8.csv',jj);
    fname4 = sprintf('experiment%dBot46Stepping8.csv',jj);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);
    data46 = csvread(fname4,1,0);
    
    maxSizeBot3 = max(maxSizeBot3,length(data3));
    maxSizeBot30 = max(maxSizeBot30,length(data30));
    maxSizeBot32 = max(maxSizeBot32,length(data32));
    maxSizeBot46 = max(maxSizeBot46,length(data46));
end
    
    refVelBot3Stepping4Bot = nan(maxSizeBot3,numTrials);
    refVelBot30Stepping4Bot = nan(maxSizeBot30,numTrials);
    refVelBot32Stepping4Bot = nan(maxSizeBot32,numTrials);
    refVelBot46Stepping4Bot = nan(maxSizeBot46,numTrials);
    
    forceBot3Stepping4Bot = nan(maxSizeBot3,numTrials);
    forceBot30Stepping4Bot = nan(maxSizeBot30,numTrials);
    forceBot32Stepping4Bot = nan(maxSizeBot32,numTrials);
    forceBot46Stepping4Bot = nan(maxSizeBot46,numTrials);
    
    refTimeBot3Stepping4Bot = nan(maxSizeBot3,numTrials);
    refTimeBot30Stepping4Bot = nan(maxSizeBot30,numTrials);
    refTimeBot32Stepping4Bot = nan(maxSizeBot32,numTrials);
    refTimeBot46Stepping4Bot = nan(maxSizeBot46,numTrials);
    
for ii = 1:numTrials

    fname1 = sprintf('experiment%dBot30Stepping8.csv',ii);
    fname2 = sprintf('experiment%dBot32Stepping8.csv',ii);
    fname3 = sprintf('experiment%dBot03Stepping8.csv',ii);
    fname4 = sprintf('experiment%dBot46Stepping8.csv',ii);
    data30 = csvread(fname1,1,0);
    data32 = csvread(fname2,1,0);
    data3 = csvread(fname3,1,0);
    data46 = csvread(fname4,1,0);

    %Separate Vel,Force,Time Things
    refVelBot3Stepping4Bot(1:length(data3),ii) = data3(:,1);
    refVelBot30Stepping4Bot(1:length(data30),ii) = data30(:,1);
    refVelBot32Stepping4Bot(1:length(data32),ii) = data32(:,1);
    refVelBot46Stepping4Bot(1:length(data46),ii) = data46(:,1);
    
    forceBot3Stepping4Bot(1:length(data3),ii) = data3(:,3);
    forceBot30Stepping4Bot(1:length(data30),ii) = data30(:,3);
    forceBot32Stepping4Bot(1:length(data32),ii) = data32(:,3);
    forceBot46Stepping4Bot(1:length(data46),ii) = data46(:,3);
    
    refTimeBot3Stepping4Bot(1:length(data3),ii) = data3(:,2)/1000;
    refTimeBot30Stepping4Bot(1:length(data30),ii) = data30(:,2)/1000;
    refTimeBot32Stepping4Bot(1:length(data32),ii) = data32(:,2)/1000;
    refTimeBot46Stepping4Bot(1:length(data46),ii) = data46(:,2)/1000;

end

if saveFiles == 1
    save('RobotDataExperimentStepping4Bot.mat','refVelBot30Stepping4Bot','refVelBot32Stepping4Bot','refVelBot3Stepping4Bot',...
        'refVelBot46Stepping4Bot','forceBot30Stepping4Bot','forceBot32Stepping4Bot','forceBot3Stepping4Bot',...
        'forceBot46Stepping4Bot','refTimeBot30Stepping4Bot','refTimeBot32Stepping4Bot','refTimeBot3Stepping4Bot',...
        'refTimeBot46Stepping4Bot');
end