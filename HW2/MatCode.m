% EECS 395/495: Introduction to Computational Photography
% HW2: Measuring the sensor noise of your Tegra tablet

% Student Name: Haikun Liu
% Student Number: 2903021
% NetID: hlg483

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1. Loading the raw images into MATLAB

clear;
figureNo = 0;
NumOfIms = 49; %There are totally 49 RAW pictures captured for each sensitivity value

% Note that: The highest camera sensiticity value is 1600 and the lowest
% camera sensitivity value is 100

% High sensitivity
str_Path = '/Users/HKLHK/Desktop/2015 Fall Quarter (NU)/EECS495 Introduction to Computational Photography/HW/HW2/high Sensitivity 1600 (2)/';
for i = 1: NumOfIms
str_Load = strcat(str_Path, num2str(i-1), '.dng');%set image path for high sensiticity pictures
t = Tiff(str_Load, 'r');
I = read(t); %read RAW images
Im = imcrop(I,[1275,635,949 699]);%crop picture into smaller size image 950*700 centered at (1275,635)
ImageHS(:,i) = double(reshape(Im, [ ], 1));%for each image, rearrange its pixels into a column vector. (make it easy for future manipulation)
end

% Middle sensitivity = (High sensitivity + Low sensiticity)/2
figureNo = 0;
NumOfIms = 49;%There are totally 49 RAW pictures captured for each sensitivity value
str_Path = '/Users/HKLHK/Desktop/2015 Fall Quarter (NU)/EECS495 Introduction to Computational Photography/HW/HW2/Middle Sensitivity 850 (2)/';
for i = 1: NumOfIms
str_Load = strcat(str_Path, num2str(i-1), '.dng');%set image path for high sensiticity pictures
t = Tiff(str_Load, 'r');
I = read(t); %read RAW images
Im = imcrop(I,[1275,635,949 699]);%crop picture into smaller size image 950*700 centered at (1275,635)
ImageMS(:,i) = double(reshape(Im, [ ], 1));%for each image, rearrange its pixels into a column vector. (make it easy for future manipulation)
end

% Low sensiticity
figureNo = 0;
NumOfIms = 49;%There are totally 49 RAW pictures captured for each sensitivity value
str_Path = '/Users/HKLHK/Desktop/2015 Fall Quarter (NU)/EECS495 Introduction to Computational Photography/HW/HW2/Low Sensitivity 100 (2)/';
for i = 1: NumOfIms
str_Load = strcat(str_Path, num2str(i-1), '.dng');%set image path for low sensiticity pictures
t = Tiff(str_Load, 'r');
I = read(t); %read RAW images

Im = imcrop(I,[1275,635,949 699]);%crop picture into smaller size image 950*700 centered at (1275,635)
ImageLS(:,i) = double(reshape(Im, [ ], 1));%store pixel value into column vetor
end

clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 2. Plot a histogram of pixel values for a given pixel 

% High Sensitivity

figureNo = figureNo+1;
figure(figureNo);

subplot(2,2,1);
hist(ImageHS(221667,:)); %Pixel 221667
title('Histogram of Pixel 221667 (High Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');
subplot(2,2,2);
hist(ImageHS(250167,:));%Pixel 250167
title('Histogram of Pixel 250167 (High Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');
subplot(2,2,3);
hist(ImageHS(221967,:));%Pixel 221967
title('Histogram of Pixel 221967 (High Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');
subplot(2,2,4);
hist(ImageHS(250467,:));%Pixel 250467
title('Histogram of Pixel 250467 (High Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');

% Middle Sensitivity

figureNo = figureNo+1;
figure(figureNo);

subplot(2,2,1);
hist(ImageMS(221667,:)); %Pixel 221667
title('Histogram of Pixel 221667 (Middle Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');
subplot(2,2,2);
hist(ImageMS(250167,:));%Pixel 250167
title('Histogram of Pixel 250167 (Middle Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');
subplot(2,2,3);
hist(ImageMS(221967,:));%Pixel 221967
title('Histogram of Pixel 221967 (Middle Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');
subplot(2,2,4);
hist(ImageMS(250467,:));%Pixel 250467
title('Histogram of Pixel 250467 (Middle Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');

% Low Sensitivity

figureNo = figureNo+1;
figure(figureNo);


subplot(2,2,1);
hist(ImageLS(221667,:));%Pixel 221667
title('Histogram of Pixel 221667 (Low Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');
subplot(2,2,2);
hist(ImageLS(250167,:));%Pixel 250167
title('Histogram of Pixel 250167 (Low Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');
subplot(2,2,3);
hist(ImageLS(221967,:));%Pixel 221967
title('Histogram of Pixel 221967 (Low Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');
subplot(2,2,4);
hist(ImageLS(250467,:));%Pixel 250467
title('Histogram of Pixel 250467 (Low Sensitivity)');
xlabel('Pixel Value');
ylabel('# of Pixels');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 3. Original/Mean/Variance image plots

% High Sensitivity
for j = 1 : size(ImageHS(:,1))
    AvgImHS(j,1) = mean(ImageHS(j,:)); % Calculate mean image for 49 cropped images
    VarImHS(j,1) = var(ImageHS(j,:)); % Calculate variance image for 49 cropped images
end
AvgImageHS = uint16(reshape(AvgImHS, [ ], 950)); % Rearrange the mean image into 950*700 picture
VarImageHS = uint16(reshape(VarImHS, [ ], 950)); % Rearrange the variance image into 950*700 picture

figureNo = figureNo+1;
figure(figureNo);

subplot(2,2,1);
imagesc(AvgImageHS); % Plot mean picture
title('Average Image (High Sensitivity)');
subplot(2,2,2);
imagesc(VarImageHS); % plot variance picture
title('Varience Image (High Sensitivity)');
subplot(2,2,3);
imagesc(uint16(reshape(ImageHS(:,15), [ ], 950))); % plot #15 original picture
title('Original Image 15 (High Sensitivity)');
subplot(2,2,4);
imagesc(uint16(reshape(ImageHS(:,30), [ ], 950))); % plot #30 original picture
title('Original Image 30 (High Sensitivity)');

% Middle Sensitivity
for j = 1 : size(ImageMS(:,1))
    AvgImMS(j,1) = mean(ImageMS(j,:)); % Calculate mean image for 49 cropped images
    VarImMS(j,1) = var(ImageMS(j,:)); % Calculate variance image for 49 cropped images
end
AvgImageMS = uint16(reshape(AvgImMS, [ ], 950)); % Rearrange the mean image into 950*700 picture
VarImageMS = uint16(reshape(VarImMS, [ ], 950)); % Rearrange the variance image into 950*700 picture

figureNo = figureNo+1;
figure(figureNo);

subplot(2,2,1);
imagesc(AvgImageMS); % Plot mean picture
title('Average Image (Middle Sensitivity)');
subplot(2,2,2);
imagesc(VarImageMS); % plot variance picture
title('Varience Image (Middle Sensitivity)');
subplot(2,2,3);
imagesc(uint16(reshape(ImageMS(:,15), [ ], 950))); % plot #15 original picture
title('Original Image 15 (Middle Sensitivity)');
subplot(2,2,4);
imagesc(uint16(reshape(ImageMS(:,30), [ ], 950))); % plot #30 original picture
title('Original Image 30 (Middle Sensitivity)');


% Low Sensitivity
for j = 1 : size(ImageLS(:,1))
    AvgImLS(j,1) = mean(ImageLS(j,:));
    VarImLS(j,1) = var(ImageLS(j,:));
end
AvgImageLS = uint16(reshape(AvgImLS, [ ], 950));
VarImageLS = uint16(reshape(VarImLS, [ ], 950));

figureNo = figureNo+1;
figure(figureNo);

subplot(2,2,1);
imagesc(AvgImageLS); % Plot mean picture
title('Average Image (Low Sensitivity)');
subplot(2,2,2);
imagesc(VarImageLS); % plot variance picture
title('Varience Image (Low Sensitivity)');
subplot(2,2,3);
imagesc(uint16(reshape(ImageLS(:,15), [ ], 950))); % plot #15 original picture
title('Original Image 15 (Low Sensitivity)');
subplot(2,2,4);
imagesc(uint16(reshape(ImageLS(:,30), [ ], 950))); % plot #30 original picture
title('Original Image 30 (Low Sensitivity)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 4. Plot the variance as a function of the mean for each experiment

%High Sensitivy
AvgvarHS = [round(AvgImHS),VarImHS]; % Round the mean values to the nearest integer
AvgvarHS = sortrows(AvgvarHS,1); % Sort the (mean,variance) pairs accroding to ascending mean values
AvguniHS = unique(round(AvgImHS)); % Discard repeated mean values

for k = 1:length(AvguniHS)
    S = find(AvgvarHS(:,1)==AvguniHS(k));
    VaruniHS(k,1) = mean(AvgvarHS(S,2)); % Calculate the average variance for each mean value
end
AvgvaruniHS(:,1) = AvguniHS;
AvgvaruniHS(:,2) = VaruniHS;


for t = 1:(sum((AvgvaruniHS(:,1)<=14000))) % Discard saturated pixels, where pixel value > 14000)
    tAvgvaruniHS(t,:) = AvgvaruniHS(t,:);
end


%tAvgvaruniHS = AvgvaruniHS;

figureNo = figureNo+1;
figure(figureNo);

plot(tAvgvaruniHS(:,1),tAvgvaruniHS(:,2)); % Plot the variance as a function of the mean
title('Mean vs. Variance (High Sensitivity)');
xlabel('Mean Pixel Value');
ylabel ('Variance Value');

%Middle Sensitivy
AvgvarMS = [round(AvgImMS),VarImMS]; % Round the mean values to the nearest integer
AvgvarMS = sortrows(AvgvarMS,1); % Sort the (mean,variance) pairs accroding to ascending mean values
AvguniMS = unique(round(AvgImMS)); % Discard repeated mean values

for k = 1:length(AvguniMS)
    S = find(AvgvarMS(:,1)==AvguniMS(k));
    VaruniMS(k,1) = mean(AvgvarMS(S,2)); % Calculate the average variance for each mean value
end
AvgvaruniMS(:,1) = AvguniMS;
AvgvaruniMS(:,2) = VaruniMS;


for t = 1:(sum((AvgvaruniMS(:,1)<=12500))) % Discard saturated pixels, where pixel value > 12500)
    tAvgvaruniMS(t,:) = AvgvaruniMS(t,:);
end


%tAvgvaruniMS = AvgvaruniMS;

figureNo = figureNo+1;
figure(figureNo);

plot(tAvgvaruniMS(:,1),tAvgvaruniMS(:,2)); % Plot the variance as a function of the mean
title('Mean vs. Variance (Middle Sensitivity)');
xlabel('Mean Pixel Value');
ylabel ('Variance Value');


%Low Sensitivy
AvgvarLS = [round(AvgImLS),VarImLS]; % Round the mean values to the nearest integer
AvgvarLS = sortrows(AvgvarLS,1); % Sort the (mean,variance) pairs accroding to ascending mean values
AvguniLS = unique(round(AvgImLS)); % Discard repeated mean values

for k = 1:length(AvguniLS)
    S = find(AvgvarLS(:,1)==AvguniLS(k));
    VaruniLS(k,1) = mean(AvgvarLS(S,2)); % Calculate the average variance for each mean value 
end
AvgvaruniLS(:,1) = AvguniLS;
AvgvaruniLS(:,2) = VaruniLS;

for t = 1:(sum((AvgvaruniLS(:,1)<=6000))) % Discard saturated pixels, where pixel value > 6000)
    tAvgvaruniLS(t,:) = AvgvaruniLS(t,:);
end

%tAvgvaruniMS = AvgvaruniMS;

figureNo = figureNo+1;
figure(figureNo);

plot(tAvgvaruniLS(:,1),tAvgvaruniLS(:,2)); % Plot the variance as a function of the mean
title('Mean vs. Variance (Low Sensitivity)');
xlabel('Mean Pixel Value');
ylabel ('Variance Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 5. Fit a line to the plotted data

%1st order line fitting

%High Sensiticity
figureNo = figureNo+1;
figure(figureNo);
scatter(tAvgvaruniHS(:,1),tAvgvaruniHS(:,2),'b'); % scatter plot (mean, variance) pairs
hold on;
grid on;
LineHS = polyfit(tAvgvaruniHS(:,1),tAvgvaruniHS(:,2),1); % Linear Line fitting
x = (0:14000);
plot(x,LineHS(1)*x + LineHS(2),'r'); % plot linear regression line
title('Mean vs. Variance (High Sensitivity)');
xlabel('Mean Pixel Value');
ylabel('Variance Value');
legend('Mean vs. Variance Data','Linear Regression');

%Middle Sensiticity
figureNo = figureNo+1;
figure(figureNo);
scatter(tAvgvaruniMS(:,1),tAvgvaruniMS(:,2),'b'); % scatter plot (mean, variance) pairs
hold on;
grid on;
LineMS = polyfit(tAvgvaruniMS(:,1),tAvgvaruniMS(:,2),1); % Linear Line fitting
x = (0:12500);
plot(x,LineMS(1)*x + LineMS(2),'r'); % plot linear regression line
title('Mean vs. Variance (Middle Sensitivity)');
xlabel('Mean Pixel Value');
ylabel('Variance Value');
legend('Mean vs. Variance Data','Linear Regression');


%Low Sensiticity
figureNo = figureNo+1;
figure(figureNo);
scatter(tAvgvaruniLS(:,1),tAvgvaruniLS(:,2),'b'); % scatter plot (mean, variance) pairs
hold on;
grid on;
LineLS = polyfit(tAvgvaruniLS(:,1),tAvgvaruniLS(:,2),1); % Linear Line fitting
x = (0:6000);
plot(x,LineLS(1)*x + LineLS(2),'r'); % plot linear regression line
title('Mean vs. Variance (Low Sensitivity)');
xlabel('Mean Pixel Value');
ylabel('Variance Value');
legend('Mean vs. Variance Data','Linear Regression');

%{
%2nd order line fitting 
%High Sensiticity
figureNo = figureNo+1;
figure(figureNo);
scatter(tAvgvaruniHS(:,1),tAvgvaruniHS(:,2),'b'); % scatter plot (mean, variance) pairs
hold on;
grid on;
LineHS2 = polyfit(tAvgvaruniHS(:,1),tAvgvaruniHS(:,2),2); % Linear Line fitting
x = (224:14000);
plot(x,LineHS2(1)*(x.^2) + LineHS2(2)*x + LineHS2(3),'r'); % plot linear regression line
title('Mean vs. Variance (High Sensitivity)');
xlabel('Mean Pixel Value');
ylabel('Variance Value');
legend('Mean vs. Variance Data','Quadratic Regression');


%Low Sensiticity
figureNo = figureNo+1;
figure(figureNo);
scatter(tAvgvaruniLS(:,1),tAvgvaruniLS(:,2),'b'); % scatter plot (mean, variance) pairs
hold on;
grid on;
LineLS2 = polyfit(tAvgvaruniLS(:,1),tAvgvaruniLS(:,2),2); % Linear Line fitting
x = (224:14000);
plot(x,LineLS2(1)*(x.^2) + LineLS2(2)*x + LineLS2(3),'r'); % plot linear regression line
title('Mean vs. Variance (Low Sensitivity)');
xlabel('Mean Pixel Value');
ylabel('Variance Value');
legend('Mean vs. Variance Data','Quadratic Regression');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 6. Use the fitted line data to calculate the camera gain g, the read noise variance (measured in photo-electrons), and the ADC noise variance (measured in gray levels)
LineHS(1)
LineHS(2)
LineLS(1)
LineLS(2)
NoiseRead = (LineMS(2) - LineLS(2))/(LineMS(1)^2 - LineLS(1)^2)
NoiseADC = LineMS(2) - (LineMS(1)^2) * NoiseRead

%{
ans =

   41.6896


ans =

  -9.7477e+04


ans =

    6.5435


ans =

  -6.5114e+03


NoiseRead =

  -53.6605


NoiseADC =

  -4.2138e+03
  
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 7. Plot the SNR as a function of mean pixel value

%High Sensiticity

SNRHS = 20*log10(tAvgvaruniHS(:,1)./sqrt(tAvgvaruniHS(:,2))); % Calculate SNR in logarithmic domain

figureNo = figureNo+1;
figure(figureNo);

plot(tAvgvaruniHS(:,1),SNRHS); % Plot the SNR as a function of mean pixel value
title('Mean vs. Signal-to-Noise Ratio (High Sensitivity)');
xlabel('Mean Pixel Value');
ylabel ('Signal-to-Noise Ratio (dB)')

%Middle Sensiticity

SNRMS = 20*log10(tAvgvaruniMS(:,1)./sqrt(tAvgvaruniMS(:,2))); % Calculate SNR in logarithmic domain

figureNo = figureNo+1;
figure(figureNo);

plot(tAvgvaruniMS(:,1),SNRMS); % Plot the SNR as a function of mean pixel value
title('Mean vs. Signal-to-Noise Ratio (Middle Sensitivity)');
xlabel('Mean Pixel Value');
ylabel ('Signal-to-Noise Ratio (dB)')

%Low Sensiticity
SNRLS = 20*log10(tAvgvaruniLS(:,1)./sqrt(tAvgvaruniLS(:,2))); % Calculate SNR in logarithmic domain

figureNo = figureNo+1;
figure(figureNo);

plot(tAvgvaruniLS(:,1),SNRLS); % Plot the SNR as a function of mean pixel value
title('Mean vs. Signal-to-Noise Ratio (Low Sensitivity)');
xlabel('Mean Pixel Value');
ylabel ('Signal-to-Noise Ratio (dB)')



