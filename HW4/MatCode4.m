%%%%%%
clear;
clc;
str_Path = '/Users/HKLHK/Desktop/2015 Fall Quarter (NU)/EECS495 Introduction to Computational Photography/HW/HW4/Image/';
prevExposure = 67200;
j=1;
while prevExposure <= 275251200
    str_Load = strcat(str_Path, num2str(prevExposure), '.jpg');
    Image = imread(str_Load);
    Info = imfinfo(str_Load);
    I(:,:,j) = double(reshape(Image, [ ], 3));
    eT(j) = Info.DigitalCamera.ExposureTime;
    prevExposure = prevExposure + prevExposure;
    j = j + 1;
end

%%%%%

[numPi color numIm] = size(I);
j=1;
for i=1:numIm
    if j<=5
       if (sum(sum(all(I(:,:,(numIm-i+1))==255,2)))/numPi)<=0.2
          Im(:,:,5-j+1) =  I(:,:,(numIm-i+1));
          exposureTime(5-j+1) = eT(numIm-i+1);
          j = j+1;
       end
    else
        clear I eT;
        break;
    end
end

%%%%% find the camera response curves for the shield tablet 
numPrt = 1000;
id = randperm(numPi,numPrt);
Z_red = zeros(numPrt,length(exposureTime));
Z_green = zeros(numPrt,length(exposureTime));
Z_blue = zeros(numPrt,length(exposureTime));
for i=1:length(exposureTime)
    Z_red(:,i)= Im(id,1,i);
    Z_green(:,i) = Im(id,2,i);
    Z_blue(:,i) = Im(id,3,i);
end

l = 5; % [.1, 5]

[g_red,lE_red]=gsolve(Z_red,log(exposureTime),l);
[g_green,lE_green]=gsolve(Z_green,log(exposureTime),l);
[g_blue,lE_blue]=gsolve(Z_blue,log(exposureTime),l);


for i=1:length(exposureTime)
    X_red(:,i)= lE_red + log(exposureTime(i));
    X_green(:,i) = lE_green + log(exposureTime(i));
    X_blue(:,i) = lE_blue + log(exposureTime(i));
end

figure;
for i=1:length(exposureTime)
    scatter(X_red(:,i),Z_red(:,i),'.','b');
    hold on;
end
plot(g_red,0:255,'r','linewidth',2);
title('Response curve for red channel');
xlabel('log exposure');
ylabel('pixel value(Z)')
legend('fitted data','response curve')
hold off;

figure;
for i=1:length(exposureTime)
    scatter(X_green(:,i),Z_green(:,i),'.','b');
    hold on;
end
plot(g_green,0:255,'r','linewidth',2);
title('Response curve for green channel');
xlabel('log exposure');
ylabel('pixel value(Z)')
legend('fitted data','response curve')
hold off;

figure;
for i=1:length(exposureTime)
    scatter(X_blue(:,i),Z_blue(:,i),'.','b');
    hold on;
end
plot(g_blue,0:255,'r','linewidth',2);
title('Response curve for blue channel');
xlabel('log exposure');
ylabel('pixel value(Z)')
legend('fitted data','response curve')
hold off;

figure;
plot(g_red,0:255,'r','linewidth',2);
hold on;
plot(g_green,0:255,'g','linewidth',2);
hold on;
plot(g_blue,0:255,'b','linewidth',2);
title('Response curves');
xlabel('log exposure');
ylabel('pixel value(Z)')
legend('Red','Green','Blue')
hold off;

%%%%%Recover the HDR radiance map of the scene
[numPi color numIm] = size(Im);
lnE = zeros(numPi,3);
for i=1:numPi
    g(:,1) = g_red(Im(i,1,:)+1);
    lnE(i,1) = sum(g-log(exposureTime)')/numIm;
    
    g(:,1) = g_green(Im(i,2,:)+1);
    lnE(i,2) = sum(g-log(exposureTime)')/numIm;
    
    g(:,1) = g_blue(Im(i,3,:)+1);
    lnE(i,3) = sum(g-log(exposureTime)')/numIm;
end

re_lnE = zeros(1944,2592,3);
for i=1:color
    re_lnE(:,:,i) = reshape(lnE(:,i), [], 2592);
end
[X,Y] = meshgrid(1:2592,1:1944);

figure;
h=pcolor(X,Y,flip(re_lnE(:,:,1)));
set(h,'edgecolor','none','facecolor','interp');
colorbar;
title('Recovered radiance map: Red (ln scale)');
figure;
h=pcolor(X,Y,flip(re_lnE(:,:,2)));
set(h,'edgecolor','none','facecolor','interp');
colorbar;
title('Recovered radiance map: Green (ln scale)');
figure;
h=pcolor(X,Y,flip(re_lnE(:,:,3)));
set(h,'edgecolor','none','facecolor','interp');
colorbar;
title('Recovered radiance map: Blue (ln scale)');


E = exp(lnE);

re_E = zeros(1944,2592,3);
for i=1:color
    re_E(:,:,i) = reshape(E(:,i), [], 2592);
end
figure;
imshow(uint8(re_E));


%%%%%Implement a tone mapping algorithm to display HDR image 

for c=1:color
    E_norm(:,c) = (E(:,c) - min(E(:,c)))./(max(E(:,c)) - min(E(:,c)));
end

re_E_norm = zeros(1944,2592,3);
for i=1:color
    re_E_norm(:,:,i) = reshape(256*E_norm(:,i)-1, [], 2592);
end

figure;
subplot(1,2,1);
imshow(uint8(re_E_norm));
subplot(1,2,2);
hist(E_norm);
legend('r','g','b');

gamma = 0.2;
E_ga = E_norm.^gamma;

for c=1:color
    E_gamma(:,c) = (E_ga(:,c) - min(E_ga(:,c)))./(max(E_ga(:,c)) - min(E_ga(:,c)));
end

re_E_gamma = zeros(1944,2592,3);
for i=1:color
    re_E_gamma(:,:,i) = reshape(256*E_gamma(:,i)-1, [], 2592);
end

figure;
subplot(1,2,1);
imshow(uint8(re_E_gamma));
subplot(1,2,2);
hist(E_gamma);
legend('r','g','b');

L = rgb2gray(E_gamma);

L_avg = exp(mean(mean(log(L))));

a= 0.7;

T = a/L_avg *L;


L_tone = T.*(1+T./((max(max(T)))^2))./(1+T);

M = L_tone./L;

RGB_new = (265*E_gamma-1).*M;

re_RGB_new = zeros(1944,2592,3);
for i=1:color
    re_RGB_new(:,:,i) = reshape(RGB_new(:,i), [], 2592);
end

re_RGB_new = uint8(re_RGB_new);

figure;
imshow(re_RGB_new);