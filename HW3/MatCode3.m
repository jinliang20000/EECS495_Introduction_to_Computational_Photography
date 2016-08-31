% EECS 395/495: Introduction to Computational Photography
% HW3: Flash/No Flash Photography

% Student Name: Haikun Liu
% Student Number: 2903021
% NetID: hlg483

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('/Users/HKLHK/Documents/MATLAB/CP HW3/');

%Load images
[Im_w,map] = imread('/Users/HKLHK/Desktop/2015 Fall Quarter (NU)/EECS495 Introduction to Computational Photography/HW/HW3/Images/white.jpg');
Im_w = double(Im_w);
whiteIm = imcrop(Im_w,map,[95 145 1023 767]);%crop the photo center at (1400,700) with size 1024*768

[Im_b,map] = imread('/Users/HKLHK/Desktop/2015 Fall Quarter (NU)/EECS495 Introduction to Computational Photography/HW/HW3/Images/black.jpg');
Im_b = double(Im_b);
blackIm = imcrop(Im_b,map,[95 145 1023 767]);%crop the photo center at (1400,700) with size 1024*768


%Finding sigma_s and sigma_r for no-flash image

sigma_s_t = [1 1 1 4 4 4 16 16 16];
sigma_r_t = [0.05 0.10 0.20 0.05 0.10 0.20 0.05 0.10 0.20];
figure;
for j=1:9

for i=1:3
    blackIm_Denoise_t(:,:,i) = bilateralFilter(blackIm(:,:,i),sigma_s_t(j),sigma_r_t(j)*max(max(blackIm(:,:,i))));%denoise each color channel separately
end
subplot(3,3,j);
image(uint8(blackIm_Denoise_t));
xlabel(['sigma_s = ' num2str(sigma_s_t(j)) ', sigma_r =' num2str(sigma_r_t(j))]);
end

%choosing sigma_s = 16 and sigma_r = 0.1 for no-flase image

sigma_s_b = 16;
sigma_r_b = 0.1;

for i=1:3
    blackIm_Denoise(:,:,i) = bilateralFilter(Im_b(:,:,i),sigma_s_b,sigma_r_b*max(max(Im_b(:,:,i))));
end



%denoise parameter testing for flash image
sigma_s_t = [1 1 1 4 4 4 16 16 16];
sigma_r_t = [0.05 0.10 0.20 0.05 0.10 0.20 0.05 0.10 0.20];
figure;
for j=1:9

for i=1:3
    whiteIm_Denoise_t(:,:,i) = bilateralFilter(whiteIm(:,:,i),sigma_s_t(j),sigma_r_t(j)*max(max(whiteIm(:,:,i))));
end
subplot(3,3,j);
image(uint8(whiteIm_Denoise_t));
xlabel(['sigma_s = ' num2str(sigma_s_t(j)) ', sigma_r =' num2str(sigma_r_t(j))]);
end

%choosing sigma_s = 16 and sigma_r = 0.1 for flase image

sigma_s_w = 16;
sigma_r_w = 0.1;

for i=1:3
    whiteIm_Denoise(:,:,i) = bilateralFilter(Im_w(:,:,i),sigma_s_w,sigma_r_w*max(max(Im_w(:,:,i))));
end


%fuse image

epsilon = 0.02;
fusedIm = blackIm_Denoise.*((Im_w+epsilon)./(whiteIm_Denoise+epsilon));

%compare plots
figure;
image(uint8(fusedIm));
figure;
image(uint8(blackIm_Denoise))
figure;
image(uint8(Im_b));




