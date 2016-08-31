clear;
clc;
str_Load = '/Users/HKLHK/Desktop/Northwestern University Q1/EECS495 Introduction to Computational Photography/HW/HW6/VID_20151127_160303.mp4';
v = VideoReader(str_Load);
i=1;
video = read(v);
[row column color numIm] = size(video);
sampled_video=video(:,:,:,1:round(numIm/60):end);
[row column color numIm] = size(sampled_video);
for k = 1: numIm
    gray_video(:,:,k) = rgb2gray(sampled_video(:,:,:,k));
end

% for i=1:row
%     for j=1:column
%         var_gray_video(i,j)=var(double(gray_video(i,j,:)));
%     end
% end
% 
% normalized_var_gray_video = uint8(round(var_gray_video/max(max(var_gray_video))*255));


for i=1:row
    for j=1:column
        var_gray_video(i,j)=sum(double(gray_video(i,j,:)));
    end
end
normalized_var_gray_video = uint8(round(var_gray_video/max(max(var_gray_video))*255));

figure;
[template, rect_template] = imcrop(gray_video(:,:,1));

[search, rect_search] = imcrop(normalized_var_gray_video);

clear search;

for k = 1: numIm
    search(:,:,k) = imcrop(gray_video(:,:,k),rect_search);
end 

XWindowPixelShift = [];
YWindowPixelShift = [];
for k = 1: numIm
    [xoffSet,yoffSet] = correlation( template(:,:),search(:,:,k) );
    XWindowPixelShift(k) = xoffSet;
    YWindowPixelShift(k) = yoffSet;
end 

XPixShift=XWindowPixelShift-XWindowPixelShift(1);
YPixShift=YWindowPixelShift-YWindowPixelShift(1);
figure;
scatter(XPixShift,YPixShift,'*');
title('Pixel shift for each frame of video');
ylabel('Y Pixel Shift');
xlabel('X Pixel Shift');
axis('square');
for k = 1: numIm
    T(:,:,:,k) = imtranslate(sampled_video(:,:,:,k),[-XPixShift(k), -YPixShift(k)],'FillValues',255);
end

P = uint8(round(sum(T,4)/numIm));
Crop_P = imcrop(P,[max(abs(YPixShift)) max(abs(XPixShift)) column row]);
figure;
imshow(Crop_P);

%Finished%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following code is for Part 4%

figure;
[template, rect_template] = imcrop(gray_video(:,:,1));

[search, rect_search] = imcrop(normalized_var_gray_video);

clear search;

for k = 1: numIm
    search(:,:,k) = imcrop(gray_video(:,:,k),rect_search);
end 

XWindowPixelShift = [];
YWindowPixelShift = [];
for k = 1: numIm
    [xoffSet,yoffSet] = correlation( template(:,:),search(:,:,k) );
    XWindowPixelShift(k) = xoffSet;
    YWindowPixelShift(k) = yoffSet;
end 

XPixShift=XWindowPixelShift-XWindowPixelShift(1);
YPixShift=YWindowPixelShift-YWindowPixelShift(1);
for k = 1: numIm
    T(:,:,:,k) = imtranslate(sampled_video(:,:,:,k),[-XPixShift(k), -YPixShift(k)],'FillValues',255);
end

P = uint8(round(sum(T,4)/numIm));
Crop_P = imcrop(P,[max(abs(YPixShift)) max(abs(XPixShift)) column row]);
figure;
imshow(Crop_P);



