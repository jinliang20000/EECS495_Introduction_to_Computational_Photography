clear;
clc;
figureNo = 1;
NumOfIms = 5;
str_Path = '/Users/HKLHK/Desktop/Northwestern University Q1/EECS495 Introduction to Computational Photography/HW/HW5/Image/';

for i = 1: NumOfIms
    str_Load = strcat(str_Path, num2str(i), '.jpg');
    Image = imread(str_Load);
    I(:,:,:,i) = Image;
end

f = 2.95/1000;
v=[1.5 2.25 3.375 5.0625 7.59375];
%v=flip([1.5 2.25 3.375 5.0625 7.59375]);
u = 1./(1/f-1./v);

m = u(end)./u;
mr = m;
mi = 2-m;

[row column color numIm] = size(I);

for k=1:numIm
    for i=1:row
        for j = 1:column
            
            if i<=row/2
                mx = round(i*mi(k));
                if(mx<=0)
                    mx=1;
                end
            else
                mx = round(i*mr(k));
                if(mx>row)
                    mx=row;
                end
            end
            
            if j<=column/2
                my = round(j*mi(k));
                if(my<=0)
                    my=1;
                end
            else
                my = round(j*mr(k));
                if(my>column)
                    my=column;
                end
            end
  
           I_apostrophe(i,j,:,k) = I(mx,my,:,k); 
        end
    end
end

for k = 1: numIm
    I_gray(:,:,k) = rgb2gray(I_apostrophe(:,:,:,k));
end

Laplacian = [1,4,1;4,-20,4;1,4,1]/6;

for k=1:numIm
    LaplacianImage = imfilter(I_gray(:,:,k),Laplacian,'replicate');
    M(:,:,k) = LaplacianImage.*LaplacianImage;
end

% kernelSize = 5;
% 
% for k=1:numIm
%     M(:,:,k) = squaredLaplacian(I_gray(:,:,k),kernelSize);
% end

for i=1:row
    for j = 1:column
       Depthmap(i,j)=max(M(i,j,:));
    end
end

figure;
imshow(Depthmap);

for i=1:row
    for j = 1:column
       index = find(M(i,j,:)==max(M(i,j,:)));
       if length(index)>1
           D(i,j) = index(1);
       else
           D(i,j) = index;
       end
    end
end

DepthIndexMap = uint8(D.*(255/max(max(D))));
figure;
imshow(DepthIndexMap);

for c = 1:color
    for i=1:row
        for j = 1:column
            A(i,j,c) = I_apostrophe(i,j,c,D(i,j));
        end
    end
end

% A = double(A);
% 
% for c=1:color
%     A(:,:,c) = bilateralFilter(A(:,:,c),1.5,0.1*max(max(A(:,:,c))));
% end
% A = uint8(A);
figure;
imshow(A)


