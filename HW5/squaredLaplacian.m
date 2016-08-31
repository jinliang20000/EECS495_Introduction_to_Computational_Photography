function [ focusMeasure ] = squaredLaplacian(Image,kernelSize)
k=kernelSize;
[M,N] = size(Image);
DH = Image;
DV = Image;
DH(1:M-k,:) = diff(Image,k,1);
DV(:,1:N-k) = diff(Image,k,2);
FM = max(DH, DV);        
focusMeasure = FM.*FM;
end

