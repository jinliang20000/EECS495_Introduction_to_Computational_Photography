function [ xoffSet,yoffSet ] = correlation( template,searchwindow )
c = normxcorr2(template,searchwindow);
[ypeak, xpeak] = find(c==max(c(:)));
yoffSet = ypeak-size(template,1);
xoffSet = xpeak-size(template,2);
end

