function [fwhm, hhx, hhy, fwhma]=hhfw(varargin)
%Usage: [fwhm, hhx, hhy, fwhma]=hhfw(ydata, xdata)
%
% ydata:    the data for which HHFW will be calculated
% xdata:    the x-axis (optional); if not provided, xdata=1:length(y)
%
% fwhm:     full width at half max
% hhx:      1 x 2, x-axis values for half-height
% hhy:      1 x 2, repeats of half max value of y data
% fwhma:    hhfw integral area
%
debugflag=0;

ydata=varargin{1};

if length(varargin)>1
    xdata=varargin{2};
else
    xdata=1:length(ydata);
end
dx=xdata(2)-xdata(1);

%find max in this window:
maxy    = max(ydata);
hhy     = [1 1]*maxy/2; % easier to plot, later
maxi    = find(ydata == maxy);
[~,idx] = min(abs(maxi - ceil(numel(ydata)/2)));
maxi    = maxi(idx); % take max closest to center

%find left half-height:
li=maxi;
while (ydata(li)>hhy(1))&&(li>1)
    li=li-1;
end
m=(ydata(li+1)-ydata(li))/(xdata(li+1)-xdata(li));
hhx(1)=(hhy(1)-ydata(li))/m+xdata(li);
if li==1
    warning('Error! No left half value found.')
end

%find right half-height
ri=maxi;
while (ydata(ri)>hhy(2))&&(ri<length(xdata))
    ri=ri+1;
end
m=(ydata(ri)-ydata(ri-1))/(xdata(ri)-xdata(ri-1));
hhx(2)=(hhy(2)-ydata(ri-1))/m+xdata(ri-1);
if ri==length(xdata)
    warning('Error! No right half value found.')
end

%Calculate full width at half max:
fwhm=hhx(2)-hhx(1);

%[li maxi ri]
%Calculate HHFW area:
if (li+1)<(ri-1)
    fwhma=trapz(xdata(li+1:ri-1),ydata(li+1:ri-1));
else
    fwhma=0;
end
%Add partial index at left:
tmp = 0.5*(hhy(1)+ydata(li+1))*(xdata(li+1)-hhx(1));
fwhma=fwhma+tmp;
%Add partial index at right:
tmp = 0.5*(ydata(ri-1)+hhy(2))*(hhx(2)-xdata(ri-1));
fwhma=fwhma+tmp;

if debugflag
    figure(999)
    plot(xdata,ydata,'k-',hhx,hhy,'ro')
    title(sprintf('HHFW Area: %f',fwhma))
end
