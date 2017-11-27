%% Isotemporal Analysis for Emission Peak Images on Streak Cameras
clear all
close all
%% Open the hdf image
stordir = ['/Users/bhend/Desktop/test dump/'] %directory to dump images

shotnum = '25386';
filename = strcat('/Users/mhuf/Documents/GitHub/ComboCode/sop_s25386.hdf'); %path to silicon SOP hdf
shot = readhdf(filename); %Reads the hdf file
streak = shot.Streak_array.Streak_array; %branch which contains SOP streak image
delt = shot.Streak_array.deltat_per_px_in_ns; %ratio of nanoseconds per pixel
t0 = shot.Streak_array.t0_of_first_px_in_ns; %t=0 as defined by calibration
image = imrotate(streak,90); %rotate image by 90 degrees
image = medfilt2(image, [2 2]); %Filter SOP streak image by taking the median in 2x2 squares and replacing the values with these medians
maxsig = 30*median(median(image)); %New intensity threshold maximum set by the overall median intensity of the image

figure(1)
imagesc(image,[0,maxsig]) %filtered SOP image with new pixel intensity threshold
%% Lineout
win = figure; imagesc(image,[0,maxsig]); title('Select ROI')
box = imrect(gca, [100, 800, 200, 200]) %selecting region of interest (gca is last figure generated)
pos = wait(box); %collects [xmin ymin width height] wrt top left corner
xmin = round(pos(1,1)); xmax = round(pos(1,1)+pos(1,3)); ymin = round(pos(1,2)); ymax = round(pos(1,2)+pos(1,4)); %New axis parameters as defined by ROI box
close(win)

%% Plot
yaxis = (ymin:ymax); %Create new x and y axes defined by the ROI
xaxis = (xmin:xmax);
xaxis = xaxis.'; %convert row vector to column vector
len = length(yaxis);
roi = image(ymin:ymax,xmin:xmax);

avglineout = mean(roi,1); %Take the mean of the lineouts vertically
avglineout = avglineout.'; %convert row vector to column vector
[amp, peak] = max(avglineout); %Finding the center of the emission peak
region = (xmin:peak+20); %create region for fitting
peakregion = avglineout(region);

%Convert pixels to time
timeaxis = (xaxis - 1)*delt + t0; %First pixel =1 so we need to shift it to 0 before converting to time

figure(2)
subplot(1,2,1); imagesc(timeaxis,yaxis,roi); xlabel('Time (nanoseconds)'); ylabel('Target Position (pixels)'); title('Region of Interest')
subplot(1,2,2); plot(timeaxis,avglineout,'b'); xlabel('Time (nanoseconds)'); ylabel('Counts'); title('Position Averaged Intensity Lineout'); xlim([min(timeaxis) max(timeaxis)]);
%% Fitting to Emission Peak
range = xmin:xmax;
lorentz = @(A,x0,B,C,x) A./((x-x0).^2 + B) + C;
ft1 = fittype(lorentz,'independent','x');
guess = [.4, 2.5, .4, 0];
options = fitoptions('method','nonlinearleastsquares','startpoint',guess,'robust','on');
f = fit(timeaxis(region),avglineout(region),ft1,options);

figure(3)
plot(timeaxis,avglineout,'b',timeaxis,f(timeaxis),'r'); xlabel('Time (nanoseconds)'); ylabel('Counts'); title('Fitted Lineout'); legend(strcat('t_{peak} = ',num2str(f.x0),' ns')); xlim([min(timeaxis) max(timeaxis)]);