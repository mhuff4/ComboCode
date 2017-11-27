%% Nonlinear sweep correction routine
% use peaks (pk) obtained from fidupks3w routine
function [img3i,deltat_new] = timeaxiscorr1(img3,pk,show,ncomb)
ypk=diff(pk); %differences between comb peak locations
xpk=pk(1:end-1)+ypk/2; %location of midway points between comb peaks
oldx=1:size(img3,2);

if ncomb == 0.5000
    dtpx = 4*0.506087; 
elseif ncomb == 1
    dtpx=2*0.506087;
elseif ncomb == 0
    dtpx=0.548;
elseif ncomb==2;
    dtpx=0.506087;
end

% yfit is px spacing per picket interval for every pixel
% [yfit2] = pchip(xpk,ypk,oldx); 
% yfit=yfit2';
c=fit(xpk',ypk','smoothingspline');
yfit=c(oldx);
yfit(1:(floor(pk(1))-1))=yfit(floor(pk(1)));
yfit((floor(pk(end))+1):end)=yfit(floor(pk(end)));

if show 
    figure,plot(oldx,yfit,'-g',xpk,ypk,'.r',xpk,mean(ypk)*ones(size(ypk)),':k');
    xlim([min(oldx) max(oldx)]);
    title('Time axis correction')
    xlabel('Time (px)')
    ylabel('px spacing per picket interval')
end
oldtimeax=cumsum(dtpx./yfit);  %0.506./yfit gives delta-ps for the given pixel. sum gives time axis.
oldtimeax=oldtimeax-oldtimeax(1); %makes time start at 0
deltat_new=mean(dtpx./ypk); %CHANGED FROM YFIT 2/3/16 average time/pixel
newtimeax=(0:size(img3,2)-1)*deltat_new; %linearized new time axis using average deltat_new as slope

img3a=interp1(oldtimeax,img3',newtimeax,'pchip')';%interpolates to find img3a, the values of the underlying function img3' at
%the points in newtimeax. 
if ncomb==1
    img3i=img3a;
else
%adjust intensity by ratio of dwell times
Iratio=yfit/mean(ypk)'; %CHANGED FROM YFIT 2/3/16 %larger spacing between peaks=camera spent less time on those pixels -> intensity needs to increase
img3i=img3a.*repmat(Iratio,1,size(img3a,1))';
% figure,plot(oldtimeax,mean(img3(350:450,:)),'b',newtimeax,mean(img3a(350:450,:)).*Iratio','r',newtimeax,mean(img3i(350:450,:)),'g')
% img3i=img3a;
end

end




