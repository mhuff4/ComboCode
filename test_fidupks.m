function [pk1,xx,deltat,fidu,fidutr,comb_corr] = test_fidupks(shotnumber,n,fidu,del,kind,show,start)
%%
fidutr=sum(fidu); 
fidutr=fidutr-min(fidutr); fidutr=fidutr/max(fidutr);
figure,plot(fidutr)
%sums columns to get peak profile. max peak value = 1 min peak value = 0.
npts=length(fidutr);
%tot number of px, sets cutoff px value for last peak
[pks,locs] = findpeaks(fidutr, 'minpeakheight', .4, 'minpeakdistance',15); %.5 and 10
%approx peak locations
pinc = round(mean(diff(locs))); %approx diff between peaks
%p1 = find(fidutr>=0.5,1,'First'); %returns first index corresponding to 0.5
ntines = length(locs);
if kind == 'a'
    if n == 0.5000
        dtpx = 4*0.506087; 
        twidth=0.32*pinc/dtpx;
    elseif n==8
        dtpx=1/8;
        twidth=0.007*pinc/dtpx;
    elseif n==1
        dtpx=2*0.506087;
        twidth=0.1*pinc/dtpx;
    elseif n==0
        dtpx = 0.548; twidth=0.22*pinc;
    else
        dtpx=0.506087;
        %twidth=0.08*pinc/dtpx; %usual
        twidth=0.035*pinc/dtpx;
    end
elseif kind == 'g'
    dtpx=800; %(um between geometric peaks )
    twidth=0.08*pinc;
else dtpx = 0.548; twidth=0.22*pinc;
end

% if show
%     figure, plot(fidutr);
%     hold on
%     plot(locs,pks,'or');
% end
nflag = 0;
if (npts/2-floor(npts/2))
    npts=npts-1;  fidutr=fidutr(1:end-1); nflag=1; 
    %makes even # of columns, ideal for FFTshift
end 
fidutr(1:5)=0;  fidutr(end-5:end)=0; 
%sets first 10 and last 10 columns to 0
xx=1:npts;
comb=zeros(1,npts);

%%
dw=2*pi/npts;
wArray=linspace(-npts/2,npts/2-1,npts)*dw; 
%discretized fourier exponent centered at 0 ranging from -pi to pi
fiduFFT=fftshift(ifft(fftshift(fidutr))); 
%decenters fidutr in prep for ifft (assumes origin at (1,1)) recenters
%image after ifft. ifft conatins the 1/N coefficient.
comb_corr=zeros(1,npts);
clear pk1
pk1(1)=locs(1); %assumed value of first peak. 
comb=exp(-((xx-pk1(1))/twidth).^2); %kk=0; sets first gaussian peak at pk1
for kk=1:ntines
    combFFT=fftshift(ifft(fftshift(comb)));
    DelayTemp=fminsearch(@(Delay) -real(sum(conj(fiduFFT).*combFFT.*exp(1i*wArray*Delay))),0);
    %cross correlation to find time delay between comb peak and kkth peak
    %in fidutr
    combFFT=combFFT.*exp(1i*wArray*DelayTemp); %shifts comb to correct position
    pk1(kk)=pk1(kk)+DelayTemp; %gives correct peak placement
    comb_corr=comb_corr+real(fftshift(fft(fftshift(combFFT))))*fidutr(round(pk1(kk))); 
    %peak multiplied by its value in fidutr to give same height
    if show, clf; plot(xx,fidutr,'k',xx,comb_corr,'r'); pause(.1);end
    if kk>1, disp(num2str(round([kk,pk1(kk-1),pk1(kk),pk1(kk)-pk1(kk-1)]))) 
             %peak #, prev peak loc, peak loc, peak diff
        if pk1(kk)<(pk1(kk-1)+5)  
           disp([num2str(kk), '  ' num2str(round(pk1(kk)))])
           kk=kk-1; pk1=pk1(1:end-1); break
        elseif pk1(kk)+pinc>npts  
           disp([num2str(kk), '  ' num2str(round(pk1(kk)))]); break
        elseif kk==ntines, break
        end
    end
    pk1(kk+1)=pk1(kk)+pinc; %estimated placement of next peak
    comb=exp(-((xx-pk1(kk+1))/twidth).^2); %redefines comb for next peak
end
%%
pk1=pk1+start-2; %remove if not using start
c_spacg=diff(pk1);
meansep=mean(c_spacg);
rms=100*std(c_spacg)/meansep;
deltat=dtpx/meansep; %average time/px
ffidu=pk1(1)*del;
if kind == 'a'
str=['avg peak separation: ' num2str(meansep,'%4.2f') ' px = ' num2str(dtpx,'%4.3f') ' ns'];
elseif kind =='g'
    str=['avg peak separation: ' num2str(meansep,'%4.4f') ' px = ' num2str(dtpx,'%4.4f') ' um'];
    str2=['magnification: ' num2str(27*meansep/dtpx)];
    disp(str2)
else
str=['avg peak separation: ' num2str(meansep,'%4.2f') ' px = ' num2str(dtpx,'%4.3f') ' ns'];
str3=['time of first fidu peak:  ' num2str(ffidu,'%4.4f')];
disp(' '); disp(str3)
end
disp(' ');  disp(str)
str1=['time axis: ' num2str(deltat*1000,'%4.2f') ' ps/px  for shot ',num2str(shotnumber)];
disp(str1)
str2=['rms time axis nonlinearity: ' num2str(rms,'%5.2f') ' % rms'] ;
disp(str2)


%%
if nflag, fidutr=[fidutr,0]; comb_corr=[comb_corr,0]; xx=[xx,xx(end)+1]; end
%comb_corr = new cleaner comb fit
