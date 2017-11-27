%% Time Calibration SOP
close all
clear all
%corrects time dimension of SOP images
%% Shot Parameters
shotnumber='25386';
ep=1;
thres=5000;
geocorr=0;
dt_offset=0;%ns 
% stordir='T:\Documents and Settings\mgreg\Desktop\IM15B\Images\'; %corrected image locatio
%stordir='\\hopi-fs\shares\Experimental\SSP_exps\SSP_EOS\Michelle\IM15B\Images\';
stordir=strcat('/Users/bhend/Desktop/Silicon EOS Data/April 2017/s',shotnumber,'/'); %corrected image locatio
%% Import Data
if ep==0
shotname = [strcat('/Users/bhend/Desktop/Silicon EOS Data/April 2017/s',shotnumber,'/sop_s',shotnumber,'.hdf')];
else
shotname = [strcat('/Users/bhend/Desktop/Silicon EOS Data/April 2017/s',shotnumber,'/sop_s',shotnumber,'.hdf')];
end
%shotname = ['T:\Documents and Settings\mgreg\Desktop\May1asbo\sop_sop-ross_s',shotnumber,'.hdf'];
%shotname=['\\redwood\sop-ross\Calibration\Spectral Response\2013_03_27_calibration\SOPsrTimebase_BP0nm_ODcode0_032713_072948_v0Raw.hdf'];
%shotname=['\\redwood\sop-ross\Calibration\Spectral Response\2013_03_27_calibration\SOPsrTimebase_BP0nm_ODcode0_032713_073049_v1Raw.hdf'];

shot=readhdf(shotname);
im=imrotate(squeeze(double(shot.Streak_array.Streak_array(:,:,1)))-squeeze(double(shot.Streak_array.Streak_array(:,:,2))),90);%shot
%im=imrotate(double(shot.fore.fore)-double(shot.back.back),90);%timebase
%im=double(shot.fore.fore)-double(shot.back.back);%spatial mag

if geocorr==1
if ep==0
%     image=image(:,30:1050);
    %refshotnum=20130327; %extracts reference shot #
    %distmssg=[{'SOP slow sweep spatial correction'};{['taken from date: ' num2str(refshotnum)]}];
    refshotnum=20140903; %extracts reference shot #
    distmssg=[{'SOP fast sweep spatial correction'};{['taken from date: ' num2str(refshotnum)]}];
    if shot.Streak_array.SweepSpeed==3    
        load('geoCal_SOP_17ns_090314.mat','TFORM')
    elseif shot.Streak_array.SweepSpeed==4 
        load('geoCal_SOP_6ns_090314.mat','TFORM')
    elseif shot.Streak_array.SweepSpeed==2 
        load('geoCal_SOP_46ns_090314.mat','TFORM')
    elseif shot.Streak_array.SweepSpeed==1 
        load('geoCal_SOP_96ns_090314.mat','TFORM')
    end
    load('geoCal_SOP_slow_032713.mat','TFORM')
    im=imtransform(im,TFORM,'xdata',[1 1100],'ydata',[1 1100]); 
    display(distmssg)
else
    im=im(140:1050,30:1050);
end
else
    refshotnum='1';
    im=im(1:1100,20:1070);
end
%% Fiducial analysis and time axis linearization
ncomb=shot.Streak_array.OCM_CombFrequency;%determines which comb was used
figure, imagesc(im,[0 thres]);
h = imrect(gca, [18 71 900 50]);
pos = wait(h);
xmin=round(pos(1,1)); xmax=round(xmin+pos(1,3)); ymin=round(pos(1,2)); ymax=round(ymin+pos(1,4));

        fidu1=im(ymin:ymax,xmin:xmax); %comb ROI
[pk,~,~,~,~,comb] = test_fidupks(shotnumber,ncomb,fidu1,1,'a',1,xmin);
        image1a = timeaxiscorr1(im,pk,1,ncomb);  % streak after time axis linearization 
        fidu1a=image1a(ymin:ymax,xmin:xmax); %comb ROI after time axis linearization
[pk1a,xt1a,deltata,fidu1a,~,comba] = test_fidupks(shotnumber,ncomb,fidu1a,1,'a',1,xmin);
        xt1a=(1:size(im,2)-1)*deltata; %new time axis for comb. converts pixels/column# to times. -1 so first px= 0.
%(pk1a(9)-pk1a(3))*deltata
figure, imagesc(image1a,[0 thres]); %image1a
hb = imrect(gca, [200 850 140 50]);
posb = wait(hb);
xmin=round(posb(1,1)); xmax=round(xmin+posb(1,3)); ymin=round(posb(1,2)); ymax=round(ymin+posb(1,4));
fidu1b=image1a(ymin:ymax,xmin:xmax); %images1a comb ROI   
[pkb,~,deltatb,~,~,combb] = test_fidupks(shotnumber,ncomb,fidu1b,deltata,'b',1,xmin);
%(pkb(7)-pkb(3))*deltata
%% T0 calculation
if str2double(shotnumber) > 76000
    dt_sop=0.249; %ns after 2/12/15
    %dt_sop=0.16; %ns OLD after 2/12/15
else
    dt_sop=0.07; %ns before 2/12/15
end
if ep==1
%     g2=figure; plot(sum(image1a(400:450,:)));
%     t0px=input('Choose t0 pixel: ');
%     close(g2);
%     dt_sop=-t0px*deltata+pkb(2)*deltata+dt_offset-0.548;
dt_sop=-0.015;
end
Tzero = pkb(2)*deltata - dt_sop -dt_offset-0.548;  %positive to the right
xt1b=xt1a-Tzero;  %new time axis
figure, imagesc(xt1b,size(image1a,1),image1a,[0 thres])
 %% HDF output
    xt11=xt1b(1); imgg=flipdim(image1a,1)'; %refimg=flipdim(refdata(:,leftL:rightL),1)';
    [save_status , filename ]= save_SOP_HDFm( imgg,shotnumber, xt11, deltata,stordir,refshotnum); 
%%
% figure,imagesc(image1a,[0 thres])
% h = imrect(gca, [18 400 900 400]);
% pos = wait(h);
% xminr=round(pos(1,1)); xmaxr=round(xminr+pos(1,3)); yminr=round(pos(1,2)); ymaxr=round(yminr+pos(1,4));
% rprof=mean(image1a(yminr:ymaxr,xminr:xmaxr),1)';
% stdprof=std(image1a(yminr:ymaxr,xminr:xmaxr))';
% tsop=((xminr:xmaxr)-1)*deltata-Tzero;
% figure,
% %subplot(1, 2, 1);plot(tsop,rprof,'k',tsop,rprof+stdprof,':r',tsop,rprof-stdprof,':r')
% %subplot(1, 2, 1);
% figure,plot(tsop,rprof)
% xlim([-1 tsop(end)])
% ylim([0 1.2*max(rprof)])
% %set(gca,'PlotBoxAspectRatio',[1 1 1])
% xlabel('Time (ns)')
% ylabel('Intensity (ADU)')
% % %%
% % [T,dT]=SOPtemps(rprof,stdprof,1/deltata,800,0.7,0.3); %Int, dInt, eta, Ws, emissivity, ND
% % %subplot(1,2,2);plot(tsop,T,'k',tsop,T+dT,':r',tsop,T-dT,':r')
% % %subplot(1,2,2);
% % figure,plot(tsop,T)
% % xlim([-1 tsop(end)])
% % ylim([0 1.2*max(T)])
% % set(gca,'PlotBoxAspectRatio',[1 1 1])
% % xlabel('Time (ns)')
% % ylabel('Temperature (eV)')
