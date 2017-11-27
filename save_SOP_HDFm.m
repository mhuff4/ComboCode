%File: save_ASBO_HDF.m
% Author : Bill donaldson
% Purpose:
%   Saves to an HDF file all information related to the ASBO file, adapted
%   from IR regen HDF file saver
% Synopsis : save_regen( process_ir , dt, ch , E_ir_in , driver , SSD, ...
%           power_setting,  pulseshape , fastframe , predict_uv , raw_ir , ref_ch, ...
%           tstamp, Preamble, xzero ,   shot)
%   Inputs: ASBOimg, refimg, gridimg , shotnumber, Asbo , xt1a
%           
% outputs:
%               save_status = 1 if successful
%               filename - namwe of file written to
% Revisions:
%            3/18/2008  now used writehdf
%            6-24-2008 save the preamble
%            01-08-2009 save xzero as T0
%            02-06-2009 save current file path as file attr
%            02-06-2009 Wolf Seka
%-
%*************************************  save_regen ********

function [save_status , filename ]= save_SOP_HDFm( SOPimg,shotnumber,xt11,deltat,stordir,refshotnum)

    % Create new file
    currently=now();
    str_date=datestr( currently , 30 );
    
    % create data structure to write to file
    %sets attributes for the hdf file
    
    pver=mfilename('fullpath'); % get the path of the current m-file
    dirver=dir( [pver '.m']);  % find when file was modifieid
    dv=dirver.date;
    p_ver=[  pver ' modified ' dv ]; 
    %Put back 'Processing_Message',T0mssg, after ...
    file_attr=struct( 'SOP_reduction', 'SOP' , 'Time_Stamp', datestr(now), ...
             'Program_version', p_ver);
    
    Streak_array=struct( 'Streak_array' , SOPimg , 'shotnumber', shotnumber , ...
        't0_of_first_px_in_ns', xt11, 'deltat_per_px_in_ns', deltat,'Dist_corr_made_with_shot', refshotnum);
  
    Data=struct( 'FILE_ATTR', file_attr , 'Streak_array' , Streak_array);
 
    
%stordir='\\HOPI-fs\shares\Experimental\SOP\HDF-files\';
%filename=[stordir 'ASBO' num2str(Asbo) '_' num2str(shotnumber) '.hdf'];
filename=[stordir 'sop_sop-ross_s',shotnumber,'_proc.hdf'];
save_status =writehdf( filename , Data ) ;         
    % File saved display message
    if save_status == 1
       
        disp( ['Data saved to: ' filename])
    else
        
        disp( ['failed saved to: ' filename])
    end

