%File: writehdf.m
% Author : Bill donaldson
% Purpose:
%   Saves to an HDF file all information in data structure data
% Synopsis : [ save_status  ]= writehdf( filename, data)
%           filename=string.hdf
%           data=struct( 'file_attr', struct( 'attr1', attr1 ...) ,
%           'data1', struct( 'data1', data1, 'attr1', attr1, 'attr2',attr2 ...)
%           'data2', struct( 'data2', data2, 'attr1', attr1, 'attr2',attr2) ...)
% Revisions:
%-
%*************************************  writehdf ********
function [ save_status  ]= writehdf( filename, data);
    global sds_id  sd_id
    % Create new file
    currently=now();
    str_date=datestr( currently , 30 );
    %fileDir='c:\Tek_IR_data\'
    fnames=fieldnames( data ); %cell array
    nfields=length( fnames );
    sd_id=hdfsd('start', filename ,'dfacc_create'); % open file for writing

    %sets standard file attributes for the hdf file
    astat=hdfsd('setattr',sd_id,'Time_Stamp', str_date);
    astat=hdfsd('setattr',sd_id,'Program_verison', 'writehdf.m 3-14-2008');

    for i=1:nfields
        fname=cell2mat( fnames(i));
        if strcmpi( fname , 'FILE_ATTR')
            % this ia a file attribute
            fattr=data.( fname );
            if isstruct(fattr )
                anames=fieldnames( fattr);
                nattr=length( anames );
                for j=1: nattr
                    aname=cell2mat(anames(j));
                    astat=hdfsd('setattr',sd_id,aname, fattr.(aname));
                end
            else
                astat=hdfsd('setattr',sd_id,fname, fattr);
            end
        else
            % this is a dataset
            dataset=data.(fname);

            if isstruct( dataset )
                ds_ids=fieldnames( dataset);
                nids=length( ds_ids );
                ds_indx= strcmpi( ds_ids , fname ) ;
                if sum( ds_indx ) == 1
                    % there should be only 1 field of the dataset with the
                    % same name as the dataset
                    ds=data.(fname).(fname);
                    %creates data set in hdf file
                    ds_type=class( ds );
                    if (ds_type == 'single')
                        ds_type='DFNT_FLOAT32';
                    else
                       ds_type='DFNT_FLOAT64';
                    end
                    ds_rank = ndims( ds );
                    ds_dims = fliplr(size( ds ));
                    sds_id = hdfsd('create',sd_id, fname ,ds_type,ds_rank,ds_dims);
                    ds_start = zeros(1,ndims( ds));
                    ds_stride = [];
                    ds_edges = fliplr(size( ds ));
                    %writes DATA to the data set
                    dstat = hdfsd('writedata', sds_id,...
                            ds_start, ds_stride, ds_edges, ds );
                else
                    disp( [ 'Dataset :' fname ' should have a sub-field called ' fname ])
                    whos( 'dataset' )
                end
                for j=1: nids
                   ds_id=cell2mat( ds_ids(j));
                   if  not( strcmp( ds_id , fname ) )
                       % this is not the dataset so save it as such
                       % this is an attribute so treat it as such
                       astat=hdfsd('setattr',sds_id, ds_id, data.(fname).(ds_id) );
                       %%%%%
                   end
                end
            else
                disp( 'datasets must be substructures of the input structure data so treat this as a file attribute')
                fattr=data.( fname );
                astat=hdfsd('setattr',sd_id,fname, fattr);
            end
            estat=hdfsd('endaccess',sds_id); % closes dataset
        end
    end
    statB=hdfsd('end',sd_id); % closes file
    % File saved display message
    if (dstat >= 0 ) & ( astat >= 0) & ( estat >=0 )
        save_status=1;
        disp( ['Data saved to: ' filename])
    else
        dstat
        astat
        estat
        statB
        save_status=0;
        disp( ['failed saved to: ' filename])
    end
