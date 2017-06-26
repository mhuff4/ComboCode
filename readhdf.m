%File: readhdf.m
% Author : Bill donaldson
% Purpose:
%   Reads data from  an HDF file into data structure data
%   Symmetric WRT writehdf_v1.m
% Synopsis : [data , status]= readhdf_v1( filename)
%           filename=string.hdf
%           data=struct( 'file_attr', struct( 'attr1', attr1 ...) ,
%           'data1', struct( 'data1', data1, 'attr1', attr1, 'attr2',attr2 ...)
%           'data2', struct( 'data2', data2, 'attr1', attr1, 'attr2',attr2) ...)
%           status = -1 if fails
% Revisions:
%       10-8-09 to put file_attr in a separate sub structure
%       6-30-2011 add addition comments to file
%       9-5-2012 replace old readhdf.m with readhdf_v1.m
%-
%*************************************  writehdf ********
    

function [hdat , status ]=readhdf(filename)


%   File information gathering
S = hdfinfo(filename);
if isfield(S,'Attributes') > 0
	nattribs=length(S.Attributes); %equivalent to former lenattributes
    attrs= S.Attributes;
	for i=1:nattribs
        % read file attributes
        attr_name=attrs(1,i).Name;

        attr_name=regexprep(attr_name,' ','_'); % remove spaces from attribute name
        attr_value=attrs(1,i).Value;
        hdat.('file_attr').( attr_name)=attr_value;
	end
end

%   information gathering for a given specified set of data


nnsets =length(S.SDS);
sd_id = hdfsd('start', filename ,'read'); % open the file fro reading
for k=1:nnsets %upper bound = n_images if all are valid images
    lenindattributes = length(S.SDS(k).Attributes); %Added to match former output
    aa=S.SDS(k).Name;
    %disp( aa)
    sds_name=regexprep(aa,' ','_'); % remove space from  dataset names
    sds_id = hdfsd('select',sd_id, (k-1) );
    for ii = 1:lenindattributes
        attrname=S.SDS(k).Attributes(ii).Name;
        attrname=regexprep(attrname,' ','_'); % remove space from attribute  names
        hdat.((sds_name)).((attrname))=S.SDS(k).Attributes(ii).Value;
    end
    
    [~, ds_ndims, ds_dims, ~, ~, stat] =hdfsd('getinfo',sds_id);

    ds_start = zeros(1,ds_ndims); % Creates the vector [0 0]
    ds_stride = []; 
    ds_edges = ds_dims; 

    [ds_data, status] =hdfsd('readdata',sds_id,ds_start,ds_stride,ds_edges);  % read daata set from file
    if status >= 0
        hdat.((sds_name)).((sds_name))=ds_data; % save dataset in structure
        stat = hdfsd('endaccess',sds_id); % CLOSE ACCES TO DATASET
    else
        disp( [ 'Error in reading ' sds_name ' Status: ' num2str(status) num2str(stat) ' ' num2str(sds_id)] )
        
    end
            
            
end
%disp( [ 'closing ' filename ])
status=hdfsd('end',sd_id); % closes file	


