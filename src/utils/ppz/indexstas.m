% The purpose of this function is to index data by station and cast
% Created by E.Black, WHOI
% Last edited 12/7/2015

% Data in the file used must be sorted by cast and station (can't be out of
% order) or else the function won't work.  This code could be used for any
% set of data but it's main purpose for GEOTRACES is to have a master
% matrix of CTD data that is sorted by station/cast that can be used later
% by any code/file.  Once this matrix is created it shouldn't be changed.
% Therefore, the indices will always work because this file is not changed.

function [st,ca,Is,Ie,id] = indexstas(colst,colca,varargin)
%colst is the station column in your dataset
%colca is the cast column in your dataset
%the output is the station, cast, index start and index end for whichever
%file is of interest, e.g. GEOTRACES ODF downcast CTD data file
st=[]; ca=[]; Is=[]; Ie=[]; id=[];
stations = unique(colst);  %list of stations

switch nargin
    case 2
        for i = 1:length(stations);
            StIndex = find(colst==stations(i)); %find data grouped by station
            casts = unique(colca(StIndex(1):StIndex(end))); %list of casts each station
            for j = 1:length(casts);
                [m,n]=size(st);
                st(m+1,1)=stations(i); %Station number
                ca(m+1,1)=casts(j); %Cast number
                IndTemp = find(colst==stations(i) & colca==casts(j));
                Is(m+1,1) = IndTemp(1,1); %Start Index
                Ie(m+1,1) = IndTemp(end,1); %End Index
            end
        end
            clear stations casts StIndex i j m n IndTemp colst colca; 
            
    case 3  
        for i = 1:length(stations);
            StIndex = find(colst==stations(i)); %find data grouped by station
            casts=unique(colca(StIndex(1):StIndex(end))); %list of casts each station
            tmpid = varargin{1};
            tmpid = tmpid(StIndex(1):StIndex(end));
            for j = 1:length(casts);
                [m,n]=size(st);
                st(m+1,1)=stations(i); %Station number
                ca(m+1,1)=casts(j); %Cast number
                idc = find(colca(StIndex(1):StIndex(end))==casts(j),1,'first');
                IndTemp = find(colst==stations(i) & colca==casts(j));
                Is(m+1,1) = IndTemp(1,1); %Start Index
                Ie(m+1,1) = IndTemp(end,1); %End Index
                id(m+1,1) = tmpid(idc); %cast type identifier
            end
        end
            clear stations casts StIndex i j m n IndTemp colst colca;
            
    otherwise
        disp('Error, incorrect or too many variables');
end        
end
