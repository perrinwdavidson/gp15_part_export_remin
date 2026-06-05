% The purpose of this function is to determine the 'particle production
% zone' or 'ppz' depth.  This is a proxy for euphotic depth and an
% explanation for its use is in Owens et al. (2015, DSRII).  We use this
% when there is no reliable PAR data, as on GEOTRACES Atlantic and Pacific

% Based on talks at the Catalina EPZT meeting the 2ndary Chl max (deep max)
% will not be considered in the calcuation.  The first depth at which 10%
% of the fluorescence maximum is reach becomes the PPZ depth.

% Created by E.Black, WHOI
% Last edited 1/12/2016

% Data in the file used must be sorted by cast and station (can't be out of
% order) or else the function won't work.  Could add a sort to code later
% if desired.

%This function should be used with indexstas.m to indicate stations and
%casts to use

function [maxf,maxd,ppz,maxd2,bkg] = findppz(st,Ie,Is,cold,colf)
%station is from indexstas.m, its the list of unique stations + casts
%cold is the column of your data that has depth (m)
%colf is the column of your data that has fluorescence (volts)
maxf=[]; maxd=[]; ppz=[]; maxd2=[]; bkg=[];
for i = 1:size(st);
    % Make depth and fluorescence vectors by station
    fdata = colf(Is(i):Ie(i),1);
    depth = cold(Is(i):Ie(i),1);
    if ~isnan(fdata) & max(depth)>100 & min(depth)<150;
        % Find general background value and subtract from data, smooth data
        pct = floor(.1*length(fdata)); 
        bkgtmp = mean(fdata(end-pct:end));
        [m,n]=size(bkg); bkg(m+1,1)=bkgtmp; %add avg background value as check
        fdataADJpre =(fdata-bkgtmp); %subtract out background, make 'adjusted' data
        fdataADJ=medfilt1(fdataADJpre); %smooth data with a median filter
        %Biogeosciences, 9, 2111–2125, 2012, Lavigne
        %Alternate option is 3 span moving average, use 'smooth(fdataADJpre,3)'
    
        % Find ppz using only the top of the profile (in case something
        % happened to the sensor at depth, ppz in oligotrophic N.Atlantic was 
        % 250 meters maximum so 300 meters should be extreme case for anywhere
        exdpth=300;
        if length(fdataADJ)<exdpth %Only use top 300 m
         [maxF,maxindex]=max(fdataADJ);
        else [maxF,maxindex]=max(fdataADJ(1:exdpth));
        end
        maxf(m+1,1)=maxF; %add max value of fluorescene (volts)
        maxd(m+1,1)=depth(maxindex); %add max depth of fluorescene (m)
        fdataPctMax = (fdataADJ./maxF)*100; %determine all in % of max
        [PPZind]=find(fdataPctMax(maxindex+1:end)<=10, 1, 'first'); %10 of max
	ppz_calc = depth(PPZind+maxindex); %PPZ depth (m)
	if isempty(ppz_calc)
		ppz_calc = -9999;
	end
        ppz(m+1,1)= ppz_calc;  % depth(PPZind+maxindex); %PPZ depth (m)  
	% disp(['i = ', num2str(i), ' | PPZ = ', num2str(depth(PPZind+maxindex))])
% Are there any additional 2ndary peaks?
        [maxF2,maxindex2]=max(fdataPctMax((PPZind+maxindex+1):end));
        if  maxF2<20; %20% limit so points very close to the peak don't count
        maxd2 = NaN; % 2ndardy max doesn't exist
        else
         %uncomment and make a variable to get 2ndary max value in volts
         %xxxx = fdataADJ(maxindex2+PPZind+maxindex+1);
            maxd2 = depth(maxindex2+PPZind+maxindex); % 2ndardy max depth (m)
        end
    
        %To check PPZ data uncomment to make figures for each station/cast with
        %smoothed data and raw data.
    %figure(i)
    % plot(fdataADJ,depth,'-r','linewidth',3);
    % hold on;
    % title(['Station ' num2str(st(i))]); 
    % plot(fdataADJpre,depth,'-k');
    fmaxd = maxd;
    ppzd = ppz;
    %plot([0 max(fdataADJpre)],[fmaxd fmaxd],'--b');
    %plot([0 max(fdataADJpre)],[ppzd ppzd],'--g');
    % ylim([0 500]);
    % x_max=get(gca,'xlim');
    % xlim([0 x_max(2)]);
    % set(gca,'Ydir','reverse')
    % hold off;
    
    else
        [m,n]=size(bkg); bkg(m+1,1)=NaN;
        maxf(m+1,1)=NaN; maxd(m+1,1)=NaN; ppz(m+1,1)=NaN; maxd2(m+1,1)=NaN;
    end
end

