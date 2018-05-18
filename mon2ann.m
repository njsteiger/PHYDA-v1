function [Xann] = mon2ann(Xmon,smon,emon)
%MON2ANN  Compute the annual average given monthly inputs.
%   [Xann] = mon2ann(Xmon,smon,emon)
%        'Xmon' is the monthly input data, assumed to start in Jan and end
%        in Dec, 'smon' is the starting month (e.g. Jan=1), 'emon' is the
%        ending month (e.g. Dec=12). 'Xmon' can have dimensions
%        (lat,lon,time) or (space,time), or (time) only.
%
%        If the ending month is before the starting month (the average is
%        taken beginning in one year and carrying over into the next, e.g.,
%        Apr--Mar next year), one year's worth of data will be dropped to
%        accomodate this choice. If 'smon' and 'emon' are empty, no
%        averaging is performed.
%
%   Nathan Steiger, October 2015

if isempty(smon) || isempty(emon)
   % No averaging
   Xann=Xmon;
   
elseif isvector(Xmon)
   % Determine number of years of data
   numMonths=length(Xmon);
   yrs=numMonths/12;
   
   % Correct number of yrs if end < start (e.g. avg over Jul-->Jun)
   if smon > emon
      emon=emon+12;
      yrs=yrs-1;
   end
   
   Xann=zeros(yrs,1);
   for q=1:yrs
      Xann(q)=nansum(Xmon(smon:emon))/length(smon:emon);
      smon=smon+12;
      emon=emon+12;
   end
   
elseif length(size(Xmon))==2
   
   % Determine number of years of data
   [rows,numMonths]=size(Xmon);
   yrs=numMonths/12;
   
   % Correct number of yrs if end < start (e.g. avg over Jul-->Jun)
   if smon > emon
      emon=emon+12;
      yrs=yrs-1;
   end
   
   Xann=zeros(rows,yrs);
   for q=1:yrs
      Xann(:,q)=nansum(Xmon(:,smon:emon),2)/length(smon:emon);
      smon=smon+12;
      emon=emon+12;
   end
   
elseif length(size(Xmon))==3
   
   % Determine number of years of data
   [rows,cols,numMonths]=size(Xmon);
   yrs=numMonths/12;
   
   % Correct number of yrs if end < start (e.g. avg over Jul-->Jun)
   if smon > emon
      emon=emon+12;
      yrs=yrs-1;
   end
   
   Xann=zeros(rows,cols,yrs);
   for q=1:yrs
      Xann(:,:,q)=nansum(Xmon(:,:,smon:emon),3)/length(smon:emon);
      smon=smon+12;
      emon=emon+12;
   end
   
else
   error('Incorrect input for Xmon...')
end


end

