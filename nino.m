function [nino_indx,nmeta] = nino(SST,lat,lon,indx,bc)
%NINO Compute Nino indices given model temperature data
%   [nino_indx] = nino(SST,lat,indx)
%       Dimensions of SST(lat,lon)
%       indx indicates which nino index to compute: 1 (Nino 1+2), 2 (Nino
%       3), 3 (Nino 3.4), 4 (Nino 4), 5 cross-equatorial gradient
%
%   Nathan Steiger, June 2016; July 2017; Nov 2017


if indx==1 % Nino 1+2 = Area average over 10S-0N and 90W-80W
    lt1=-10;lt2=0;ln1=360-90;ln2=360-80;   
    nmeta.lat1=lt1;nmeta.lat2=lt2;nmeta.lon1=ln1;nmeta.lon2=ln2;
    nmeta.indxnm='Nino 1+2';
 
    % Find nearest lat/lon equal to or outside bounds
    lt_a=find((abs(lat-lt1)-min(abs(lat-lt1)))<1e-5);
    lt_b=find((abs(lat-lt2)-min(abs(lat-lt2)))<1e-5);
    % if two equal mins, pick one outside of bounds
    lt_a=min(lt_a);lt_b=max(lt_b);
    
    ln_a=find((abs(lon-ln1)-min(abs(lon-ln1)))<1e-5);
    ln_b=find((abs(lon-ln2)-min(abs(lon-ln2)))<1e-5);
    % if two equal mins, pick one outside of bounds
    ln_a=min(ln_a);ln_b=max(ln_b);

 
    % Compute spatial mean
    nino_sm=wmean(SST(lt_a:lt_b,ln_a:ln_b,:),lat(lt_a:lt_b));
    
    clear lt_a lt_b ln_a ln_b
    
    
    if bc=='y'
        
      %--------------------------------------------------------------------------------------------------
      % Bias-correct using BEarth climatology (1951-1980), which over ocean is just hi-resolution HadISST
      %--------------------------------------------------------------------------------------------------
      tclim=ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','climatology');
      tclim=permute(tclim,[2 1 3]);
      xlon1 = double(ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','longitude'));
      xlat1 = double(ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','latitude'));
      [tclim,xlon]=regrid2gcm(tclim,xlon1);
      xlat=xlat1(:);xlon=xlon(:);

      % Find nearest lat/lon equal to or outside bounds
      lt_a=find((abs(xlat-lt1)-min(abs(xlat-lt1)))<1e-5);
      lt_b=find((abs(xlat-lt2)-min(abs(xlat-lt2)))<1e-5);
      % if two equal mins, pick one outside of bounds
      lt_a=min(lt_a);lt_b=max(lt_b);

      ln_a=find((abs(xlon-ln1)-min(abs(xlon-ln1)))<1e-5);
      ln_b=find((abs(xlon-ln2)-min(abs(xlon-ln2)))<1e-5);
      % if two equal mins, pick one outside of bounds
      ln_a=min(ln_a);ln_b=max(ln_b);


      % Compute spatial mean of climatology
      t_nino=wmean(tclim(lt_a:lt_b,ln_a:ln_b,:),xlat(lt_a:lt_b));

      nino_mon=reshape(nino_sm,12,length(nino_sm)/12);
      nino_dseas=bsxfun(@minus,nino_mon,mean(nino_mon,2)); % de-season
      nino_bc=bsxfun(@plus,nino_dseas,t_nino); % add in obs clim

      % Reshape nino variable
      nino_indx=reshape(nino_bc,length(nino_sm),1);
      nmeta.nino_indx=nino_indx;        

    elseif bc=='n'
      % No bias-correction
      nino_indx=nino_sm;
      nmeta.nino_indx=nino_indx;        
    else
      error('Wrong Nino bias-correction option...')
    end
    
elseif indx==2 % Nino 3 = Area average over 5S-5N and 150W-90W
    lt1=-5;lt2=5;ln1=360-150;ln2=360-90;   
    nmeta.lat1=lt1;nmeta.lat2=lt2;nmeta.lon1=ln1;nmeta.lon2=ln2;
    nmeta.indxnm='Nino 3';
 
    % Find nearest lat/lon equal to or outside bounds
    lt_a=find((abs(lat-lt1)-min(abs(lat-lt1)))<1e-5);
    lt_b=find((abs(lat-lt2)-min(abs(lat-lt2)))<1e-5);
    % if two equal mins, pick one outside of bounds
    lt_a=min(lt_a);lt_b=max(lt_b);
    
    ln_a=find((abs(lon-ln1)-min(abs(lon-ln1)))<1e-5);
    ln_b=find((abs(lon-ln2)-min(abs(lon-ln2)))<1e-5);
    % if two equal mins, pick one outside of bounds
    ln_a=min(ln_a);ln_b=max(ln_b);

 
    % Compute spatial mean
    nino_sm=wmean(SST(lt_a:lt_b,ln_a:ln_b,:),lat(lt_a:lt_b));
    
    clear lt_a lt_b ln_a ln_b
    
    
    if bc=='y'
        
      %--------------------------------------------------------------------------------------------------
      % Bias-correct using BEarth climatology (1951-1980), which over ocean is just hi-resolution HadISST
      %--------------------------------------------------------------------------------------------------
      tclim=ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','climatology');
      tclim=permute(tclim,[2 1 3]);
      xlon1 = double(ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','longitude'));
      xlat1 = double(ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','latitude'));
      [tclim,xlon]=regrid2gcm(tclim,xlon1);
      xlat=xlat1(:);xlon=xlon(:);

      % Find nearest lat/lon equal to or outside bounds
      lt_a=find((abs(xlat-lt1)-min(abs(xlat-lt1)))<1e-5);
      lt_b=find((abs(xlat-lt2)-min(abs(xlat-lt2)))<1e-5);
      % if two equal mins, pick one outside of bounds
      lt_a=min(lt_a);lt_b=max(lt_b);

      ln_a=find((abs(xlon-ln1)-min(abs(xlon-ln1)))<1e-5);
      ln_b=find((abs(xlon-ln2)-min(abs(xlon-ln2)))<1e-5);
      % if two equal mins, pick one outside of bounds
      ln_a=min(ln_a);ln_b=max(ln_b);


      % Compute spatial mean of climatology
      t_nino=wmean(tclim(lt_a:lt_b,ln_a:ln_b,:),xlat(lt_a:lt_b));

      nino_mon=reshape(nino_sm,12,length(nino_sm)/12);
      nino_dseas=bsxfun(@minus,nino_mon,mean(nino_mon,2)); % de-season
      nino_bc=bsxfun(@plus,nino_dseas,t_nino); % add in obs clim

      % Reshape nino variable
      nino_indx=reshape(nino_bc,length(nino_sm),1);
      nmeta.nino_indx=nino_indx;        
        
    elseif bc=='n'
      % No bias-correction
      nino_indx=nino_sm;
      nmeta.nino_indx=nino_indx;        
    else
      error('Wrong Nino bias-correction option...')
    end
    
elseif indx==3 % Nino 3.4 = Area average over 5S-5N and 170W-120W
    lt1=-5;lt2=5;ln1=360-170;ln2=360-120;   
    nmeta.lat1=lt1;nmeta.lat2=lt2;nmeta.lon1=ln1;nmeta.lon2=ln2;
    nmeta.indxnm='Nino 3.4';
 
    % Find nearest lat/lon equal to or outside bounds
    lt_a=find((abs(lat-lt1)-min(abs(lat-lt1)))<1e-5);
    lt_b=find((abs(lat-lt2)-min(abs(lat-lt2)))<1e-5);
    % if two equal mins, pick one outside of bounds
    lt_a=min(lt_a);lt_b=max(lt_b);
    
    ln_a=find((abs(lon-ln1)-min(abs(lon-ln1)))<1e-5);
    ln_b=find((abs(lon-ln2)-min(abs(lon-ln2)))<1e-5);
    % if two equal mins, pick one outside of bounds
    ln_a=min(ln_a);ln_b=max(ln_b);

 
    % Compute spatial mean
    nino_sm=wmean(SST(lt_a:lt_b,ln_a:ln_b,:),lat(lt_a:lt_b));
    
    clear lt_a lt_b ln_a ln_b
    
    
    if bc=='y'
        
      %--------------------------------------------------------------------------------------------------
      % Bias-correct using BEarth climatology (1951-1980), which over ocean is just hi-resolution HadISST
      %--------------------------------------------------------------------------------------------------
      tclim=ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','climatology');
      tclim=permute(tclim,[2 1 3]);
      xlon1 = double(ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','longitude'));
      xlat1 = double(ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','latitude'));
      [tclim,xlon]=regrid2gcm(tclim,xlon1);
      xlat=xlat1(:);xlon=xlon(:);

      % Find nearest lat/lon equal to or outside bounds
      lt_a=find((abs(xlat-lt1)-min(abs(xlat-lt1)))<1e-5);
      lt_b=find((abs(xlat-lt2)-min(abs(xlat-lt2)))<1e-5);
      % if two equal mins, pick one outside of bounds
      lt_a=min(lt_a);lt_b=max(lt_b);

      ln_a=find((abs(xlon-ln1)-min(abs(xlon-ln1)))<1e-5);
      ln_b=find((abs(xlon-ln2)-min(abs(xlon-ln2)))<1e-5);
      % if two equal mins, pick one outside of bounds
      ln_a=min(ln_a);ln_b=max(ln_b);


      % Compute spatial mean of climatology
      t_nino=wmean(tclim(lt_a:lt_b,ln_a:ln_b,:),xlat(lt_a:lt_b));

      nino_mon=reshape(nino_sm,12,length(nino_sm)/12);
      nino_dseas=bsxfun(@minus,nino_mon,mean(nino_mon,2)); % de-season
      nino_bc=bsxfun(@plus,nino_dseas,t_nino); % add in obs clim

      % Reshape nino variable
      nino_indx=reshape(nino_bc,length(nino_sm),1);
      nmeta.nino_indx=nino_indx;        
        
    elseif bc=='n'
      % No bias-correction
      nino_indx=nino_sm;
      nmeta.nino_indx=nino_indx;        
    else
      error('Wrong Nino bias-correction option...')
    end
    
elseif indx==4 % Nino 4 = Area average over 5S-5N and 160E-150W
    lt1=-5;lt2=5;ln1=160;ln2=360-150;   
    nmeta.lat1=lt1;nmeta.lat2=lt2;nmeta.lon1=ln1;nmeta.lon2=ln2;
    nmeta.indxnm='Nino 4';
 
    % Find nearest lat/lon equal to or outside bounds
    lt_a=find((abs(lat-lt1)-min(abs(lat-lt1)))<1e-5);
    lt_b=find((abs(lat-lt2)-min(abs(lat-lt2)))<1e-5);
    % if two equal mins, pick one outside of bounds
    lt_a=min(lt_a);lt_b=max(lt_b);
    
    ln_a=find((abs(lon-ln1)-min(abs(lon-ln1)))<1e-5);
    ln_b=find((abs(lon-ln2)-min(abs(lon-ln2)))<1e-5);
    % if two equal mins, pick one outside of bounds
    ln_a=min(ln_a);ln_b=max(ln_b);

 
    % Compute spatial mean
    nino_sm=wmean(SST(lt_a:lt_b,ln_a:ln_b,:),lat(lt_a:lt_b));
    
    clear lt_a lt_b ln_a ln_b
    
    
    if bc=='y'
        
      %--------------------------------------------------------------------------------------------------
      % Bias-correct using BEarth climatology (1951-1980), which over ocean is just hi-resolution HadISST
      %--------------------------------------------------------------------------------------------------
      tclim=ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','climatology');
      tclim=permute(tclim,[2 1 3]);
      xlon1 = double(ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','longitude'));
      xlat1 = double(ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','latitude'));
      [tclim,xlon]=regrid2gcm(tclim,xlon1);
      xlat=xlat1(:);xlon=xlon(:);

      % Find nearest lat/lon equal to or outside bounds
      lt_a=find((abs(xlat-lt1)-min(abs(xlat-lt1)))<1e-5);
      lt_b=find((abs(xlat-lt2)-min(abs(xlat-lt2)))<1e-5);
      % if two equal mins, pick one outside of bounds
      lt_a=min(lt_a);lt_b=max(lt_b);

      ln_a=find((abs(xlon-ln1)-min(abs(xlon-ln1)))<1e-5);
      ln_b=find((abs(xlon-ln2)-min(abs(xlon-ln2)))<1e-5);
      % if two equal mins, pick one outside of bounds
      ln_a=min(ln_a);ln_b=max(ln_b);


      % Compute spatial mean of climatology
      t_nino=wmean(tclim(lt_a:lt_b,ln_a:ln_b,:),xlat(lt_a:lt_b));

      nino_mon=reshape(nino_sm,12,length(nino_sm)/12);
      nino_dseas=bsxfun(@minus,nino_mon,mean(nino_mon,2)); % de-season
      nino_bc=bsxfun(@plus,nino_dseas,t_nino); % add in obs clim

      % Reshape nino variable
      nino_indx=reshape(nino_bc,length(nino_sm),1);
      nmeta.nino_indx=nino_indx;        
        
    elseif bc=='n'
      % No bias-correction
      nino_indx=nino_sm;
      nmeta.nino_indx=nino_indx;        
    else
      error('Wrong Nino bias-correction option...')
    end
    
elseif indx==5 % Cross equatorial SST gradient = diff of boxes at  5S-5N and 150E->160W minus 130W->80W
    lt1=-5;lt2=5;ln1=150;ln2=360-160;ln3=360-130;ln4=360-80;   
    nmeta.lat1=lt1;nmeta.lat2=lt2;nmeta.lon1=ln1;nmeta.lon2=ln2;nmeta.lon3=ln3;nmeta.lon4=ln4;
    nmeta.indxnm='Pacific del SST: west minus east';
 
    % Find nearest lat/lon equal to or outside bounds
    lt_a=find((abs(lat-lt1)-min(abs(lat-lt1)))<1e-5);
    lt_b=find((abs(lat-lt2)-min(abs(lat-lt2)))<1e-5);
    % if two equal mins, pick one outside of bounds
    lt_a=min(lt_a);lt_b=max(lt_b);
    
    ln_a=find((abs(lon-ln1)-min(abs(lon-ln1)))<1e-5);
    ln_b=find((abs(lon-ln2)-min(abs(lon-ln2)))<1e-5);
    % if two equal mins, pick one outside of bounds
    ln_a=min(ln_a);ln_b=max(ln_b);

    ln_c=find((abs(lon-ln3)-min(abs(lon-ln3)))<1e-5);
    ln_d=find((abs(lon-ln4)-min(abs(lon-ln4)))<1e-5);
    % if two equal mins, pick one outside of bounds
    ln_c=min(ln_c);ln_d=max(ln_d);

 
    % Compute spatial mean in each box, compute the SST gradient
    ebox=wmean(SST(lt_a:lt_b,ln_a:ln_b,:),lat(lt_a:lt_b));
    wbox=wmean(SST(lt_a:lt_b,ln_c:ln_d,:),lat(lt_a:lt_b));
    m_delsst=ebox-wbox;
 
    clear lt_a lt_b ln_a ln_b ln_c ln_d
    
    
    if bc=='y'
        
      %--------------------------------------------------------------------------------------------------
      % Bias-correct using BEarth climatology (1951-1980), which over ocean is just hi-resolution HadISST
      %--------------------------------------------------------------------------------------------------
      tclim=ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','climatology');
      tclim=permute(tclim,[2 1 3]);
      xlon1 = double(ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','longitude'));
      xlat1 = double(ncread('/d1/nsteiger/climate-data/berkeley-earth/Land_and_Ocean_LatLong1_1850_2016.nc','latitude'));
      [tclim,xlon]=regrid2gcm(tclim,xlon1);
      xlat=xlat1(:);xlon=xlon(:);

      % Find nearest lat/lon equal to or outside bounds
      lt_a=find((abs(xlat-lt1)-min(abs(xlat-lt1)))<1e-5);
      lt_b=find((abs(xlat-lt2)-min(abs(xlat-lt2)))<1e-5);
      % if two equal mins, pick one outside of bounds
      lt_a=min(lt_a);lt_b=max(lt_b);

      ln_a=find((abs(xlon-ln1)-min(abs(xlon-ln1)))<1e-5);
      ln_b=find((abs(xlon-ln2)-min(abs(xlon-ln2)))<1e-5);
      % if two equal mins, pick one outside of bounds
      ln_a=min(ln_a);ln_b=max(ln_b);

      ln_c=find((abs(xlon-ln3)-min(abs(xlon-ln3)))<1e-5);
      ln_d=find((abs(xlon-ln4)-min(abs(xlon-ln4)))<1e-5);
      % if two equal mins, pick one outside of bounds
      ln_c=min(ln_c);ln_d=max(ln_d);


      % Compute spatial mean in each box, compute the SST gradient
      ebox=wmean(tclim(lt_a:lt_b,ln_a:ln_b,:),xlat(lt_a:lt_b));
      wbox=wmean(tclim(lt_a:lt_b,ln_c:ln_d,:),xlat(lt_a:lt_b));
      t_delsst=ebox-wbox;


      sst_mon=reshape(m_delsst,12,length(m_delsst)/12);
      sst_dseas=bsxfun(@minus,sst_mon,mean(sst_mon,2)); % de-season
      sst_bc=bsxfun(@plus,sst_dseas,t_delsst); % add in obs clim

      % Reshape nino variable
      nino_indx=reshape(sst_bc,length(m_delsst),1);
      nmeta.nino_indx=nino_indx;        
        
    elseif bc=='n'
      % No bias-correction
      nino_indx=m_delsst;
      nmeta.nino_indx=nino_indx;        
    else
      error('Wrong Nino bias-correction option...')
    end
    
    
else
    error('Option not yet available...')
end


% % LOOK AT SPECTRA OF ENSO
%
% % detrend data
% ad=detrend(nino_indx);
%
% % get actual spectrum
% [pxx_a,f_a]=periodogram(ad,hamming(length(ad)),[],12);
%
% % NINO SPECTRA
% figure
% hold on
% plot(f_a,pxx_a,'b','linewidth',2)
% %plot(1./f_a,smooth(pntws(:,2),25),'--r','linewidth',2)
% %plot(1./f_a,smooth(mean(pxx,2),25),'r','linewidth',2)
% title('CESM','fontsize',18,'fontweight','bold')
% legend('Nino3.4')%,'RN 95%','RN mean')
% %xlim([0 ln])
% %xlim([0 0.5])
% xlabel('Frequency (1/yr)','fontsize',14)
% ylabel('Power','fontsize',14)
% set(gca,'fontsize',14)
% box on

% % COMPARE WITH REAL ENSO DATA
%
% nino34_hist=load('nino34_data.txt');
% nv=reshape(nino34_hist(:,2:13)',66*12,1);
%
% % detrend data
% ad=detrend(nv);
%
% % get actual spectrum
% % [pxx,f] = periodogram(x,window,f,fs) returns the two-sided periodogram
% % estimates at the frequencies specified in the vector, f. f must contain
% % at least two elements. The frequencies in f are in cycles per unit time.
% % The sampling frequency, fs, is the number of samples per unit time. If
% % the unit of time is seconds, then f is in cycles/second (Hz).
% [pxx_a,f_a]=periodogram(ad,hamming(length(ad)),[],12);
%
% % NINO SPECTRA
% ln=100;
% figure
% hold on
% plot(f_a,pxx_a,'b','linewidth',2)
% %plot(1./f_a,smooth(pntws(:,2),25),'--r','linewidth',2)
% %plot(1./f_a,smooth(mean(pxx,2),25),'r','linewidth',2)
% title('Real data','fontsize',18,'fontweight','bold')
% legend('Nino3.4')%,'RN 95%','RN mean')
% %xlim([0 0.5])
% %xlim([0 ln])
% xlabel('Frequency (1/yr)','fontsize',14)
% ylabel('Power','fontsize',14)
% set(gca,'fontsize',14)
% box on

end



