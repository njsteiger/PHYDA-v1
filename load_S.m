function [X,xmeta] = load_S(mdl_c,state_tp,mon_avg_o,mon_avg_f,m_yrs)
%LOAD_S  Load climate model or reanalysis data for state vector
%   [X,xmeta] = load_S(model,state_tp,mon_avg_o,mon_avg_f,m_yrs)
%     Given a model code 'mdl_c', state vector variable type 'state_tp',
%     start and end months for averaging, and the years from the model
%     'm_yrs', this function will load the data saved to the computer
%     'humid'. The outputs are the state vector 'X', the latitudes and
%     longitudes of the gridded data, and the locations of the state
%     variables 'id_X' which are in the same order as given in 'state_tp'.
%     The state variables are annually averaged according to the start
%     month and end moth inputs, unless these are empty in which case the
%     original monthly data is returned. Also available for certain
%     variables are the ocean lats/lons: olat/olon
%
%     The possible data source codes are:
%     'ccsm4_lm'
%     'ccsm4_2mpi'
%     'ccsm4_2gfdlpic'
%     'ccsm4_pic'
%     'gfdlcm3_pic'
%     'gfdlcm3_hist'
%     'gfdlesm2g_pi'
%     'gfdlesm2m_pi'
%     'mpi_lm'
%     'echam5'
%     'echam5_2cam5iso'
%     'echam5_2speedy'
%     'speedy_lm'
%     'icam5_hist'
%     All full forcing simulations from 'cesm_lmeXXX'
%
%   Nathan Steiger, October 2015; Nov 2017

X=[];%id_X=zeros(2*length(state_tp),1);
xmeta=cell(0);

if strcmp(mdl_c(1:8),'cesm_lme')
    
    ens=mdl_c(9:11);
    
    lat=ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TREFHT.085001-184912.nc'],'lat');
    lon=ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TREFHT.085001-184912.nc'],'lon');
    
    for j=1:length(state_tp)
        switch state_tp(j)
            case 't'
                % Surface temperature
                Xv_mon=ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TREFHT.085001-184912.nc'],'TREFHT');
                Xv_mon=permute(Xv_mon,[2 1 3]);

		disp('Bias correcting surface temperature to climatology')
		[X_clim] = load_obsclim('bearth',lat,lon);
		[nt,nn,nm]=size(Xv_mon);
		Xv_mon12=reshape(Xv_mon,nt,nn,12,nm/12);
		Xv_ns=bsxfun(@minus,Xv_mon12,mean(Xv_mon12,4));% remove climatology
		Xv_bc=bsxfun(@plus,Xv_ns,X_clim); % replace climatology
		Xv_mon=reshape(Xv_bc,nt,nn,nm);	

		Xv=mon2ann(Xv_mon,mon_avg_o,mon_avg_f);
                

		[rows,cols,len]=size(Xv(:,:,m_yrs));
                id_x0=size(X,1)+1;
                X=cat(1,X,reshape(Xv(:,:,m_yrs),rows*cols,len));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
		xmeta{j}.varnmlng='2 m temperature';
		xmeta{j}.lat=lat;
		xmeta{j}.lon=lon;
		xmeta{j}.id_X=[id_x0,id_x1];
            case '2'
                % Global mean 2m temperature
                Xv_mon=ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TREFHT.085001-184912.nc'],'TREFHT');
                Xv_mon=permute(Xv_mon,[2 1 3]);

		disp('GMT: Bias correcting surface temperature to climatology')
		[X_clim] = load_obsclim('bearth',lat,lon);
		[nt,nn,nm]=size(Xv_mon);
		Xv_mon12=reshape(Xv_mon,nt,nn,12,nm/12);
		Xv_ns=bsxfun(@minus,Xv_mon12,mean(Xv_mon12,4));% remove climatology
		Xv_bc=bsxfun(@plus,Xv_ns,X_clim); % replace climatology
		Xv_mon=reshape(Xv_bc,nt,nn,nm);	

		Xv=mon2ann(Xv_mon,mon_avg_o,mon_avg_f);
                
                
		gmt=wmean(Xv,lat);
                len=length(gmt(m_yrs));
                id_x0=size(X,1)+1;
                X=cat(1,X,reshape(gmt(m_yrs),1,len));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
		xmeta{j}.varnmlng='global mean 2 m temperature';
		xmeta{j}.id_X=[id_x0,id_x1];
            case 's'
                % SST
                olat=permute(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SST.085001-184912.nc'],'TLAT'),[2 1]);
                olon=permute(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SST.085001-184912.nc'],'TLONG'),[2 1]);
                Xv_mon=squeeze(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SST.085001-184912.nc'],'SST'));
                Xv=mon2ann(Xv_mon,mon_avg_o,mon_avg_f);
                Xv=permute(Xv,[2 1 3]);
                Xv(Xv==0)=NaN; % replace nans
                [rows,cols,len]=size(Xv(:,:,m_yrs));
                id_x0=size(X,1)+1;
                X=cat(1,X,reshape(Xv(:,:,m_yrs),rows*cols,len));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
		xmeta{j}.varnmlng='sea surface temperature';
		xmeta{j}.olat=olat;
		xmeta{j}.olon=olon;
		xmeta{j}.id_X=[id_x0,id_x1];
            case 'y'
 		% Sea surface salinit(y)
                olat=permute(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SSS.085001-184912.nc'],'TLAT'),[2 1]);
                olon=permute(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SSS.085001-184912.nc'],'TLONG'),[2 1]);
                Xv_mon=squeeze(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SSS.085001-184912.nc'],'SSS'));
                Xv=mon2ann(Xv_mon,mon_avg_o,mon_avg_f);
                Xv=permute(Xv,[2 1 3]);
                Xv(Xv==0)=NaN; % replace nans
                [rows,cols,len]=size(Xv(:,:,m_yrs));
                id_x0=size(X,1)+1;
                X=cat(1,X,reshape(Xv(:,:,m_yrs),rows*cols,len));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
		xmeta{j}.varnmlng='sea surface salinity';
		xmeta{j}.olat=olat;
		xmeta{j}.olon=olon;
		xmeta{j}.id_X=[id_x0,id_x1];
            case 'p'
                % Precipitation
                Xv_mon1=ncread(['./input-recon/pr/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.PRECC.085001-184912.nc'],'PRECC');
                Xv_mon2=ncread(['./input-recon/pr/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.PRECL.085001-184912.nc'],'PRECL');
                Xv_mon=Xv_mon1+Xv_mon2; % precip = convective + large scale

                Xv_mon=permute(Xv_mon,[2 1 3]);

		disp('Precip: Bias correcting precipitation to climatology')
		[X_clim] = load_obsclim('gpcp',lat,lon);
		[nt,nn,nm]=size(Xv_mon);
		% Convert precipitation from m/s to mm/day (1 day = 86400 sec)
		Xv_mon=Xv_mon*86400*1000;
		Xv_mon12=reshape(Xv_mon,nt,nn,12,nm/12);
		Xv_ns=bsxfun(@minus,Xv_mon12,mean(Xv_mon12,4));% remove climatology
		Xv_bc=bsxfun(@plus,Xv_ns,X_clim); % replace climatology
		% Replace zero or less with gamma random numbers with a mean of 0.05, a value much less than error in clim obs.
		%  In distribution plots, this gamma error makes the distributions at small values look very much like
		%  the original model distributions (when it has been converted into the appropriate units)
		smlidx=find(Xv_bc<=0);Xv_bc(smlidx)=gamrnd(1,0.05,length(smlidx),1);
		Xv_mon=reshape(Xv_bc,nt,nn,nm);	

		Xv=mon2ann(Xv_mon,mon_avg_o,mon_avg_f);


                [rows,cols,len]=size(Xv(:,:,m_yrs));
                id_x0=size(X,1)+1;
                X=cat(1,X,reshape(Xv(:,:,m_yrs),rows*cols,len));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
		xmeta{j}.varnmlng='precipitation';
		xmeta{j}.lat=lat;
		xmeta{j}.lon=lon;
		xmeta{j}.id_X=[id_x0,id_x1];
                
            case 'i'
                % SPEI 12
                disp(' !! Loading full, bias-corrected, SPEI field (NaNs included) !! ')
                load(['./input-recon/cesm_lme_' ens '_spei_output_biascorr_scl_12_krnl_e_u2_gcm.mat'],'spei_f')
                Xv=mon2ann(spei_f,mon_avg_o,mon_avg_f);
                Xv(Xv==0)=NaN;
                [rows,cols,len]=size(Xv(:,:,m_yrs));
                id_x0=size(X,1)+1;
                X=cat(1,X,reshape(Xv(:,:,m_yrs),rows*cols,len));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
		xmeta{j}.varnmlng='SPEI 12 month';
		xmeta{j}.lat=lat;
		xmeta{j}.lon=lon;
		xmeta{j}.id_X=[id_x0,id_x1];
                
            case 'r'
                disp(' !! Loading full, bias-corrected, PDSI field (NaNs included) !! ')
		load(['./input-recon/cesm_lme_' ens '_pdsi_output_biascorr_AWC_c_7_u2_gcm.mat'],'pdsi_f')             
   		%load(['./input-recon/cesm_lme_' ens '_pdsi_output_AWC_c_7_u2_gcm.mat'],'pdsi_f')
                Xv=mon2ann(pdsi_f,mon_avg_o,mon_avg_f);
                Xv(Xv==0)=NaN;
                [rows,cols,len]=size(Xv(:,:,m_yrs));
                id_x0=size(X,1)+1;
                X=cat(1,X,reshape(Xv(:,:,m_yrs),rows*cols,len));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
		xmeta{j}.varnmlng='PDSI';
		xmeta{j}.lat=lat;
		xmeta{j}.lon=lon;
		xmeta{j}.id_X=[id_x0,id_x1];
                
            case 'j'
                % SPEI 12 (land only -- for state vector...)
                
                %load(['./input-recon/cesm_lme_' ens '_spei_output_scl_12_krnl_e_u2_gcm.mat'],'spei_lndonly')
                disp('Bias-corrected SPEI')
		load(['./input-recon/cesm_lme_' ens '_spei_output_biascorr_scl_12_krnl_e_u2_gcm.mat'],'spei_lndonly','lndidx')
 		Xv=mon2ann(spei_lndonly,mon_avg_o,mon_avg_f);
                Xv(Xv==0)=NaN;
                id_x0=size(X,1)+1;
                X=cat(1,X,Xv(:,m_yrs));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
		xmeta{j}.varnmlng='SPEI 12 months (land only)';
		xmeta{j}.lat=lat;
		xmeta{j}.lon=lon;
		xmeta{j}.id_X=[id_x0,id_x1];
		xmeta{j}.lndidx=lndidx;
                
            case 'q'
                % PDSI (land-only--for state vector...)
                % single uniform available water content
                
                %load(['./input-recon/cesm_lme_' ens '_pdsi_output_AWC_c_7_u2_gcm.mat'],'pdsi_lndonly')
                disp('Bias-corrected PDSI')
		load(['./input-recon/cesm_lme_' ens '_pdsi_output_biascorr_AWC_c_7_u2_gcm.mat'],'pdsi_lndonly','lndidx')   
                Xv=mon2ann(pdsi_lndonly,mon_avg_o,mon_avg_f);
                Xv(Xv==0)=NaN;
                id_x0=size(X,1)+1;
                X=cat(1,X,Xv(:,m_yrs));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
		xmeta{j}.varnmlng='PDSI (land only)';
		xmeta{j}.lat=lat;
		xmeta{j}.lon=lon;
		xmeta{j}.id_X=[id_x0,id_x1];
		xmeta{j}.lndidx=lndidx;
                
            case 'o'
                % AMO (annual, using skin temperature on regular lat/lon grid)
                Xv_mon=ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TS.085001-184912.nc'],'TS');
               	
                Xv_mon=permute(Xv_mon,[2 1 3]);

		disp('AMO: Bias correcting surface temperature to climatology')
		[X_clim] = load_obsclim('bearth',lat,lon);
		[nt,nn,nm]=size(Xv_mon);
		Xv_mon12=reshape(Xv_mon,nt,nn,12,nm/12);
		Xv_ns=bsxfun(@minus,Xv_mon12,mean(Xv_mon12,4));% remove climatology
		Xv_bc=bsxfun(@plus,Xv_ns,X_clim); % replace climatology
		Xv_mon=reshape(Xv_bc,nt,nn,nm);	

		Xv=mon2ann(Xv_mon,mon_avg_o,mon_avg_f);

                load('./input-recon/north_atlantic_ocean_mask.mat')
                AMO = amo(Xv,atlntc_mask,lat,0);%0=no smoothing
                len=length(AMO(m_yrs));
                id_x0=size(X,1)+1;
                X=cat(1,X,reshape(AMO(m_yrs),1,len));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
		xmeta{j}.varnmlng='AMO (NASST index)';
		xmeta{j}.id_X=[id_x0,id_x1];

            case 'n'
                % Nino indices (monthly, using skin temperature on regular lat/lon grid)
                Xv_mon=ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TS.085001-184912.nc'],'TS');
                Xv_mon=permute(Xv_mon,[2 1 3]);
       
		disp('Nino indices: Bias correcting surface temperature to climatology')
		[X_clim] = load_obsclim('bearth',lat,lon);
		[nt,nn,nm]=size(Xv_mon);
		Xv_mon12=reshape(Xv_mon,nt,nn,12,nm/12);
		Xv_ns=bsxfun(@minus,Xv_mon12,mean(Xv_mon12,4));% remove climatology
		Xv_bc=bsxfun(@plus,Xv_ns,X_clim); % replace climatology
		Xv_mon=reshape(Xv_bc,nt,nn,nm);	

                % Compute all the nino indices
                %bc='y';disp('Bias-correcting Nino index');% option corrects index; field is already corrected above
                bc='n';
		nr=[];nmeta=cell(0);
                for ni=1:5
		   [nv,nmeta0]=nino(Xv_mon,lat,lon,ni,bc);
		   % Account for different year definitions
		   len=size(Xv_mon,3);
		   if mon_avg_o==4 && mon_avg_f==3
		      nv=nv(mon_avg_o:(len-(12-mon_avg_f)));
		   elseif mon_avg_o==1 && mon_avg_f==12
		      nv=nv(mon_avg_o:(len-(12-mon_avg_f)));
		   else 
		      disp('Just using tropical year for monthly Nino index for seasonal recons')
		      nv=nv(4:(len-(12-3)));
		   end
		   % Reshape index so that each year has 12 monthly values
		   n0=reshape(nv,12,length(nv)/12);
		   nr=cat(1,nr,n0);
		   nmeta{ni}=nmeta0;
                end   

                id_x0=size(X,1)+1;
                X=cat(1,X,nr(:,m_yrs));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
	        xmeta{j}.indinfo=nmeta; 	
		xmeta{j}.varnmlng='Nino indices';
		xmeta{j}.id_X=[id_x0,id_x1];


            case 'z'
                % ITCZ location (can be seasonal or annual)
                Xv_mon1=ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.PRECC.085001-184912.nc'],'PRECC');
                Xv_mon2=ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.PRECL.085001-184912.nc'],'PRECL');
                Xv_mon=Xv_mon1+Xv_mon2; % precip = convective + large scale
                Xv_mon=permute(Xv_mon,[2 1 3]);

		disp('ITCZ: Bias correcting precipitation to climatology')
		[X_clim] = load_obsclim('gpcp',lat,lon);
		[nt,nn,nm]=size(Xv_mon);
		% Convert precipitation from m/s to mm/day (1 day = 86400 sec)
		Xv_mon=Xv_mon*86400*1000;
		Xv_mon12=reshape(Xv_mon,nt,nn,12,nm/12);
		Xv_ns=bsxfun(@minus,Xv_mon12,mean(Xv_mon12,4));% remove climatology
		Xv_bc=bsxfun(@plus,Xv_ns,X_clim); % replace climatology
		% Replace zero or less with gamma random numbers with a mean of 0.05, a value much less than error in clim obs.
		%  In distribution plots, this gamma error makes the distributions at small values look very much like
		%  the original model distributions (when it has been converted into the appropriate units)
		smlidx=find(Xv_bc<=0);Xv_bc(smlidx)=gamrnd(1,0.05,length(smlidx),1);
		Xv_mon=reshape(Xv_bc,nt,nn,nm);	

		Xv=mon2ann(Xv_mon,mon_avg_o,mon_avg_f);
                
                
		% ITCZ location over the ocean regions
		
		ilons={[320 345],[130 260],[160 260],[130 170],[170 260],[65 95],[50 95],[260 320],[95 130],[345 50],[28 50]};
		ilnf={'Atlantic','Pacific','Pacific (Schneider 2014)','East Pacific','West Pacific','South Asian Monsoon (Schneider 2014)',...
			'Indian','South America','Indonesia','Africa','East Africa'};
		imeta=cell(0);ITCZ=[];
		for ii=1:length(ilons);
		   [itcz_mx,~,imeta0] = itcz_lon(Xv,lat,lon,ilons{ii});                   		
		   imeta{ii}=imeta0;
		   imeta{ii}.varnmlng=ilnf{ii};
		   ITCZ=cat(1,ITCZ,itcz_mx);
		end


                id_x0=size(X,1)+1;
                X=cat(1,X,ITCZ(:,m_yrs));
                id_x1=size(X,1);
                clear Xv Xv_mon
		xmeta{j}.varnm=state_tp(j);
	        xmeta{j}.indinfo=imeta; 	
		xmeta{j}.varnmlng='ITCZ indices';
		xmeta{j}.id_X=[id_x0,id_x1];

            otherwise
                error('Incorrect state specification.')
        end
    end

end
end

