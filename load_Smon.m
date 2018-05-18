function [X,lat,lon,id_X,olat,olon] = load_Smon(mdl_c,state_tp,mon_avg_o,mon_avg_f,m_yrs)
%LOAD_SMON  Load climate model or reanalysis data for state vector, monthly version.
%   [X,xlat,xlon,id_X] = load_Smon(model,state_tp,mon_avg_o,mon_avg_f,m_yrs)
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
%   Nathan Steiger, LDEO July 2017

X=[];id_X=zeros(2*length(state_tp),1);

if strcmp(mdl_c(1:8),'cesm_lme')
    
    ens=mdl_c(9:11);
    
    lat=ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TREFHT.085001-184912.nc'],'lat');
    lon=ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TREFHT.085001-184912.nc'],'lon');
    
    for j=1:length(state_tp)
        switch state_tp(j)
            case 's'
                % SST
                olat=permute(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SST.085001-184912.nc'],'TLAT'),[2 1]);
                olon=permute(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SST.085001-184912.nc'],'TLONG'),[2 1]);
                Xv_mon=squeeze(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SST.085001-184912.nc'],'SST'));
                Xv_mon=permute(Xv_mon,[2 1 3]);
      		% Pull out the months of the chosen years
		if mon_avg_o > mon_avg_f % need to pull out an extra year for the last year to average correctly
		   amons=bsxfun(@plus,([m_yrs m_yrs(end)+1]-1)*12+1,repmat(0:11,[length(m_yrs)+1,1])');
		else 
		   amons=bsxfun(@plus,(m_yrs-1)*12+1,repmat(0:11,[length(m_yrs),1])');
		end
		Xv=Xv_mon(:,:,amons(:));
                [rows,cols,len]=size(Xv);
                id_X(2*j-1)=size(X,1)+1;
                X=cat(1,X,reshape(Xv,rows*cols,len));
                id_X(2*j)=size(X,1);
            case 'y'
 		% Sea surface salinit(y)
                olat=permute(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SSS.085001-184912.nc'],'TLAT'),[2 1]);
                olon=permute(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SSS.085001-184912.nc'],'TLONG'),[2 1]);
                Xv_mon=squeeze(ncread(['./input-recon/b.e11.BLMTRC5CN.f19_g16.' ens '.pop.h.SSS.085001-184912.nc'],'SSS'));
                Xv_mon=permute(Xv_mon,[2 1 3]);
      		% Pull out the months of the chosen years
		if mon_avg_o > mon_avg_f % need to pull out an extra year for the last year to average correctly
		   amons=bsxfun(@plus,([m_yrs m_yrs(end)+1]-1)*12+1,repmat(0:11,[length(m_yrs)+1,1])');
		else 
		   amons=bsxfun(@plus,(m_yrs-1)*12+1,repmat(0:11,[length(m_yrs),1])');
		end
		Xv=Xv_mon(:,:,amons(:));
                [rows,cols,len]=size(Xv);
                id_X(2*j-1)=size(X,1)+1;
                X=cat(1,X,reshape(Xv,rows*cols,len));
                id_X(2*j)=size(X,1);
            otherwise
                error('Incorrect state specification.')
        end
    end
    
else
    error('Not available for another model...')
end


end

