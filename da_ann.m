%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
%  DRIVING FUNCTION FOR ANNUAL-ONLY DA RECONSTRUCTIONS
%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear

% Load parameters and data (see file for details)
prmtrs_ann_da;

disp('Performing reconstruction...')

tic

%==========================================================================
%                          RECONSTRUCTION
%==========================================================================

% Optionally use 'matfile' to save large Xa ensemble variable to hard drive
if ens_save =='a'
    Xa_ens=zeros(svl,length(p_yrs),reconYrs);
elseif ens_save =='y';
    mf=matfile([outpath fln],'writable',true);
    mf.Xa_ens(svl,length(p_yrs),reconYrs)=0;
elseif ens_save =='s' % Save subsample of posterior
    Xa_ens=zeros(svl,sub_ens,reconYrs);
    se=randsample(length(p_yrs),sub_ens);
end

Xa_m=zeros(svl,reconYrs);
Xa_sigma=zeros(svl,reconYrs);
Xa_prct=zeros(svl,length(xa_prct),reconYrs);

% Extract proxy and proxy estimate information
HXb0=zeros(length(proxy),length(p_yrs));
y0=zeros(length(proxy),reconYrs);
R0v=zeros(length(proxy),length(clb_prd(1):clb_prd(2)));
for i=1:length(proxy)
   HXb0(i,:)=proxy{i}.HXb;
   R0v(i,:)=proxy{i}.Rvec; 
   % Extract the correct years of the proxy data
   ay=proxy{i}.age_yng;
   [~,~,ia] = intersect(r_o:r_f,ay);
   y0(i,:)=proxy{i}.data(ia);
end


%%%%%% Option of MC iterations commented out %%%%%%
% Loop over MC iterations
%for ii=1:itr
    
% Randomly sample 'swtchbld'% of proxies or use all proxies but random order
%rsp=randsample(length(proxy),round(length(proxy)*(swtchbld/100)))';
%disp('Using 100% of available proxies each MC iteration')      

% Loop over time
for t=1:reconYrs
	
   % Initialize Xb for each year (could also augment state vector here)
   Xb=Xb_o;

   % Use only non-NaN proxy values for the time step       
   ni=find(~isnan(y0(:,t)));
   % Extract proxies, HXb and R for each year
   if length(ni) > 0
      y=y0(ni,t);
      HXb=HXb0(ni,:);
      % Assume diagonal R
      R=diag(diag(nancov(R0v(ni,:)'))); 
      % Don't assume diagonal R --> Modify line 42 in M_update!! -->
      % results may be numerically unstable
      %R=nancov(R0v(ni,:)'); 
      [Xa]=M_update(Xb,HXb,y,R,infl);

   else
      % No proxies are assimilated at time 't'
      Xa=Xb;

   end

   % Compute statistics of full analysis each year
   Xa_m(:,t)=mean(Xa(1:svl,:),2);
   Xa_sigma(:,t)=std(Xa(1:svl,:),0,2);       
   Xa_prct(:,:,t)=prctile(Xa(1:svl,:),xa_prct,2);

   % Save full analysis
   if ens_save=='a'
      Xa_ens(:,:,t)=Xa(1:svl,:);
   elseif ens_save=='y'
      mf.Xa_ens(:,:,t)=Xa(1:svl,:);
   elseif ens_save=='s'
      Xa_ens(:,:,t)=Xa(1:svl,se);
   end
      
   % Note progress for reconstructions > 100 years
   if mod(t,100)==0; disp(['Year ' num2str(t) ' completed']); toc; end

end  

%disp(['DA iteration = ' num2str(ii) ' completed'])
%toc

%end

disp('DA reconstruction completed')
toc


%---------------------
% SAVE ALL VARIABLES
%---------------------
if outsave=='y'
   disp('Saving DA output...')
   if ens_save=='y'
     % append all variables to Xma and move to new directory
     allvars=whos;
     save([outpath fln],allvars.name,'-append') % use for matlab i/o file
   else
     save([outpath fln],'-v7.3')
   end
end

%toc
disp(' ')



















