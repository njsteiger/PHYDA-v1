function [gcm_syr,gcm_eyr]=modelyrs(mdl)
% model_info    Get's important info about model data

if strcmp(mdl,'ccsm4_lm')
  gcm_eyr=1850; gcm_syr=850;
elseif strcmp(mdl,'ccsm4_2mpi')
  gcm_eyr=1850; gcm_syr=850;
elseif strcmp(mdl,'ccsm4_2gfdlpic')
  gcm_eyr=1850; gcm_syr=850;
elseif strcmp(mdl,'ccsm4_pic')
  gcm_eyr=1051; gcm_syr=1; 
elseif strcmp(mdl,'mpi_lm')
  gcm_eyr=1849; gcm_syr=850;
elseif strcmp(mdl,'gfdlcm3_pic')
  gcm_eyr=800; gcm_syr=1;
elseif strcmp(mdl,'gfdlcm3_hist')
  gcm_eyr=2005; gcm_syr=1860;
elseif strcmp(mdl,'echam5')
  gcm_eyr=2011; gcm_syr=1871;
elseif strcmp(mdl,'echam5_2cam5iso')
    gcm_eyr=2011; gcm_syr=1871;
elseif strcmp(mdl,'speedy_lm')
    gcm_eyr=2004; gcm_syr=1000;
elseif strcmp(mdl,'icam5_hist')
    gcm_eyr=2014; gcm_syr=1850;
elseif strcmp(mdl,'cesm_lme002')
    gcm_eyr=1849; gcm_syr=850;
elseif strcmp(mdl,'cesm_lme003')
    gcm_eyr=1849; gcm_syr=850;
elseif strcmp(mdl,'cesm_lme004')
    gcm_eyr=1849; gcm_syr=850;
elseif strcmp(mdl,'cesm_lme005')
    gcm_eyr=1849; gcm_syr=850;
elseif strcmp(mdl,'cesm_lme006')
    gcm_eyr=1849; gcm_syr=850;
elseif strcmp(mdl,'cesm_lme007')
    gcm_eyr=1849; gcm_syr=850;
elseif strcmp(mdl,'cesm_lme008')
    gcm_eyr=1849; gcm_syr=850;
elseif strcmp(mdl,'cesm_lme009')
    gcm_eyr=1849; gcm_syr=850;
elseif strcmp(mdl,'cesm_lme010')
    gcm_eyr=1849; gcm_syr=850;
elseif strcmp(mdl,'gfdlesm2g_pi')
    gcm_eyr=500; gcm_syr=1;
elseif strcmp(mdl,'gfdlesm2m_pi')
    gcm_eyr=500; gcm_syr=1;
else
    error('unknown model type...')
end

