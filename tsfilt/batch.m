flight=40030;
inst=1;
%% do the filtering
%for ifield=1:8
%    filt_dark_flight(flight,inst,ifield);
%end
%% filter dark long and get DC template
%get_filt_dctemplate(flight,inst);
%% get bigmask, FF, Fweight
get_bigmask(flight,inst);
get_FFflight(flight,inst);
get_fweight(flight,inst);
get_linearized_map(flight,inst);
stack_preprocess(flight,inst);

% mkk & xmkk needs rerun with new mask
get_mkk(flight,inst);
get_xmkk(flight,inst);
%% plot the 2D PS before and after filtering (can skip)
%plot_filt_2DPS(flight,inst,8);
%% get the diff power spectrum and plot wrt ell and nfr
get_raw_filt_Cl(flight,inst);
plot_Cl_l(flight,inst);
plot_Cl_nfr(flight,inst);
%% Sims
%get_darkstat(flight,inst);
%rnrealization_valid(flight,inst,ifield);% can skip
%sim_nbias(flight,inst,G1);
%% fit G1 from diff
get_chi2_fitg1diff(flight,inst);
plot_chi2g1fit(flight,inst);
plot_Cl_ph(flight,inst);

%get_chi2_fitg1dfull(flight,inst); %% code not done yet
%%
get_bl(flight,inst);

%G1=-3.6;
%get_flight_PS(flight,inst,G1);
%plot_flight_auto(flight,inst);