function [err, timepoints, species_out, observables_out] = TSP1_model_withVEGF_V3_resubmission( timepoints, species_init, parameters, suppress_plot )
%TSP1_MODEL_WITHVEGF_V3_RESUBMISSION Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'TSP1_model_withVEGF_V3_resubmission' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. TSP1_MODEL_WITHVEGF_V3_RESUBMISSION returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = TSP1_model_withVEGF_V3_resubmission( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of 132 initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of 88 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 0.51931, 0.00027, 0.002446, 7.5e-7, 1.3e-5, 1.3e-5, 1534, 105, 9.97E-06, 1.00E-05, 1.00E-05, 6.02E+23, 10, 10, 330, 330, 0.387, 0.0324, 0.00028, 5.00E+05, 5.00E+05, 5.00E+05, 5.00E+05, 8600, 2.10E+05, 5.00E+05, 10000, 1.00E+05, 3.20E+06, 1.00E+06, 3.00E+07, 1.00E+07, 3.00E+07, 1.00E+07, 0.00049, 5.30E-08, 1.00E+14, 3.10E+13, 1.00E+14, 0.00028, 0.115, 0.00069, 0.1, 0.005, 0.005, 0.0025, 0.05, 0.001, 0.0022303, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.01, 0.000193, 0.00033, 0.0012, 0.00386, 0.0019, 3.10E+11, 3.10E+13, 3.10E+11, 631, 2500, 10000, 5000, 10000, 1250, 5000, 2500, 5000, 1100, 550, 39500, 39500, 3750, 300, 20000, 20000 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 88  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 88].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 132  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 132].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,10*3600,600)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  ||  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  ||  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'ECM_Vol_tis_dis', 'EBM_Vol_tis_dis', 'PBM_Vol_tis_dis', 'ECM_conc', 'EBM_conc', 'PBM_conc', 'tumorSA_Vol_tis_dis', 'VesselSA_Vol_tis_dis', 'tumorCellSurfArea_tis_dis', 'VesselCellSurfArea_blood', 'VesselCellSurfArea_tis_dis', 'Avogadro', 'ratio', 'qTSP1', 'qMMP3', 'qproMMP9', 'qV165_tumor', 'qV165_disEC', 'sR_receptors', 'kon_TSP1_CD36', 'kon_TSP1_GAG', 'kon_TSP1_CD47', 'kon_TSP1_VEGF', 'kon_V165_GAG', 'kon_TSP1_LRP1', 'kon_TSP1_B1', 'kon_MMP3_proMMP9', 'kon_TSP1_MMP3', 'kon_V165_N1', 'kon_V165_N2', 'kon_V165_R1', 'kon_V165_R2', 'kon_V121_R1', 'kon_V121_R2', 'koff_MMP9_LRP1', 'Kd_MMP9_LRP1', 'kc_V165N_R2_tumor', 'kc_V165R2_N_tumor', 'kc_R1_N_tumor', 'k_int_receptors', 'koff_TSP1_CD36', 'koff_V165_GAG', 'koff_TSP1_GAG', 'koff_TSP1_CD47', 'koff_TSP1_VEGF', 'koff_TSP1_LRP1', 'koff_TSP1_B1', 'koff_MMP3_proMMP9', 'koff_MMP3_TSP1', 'koff_V165N1_R2', 'koff_V121_R1', 'koff_V165R2_N1', 'koff_V121_R2', 'koff_V165N2_R2', 'koff_V165R2_N2', 'koff_V165_N1', 'koff_V165_N2', 'koff_V165_R1', 'koff_V165_R2', 'kdissoc_CD36_B1', 'kdissoc_CD47_VEGFR2', 'kdissoc_CD36_VEGFR2', 'kdissoc_R1_N', 'kdeg_VEGF', 'kdeg_TSP1', 'kdeg_MMP', 'k_TSP1cleave', 'k_act_MMP3_proMMP9', 'kcCD36_R2', 'kcCD36_B1', 'kcCD47_R2', 'kp_mmp', 'CD36_number_tum', 'CD47_number_tum', 'LRP1_number_tum', 'B1_number_tum', 'CD36_number_tis_dis_A', 'CD47_number_tis_dis_A', 'LRP1_number_tis_dis_A', 'B1_number_tis_dis_A', 'R1_number_tum', 'R2_number_tum', 'N1_number_tum', 'N2_number_tum', 'R1_number_tis_dis_A', 'R2_number_tis_dis_A', 'N1_number_tis_dis_A', 'N2_number_tis_dis_A' };



%% Integrate Network Model
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
opts = odeset( 'RelTol',   1e-12,   ...
               'AbsTol',   1e-18,   ...
               'Stats',    'off',  ...
               'BDF',      'off',    ...
               'MaxOrder', 5   );


% define derivative function
rhs_fcn = @(t,y)( calc_species_deriv( t, y, expressions ) );

% simulate model system (stiff integrator)
try 
    [~, species_out] = ode15s( rhs_fcn, timepoints, species_init', opts );
    if(length(timepoints) ~= size(species_out,1))
        exception = MException('ODE15sError:MissingOutput','Not all timepoints output\n');
        throw(exception);
    end
catch
    err = 1;
    fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
    return;
end

% calculate observables
observables_out = zeros( length(timepoints), 32 );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'TSP1_free', 'TSP1_VEGF', 'VEGF', 'MMP3', 'MMP9', 'proMMP9', 'MMP3_TSP1', 'V165', 'V121', 'V114', 'TSP1_cleaved', 'TSP1_total', 'TSP1_VEGFbound', 'TSP1_RECbound', 'TSP1_ECMbound', 'TSP1_MMP3bound', 'VEGF_total', 'VEGF_RECbound', 'VEGF_ECMbound', 'TSP1_CD36', 'TSP1_CD47', 'TSP1_LRP1', 'TSP1_B1', 'TSP1_B1_CD36', 'VEGF_R1', 'VEGF_R2', 'VEGF_N12', 'VEGF_R1_N1N2', 'VEGF_R2_N1N2', 'VEGF_R2_other', 'Total_REC_TSP1', 'Total_REC_VEGF' };

    % construct figure
    plot(timepoints,observables_out);
    title('TSP1_model_withVEGF_V3_resubmission observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number or concentration','fontSize',12,'Interpreter','none');

end


%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%


% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,132);
    species_init(1) = 0;
    species_init(2) = 1;
    species_init(3) = (params(4)*params(1))/1000;
    species_init(4) = (params(5)*params(2))/1000;
    species_init(5) = (params(6)*params(3))/1000;
    species_init(6) = ((params(73)/params(9))/params(12))*params(7);
    species_init(7) = ((params(76)/params(9))/params(12))*params(7);
    species_init(8) = ((params(82)/params(9))/params(12))*params(7);
    species_init(9) = ((params(74)/params(9))/params(12))*params(7);
    species_init(10) = ((params(75)/params(9))/params(12))*params(7);
    species_init(11) = ((params(77)/params(11))/params(12))*params(8);
    species_init(12) = ((params(80)/params(11))/params(12))*params(8);
    species_init(13) = ((params(86)/params(11))/params(12))*params(8);
    species_init(14) = ((params(78)/params(11))/params(12))*params(8);
    species_init(15) = ((params(79)/params(11))/params(12))*params(8);
    species_init(16) = ((params(81)/params(9))/params(12))*params(7);
    species_init(17) = ((params(83)/params(9))/params(12))*params(7);
    species_init(18) = ((params(84)/params(9))/params(12))*params(7);
    species_init(19) = ((params(87)/params(11))/params(12))*params(8);
    species_init(20) = ((params(85)/params(11))/params(12))*params(8);
    species_init(21) = ((params(88)/params(11))/params(12))*params(8);
    species_init(22) = 0;
    species_init(23) = 0;
    species_init(24) = 0;
    species_init(25) = 0;
    species_init(26) = 0;
    species_init(27) = 0;
    species_init(28) = 0;
    species_init(29) = 0;
    species_init(30) = 0;
    species_init(31) = 0;
    species_init(32) = 0;
    species_init(33) = 0;
    species_init(34) = 0;
    species_init(35) = 0;
    species_init(36) = 0;
    species_init(37) = 0;
    species_init(38) = 0;
    species_init(39) = 0;
    species_init(40) = 0;
    species_init(41) = 0;
    species_init(42) = 0;
    species_init(43) = 0;
    species_init(44) = 0;
    species_init(45) = 0;
    species_init(46) = 0;
    species_init(47) = 0;
    species_init(48) = 0;
    species_init(49) = 0;
    species_init(50) = 0;
    species_init(51) = 0;
    species_init(52) = 0;
    species_init(53) = 0;
    species_init(54) = 0;
    species_init(55) = 0;
    species_init(56) = 0;
    species_init(57) = 0;
    species_init(58) = 0;
    species_init(59) = 0;
    species_init(60) = 0;
    species_init(61) = 0;
    species_init(62) = 0;
    species_init(63) = 0;
    species_init(64) = 0;
    species_init(65) = 0;
    species_init(66) = 0;
    species_init(67) = 0;
    species_init(68) = 0;
    species_init(69) = 0;
    species_init(70) = 0;
    species_init(71) = 0;
    species_init(72) = 0;
    species_init(73) = 0;
    species_init(74) = 0;
    species_init(75) = 0;
    species_init(76) = 0;
    species_init(77) = 0;
    species_init(78) = 0;
    species_init(79) = 0;
    species_init(80) = 0;
    species_init(81) = 0;
    species_init(82) = 0;
    species_init(83) = 0;
    species_init(84) = 0;
    species_init(85) = 0;
    species_init(86) = 0;
    species_init(87) = 0;
    species_init(88) = 0;
    species_init(89) = 0;
    species_init(90) = 0;
    species_init(91) = 0;
    species_init(92) = 0;
    species_init(93) = 0;
    species_init(94) = 0;
    species_init(95) = 0;
    species_init(96) = 0;
    species_init(97) = 0;
    species_init(98) = 0;
    species_init(99) = 0;
    species_init(100) = 0;
    species_init(101) = 0;
    species_init(102) = 0;
    species_init(103) = 0;
    species_init(104) = 0;
    species_init(105) = 0;
    species_init(106) = 0;
    species_init(107) = 0;
    species_init(108) = 0;
    species_init(109) = 0;
    species_init(110) = 0;
    species_init(111) = 0;
    species_init(112) = 0;
    species_init(113) = 0;
    species_init(114) = 0;
    species_init(115) = 0;
    species_init(116) = 0;
    species_init(117) = 0;
    species_init(118) = 0;
    species_init(119) = 0;
    species_init(120) = 0;
    species_init(121) = 0;
    species_init(122) = 0;
    species_init(123) = 0;
    species_init(124) = 0;
    species_init(125) = 0;
    species_init(126) = 0;
    species_init(127) = 0;
    species_init(128) = 0;
    species_init(129) = 0;
    species_init(130) = 0;
    species_init(131) = 0;
    species_init(132) = 0;

end


% user-defined functions



% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,164);
    expressions(1) = parameters(1);
    expressions(2) = parameters(2);
    expressions(3) = parameters(3);
    expressions(4) = parameters(4);
    expressions(5) = parameters(5);
    expressions(6) = parameters(6);
    expressions(7) = parameters(7);
    expressions(8) = parameters(8);
    expressions(9) = parameters(9);
    expressions(10) = parameters(10);
    expressions(11) = parameters(11);
    expressions(12) = parameters(12);
    expressions(13) = parameters(13);
    expressions(14) = parameters(14);
    expressions(15) = parameters(15);
    expressions(16) = parameters(16);
    expressions(17) = parameters(17);
    expressions(18) = parameters(18);
    expressions(19) = parameters(19);
    expressions(20) = parameters(20);
    expressions(21) = parameters(21);
    expressions(22) = parameters(22);
    expressions(23) = parameters(23);
    expressions(24) = parameters(24);
    expressions(25) = parameters(25);
    expressions(26) = parameters(26);
    expressions(27) = parameters(27);
    expressions(28) = parameters(28);
    expressions(29) = parameters(29);
    expressions(30) = parameters(30);
    expressions(31) = parameters(31);
    expressions(32) = parameters(32);
    expressions(33) = parameters(33);
    expressions(34) = parameters(34);
    expressions(35) = parameters(35);
    expressions(36) = parameters(36);
    expressions(37) = parameters(37);
    expressions(38) = parameters(38);
    expressions(39) = parameters(39);
    expressions(40) = parameters(40);
    expressions(41) = parameters(41);
    expressions(42) = parameters(42);
    expressions(43) = parameters(43);
    expressions(44) = parameters(44);
    expressions(45) = parameters(45);
    expressions(46) = parameters(46);
    expressions(47) = parameters(47);
    expressions(48) = parameters(48);
    expressions(49) = parameters(49);
    expressions(50) = parameters(50);
    expressions(51) = parameters(51);
    expressions(52) = parameters(52);
    expressions(53) = parameters(53);
    expressions(54) = parameters(54);
    expressions(55) = parameters(55);
    expressions(56) = parameters(56);
    expressions(57) = parameters(57);
    expressions(58) = parameters(58);
    expressions(59) = parameters(59);
    expressions(60) = parameters(60);
    expressions(61) = parameters(61);
    expressions(62) = parameters(62);
    expressions(63) = parameters(63);
    expressions(64) = parameters(64);
    expressions(65) = parameters(65);
    expressions(66) = parameters(66);
    expressions(67) = parameters(67);
    expressions(68) = parameters(68);
    expressions(69) = parameters(69);
    expressions(70) = parameters(70);
    expressions(71) = parameters(71);
    expressions(72) = parameters(72);
    expressions(73) = parameters(73);
    expressions(74) = parameters(74);
    expressions(75) = parameters(75);
    expressions(76) = parameters(76);
    expressions(77) = parameters(77);
    expressions(78) = parameters(78);
    expressions(79) = parameters(79);
    expressions(80) = parameters(80);
    expressions(81) = parameters(81);
    expressions(82) = parameters(82);
    expressions(83) = parameters(83);
    expressions(84) = parameters(84);
    expressions(85) = parameters(85);
    expressions(86) = parameters(86);
    expressions(87) = parameters(87);
    expressions(88) = parameters(88);
    expressions(89) = ((expressions(1)+expressions(2))+expressions(3));
    expressions(90) = ((((expressions(19)*expressions(73))/expressions(9))/expressions(12))*expressions(7));
    expressions(91) = ((((expressions(19)*expressions(76))/expressions(9))/expressions(12))*expressions(7));
    expressions(92) = ((((expressions(19)*expressions(82))/expressions(9))/expressions(12))*expressions(7));
    expressions(93) = ((((expressions(19)*expressions(74))/expressions(9))/expressions(12))*expressions(7));
    expressions(94) = ((((expressions(19)*expressions(75))/expressions(9))/expressions(12))*expressions(7));
    expressions(95) = ((((expressions(19)*expressions(83))/expressions(9))/expressions(12))*expressions(7));
    expressions(96) = ((((expressions(19)*expressions(84))/expressions(9))/expressions(12))*expressions(7));
    expressions(97) = ((((expressions(19)*expressions(81))/expressions(9))/expressions(12))*expressions(7));
    expressions(98) = ((((expressions(19)*expressions(80))/expressions(11))/expressions(12))*expressions(8));
    expressions(99) = ((((expressions(19)*expressions(77))/expressions(11))/expressions(12))*expressions(8));
    expressions(100) = ((((expressions(19)*expressions(78))/expressions(11))/expressions(12))*expressions(8));
    expressions(101) = ((((expressions(19)*expressions(79))/expressions(11))/expressions(12))*expressions(8));
    expressions(102) = ((((expressions(19)*expressions(85))/expressions(11))/expressions(12))*expressions(8));
    expressions(103) = ((((expressions(19)*expressions(87))/expressions(11))/expressions(12))*expressions(8));
    expressions(104) = ((((expressions(19)*expressions(86))/expressions(11))/expressions(12))*expressions(8));
    expressions(105) = ((((expressions(19)*expressions(88))/expressions(11))/expressions(12))*expressions(8));
    expressions(106) = ((expressions(20)*1000)/expressions(89));
    expressions(107) = ((expressions(22)*1000)/expressions(89));
    expressions(108) = ((expressions(21)*1000)/expressions(89));
    expressions(109) = ((expressions(23)*1000)/expressions(89));
    expressions(110) = ((expressions(24)*1000)/expressions(89));
    expressions(111) = ((expressions(25)*1000)/expressions(89));
    expressions(112) = ((expressions(26)*1000)/expressions(89));
    expressions(113) = ((expressions(28)*1000)/expressions(89));
    expressions(114) = ((expressions(27)*1000)/expressions(89));
    expressions(115) = ((expressions(29)*1000)/expressions(89));
    expressions(116) = ((expressions(30)*1000)/expressions(89));
    expressions(117) = ((expressions(31)*1000)/expressions(89));
    expressions(118) = ((expressions(32)*1000)/expressions(89));
    expressions(119) = ((expressions(33)*1000)/expressions(89));
    expressions(120) = ((expressions(34)*1000)/expressions(89));
    expressions(121) = (((expressions(35)/expressions(36))*1000)/expressions(89));
    expressions(122) = (expressions(70)/expressions(7));
    expressions(123) = (expressions(71)/expressions(7));
    expressions(124) = (expressions(69)/expressions(7));
    expressions(125) = (expressions(39)/expressions(7));
    expressions(126) = (expressions(37)/expressions(7));
    expressions(127) = (expressions(38)/expressions(7));
    expressions(128) = (expressions(70)/expressions(8));
    expressions(129) = (expressions(69)/expressions(8));
    expressions(130) = (expressions(71)/expressions(8));
    expressions(131) = (expressions(37)/expressions(8));
    expressions(132) = (expressions(38)/expressions(8));
    expressions(133) = (expressions(39)/expressions(8));
    expressions(134) = (((expressions(16)/expressions(12))/expressions(11))*expressions(8));
    expressions(135) = (((expressions(14)/expressions(12))/expressions(11))*expressions(8));
    expressions(136) = (((expressions(17)/expressions(12))/expressions(9))*expressions(7));
    expressions(137) = (((((expressions(18)*10)/90)/expressions(12))/expressions(11))*expressions(8));
    expressions(138) = (((expressions(18)/expressions(12))/expressions(11))*expressions(8));
    expressions(139) = (((expressions(15)/expressions(12))/expressions(11))*expressions(8));
    expressions(140) = (((expressions(17)/expressions(12))/expressions(9))*expressions(7));
    expressions(141) = ((((expressions(14)/expressions(13))/expressions(12))/expressions(9))*expressions(7));
    expressions(142) = ((expressions(72)*1000)/expressions(89));
    expressions(143) = ((expressions(4)*expressions(1))/1000);
    expressions(144) = ((expressions(5)*expressions(2))/1000);
    expressions(145) = ((expressions(6)*expressions(3))/1000);
    expressions(146) = (((expressions(73)/expressions(9))/expressions(12))*expressions(7));
    expressions(147) = (((expressions(76)/expressions(9))/expressions(12))*expressions(7));
    expressions(148) = (((expressions(82)/expressions(9))/expressions(12))*expressions(7));
    expressions(149) = (((expressions(74)/expressions(9))/expressions(12))*expressions(7));
    expressions(150) = (((expressions(75)/expressions(9))/expressions(12))*expressions(7));
    expressions(151) = (((expressions(77)/expressions(11))/expressions(12))*expressions(8));
    expressions(152) = (((expressions(80)/expressions(11))/expressions(12))*expressions(8));
    expressions(153) = (((expressions(86)/expressions(11))/expressions(12))*expressions(8));
    expressions(154) = (((expressions(78)/expressions(11))/expressions(12))*expressions(8));
    expressions(155) = (((expressions(79)/expressions(11))/expressions(12))*expressions(8));
    expressions(156) = (((expressions(81)/expressions(9))/expressions(12))*expressions(7));
    expressions(157) = (((expressions(84)/expressions(9))/expressions(12))*expressions(7));
    expressions(158) = (((expressions(83)/expressions(9))/expressions(12))*expressions(7));
    expressions(159) = (((expressions(87)/expressions(11))/expressions(12))*expressions(8));
    expressions(160) = (((expressions(85)/expressions(11))/expressions(12))*expressions(8));
    expressions(161) = (((expressions(88)/expressions(11))/expressions(12))*expressions(8));
    expressions(162) = (expressions(141)+expressions(135));
    expressions(163) = (expressions(136)+expressions(138));
    expressions(164) = (expressions(140)+expressions(137));
   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,32);
    observables(1) = species(22);
    observables(2) = species(30) +species(31) +species(32);
    observables(3) = species(23) +species(24) +species(25);
    observables(4) = species(28);
    observables(5) = species(27);
    observables(6) = species(26);
    observables(7) = species(36);
    observables(8) = species(23);
    observables(9) = species(24);
    observables(10) = species(25);
    observables(11) = species(29);
    observables(12) = species(22) +species(29) +species(30) +species(31) +species(32) +species(33) +species(34) +species(36) +species(39) +species(40) +species(41) +species(61) +species(62) +species(65) +species(66) +species(67) +species(68) +species(73) +species(74) +species(75) +species(76) +species(77) +species(78) +species(95) +species(96) +species(98) +species(99) +species(101) +species(102) +species(103) +species(104) +species(105) +species(106) +species(107) +species(108) +species(110) +species(112) +species(113) +species(116) +species(118) +species(121) +species(123) +species(124) +species(125) +species(126) +species(127) +species(128) +species(129) +species(130) +species(131) +species(132);
    observables(13) = species(30) +species(31) +species(32) +species(73) +species(74) +species(75) +species(76) +species(77) +species(78);
    observables(14) = species(33) +species(34) +species(61) +species(62) +species(65) +species(66) +species(67) +species(68) +species(95) +species(96) +species(98) +species(99) +species(101) +species(102) +species(103) +species(104) +species(105) +species(106) +species(107) +species(108) +species(110) +species(112) +species(113) +species(116) +species(118) +species(121) +species(123) +species(124) +species(125) +species(126) +species(127) +species(128) +species(129) +species(130) +species(131) +species(132);
    observables(15) = species(39) +species(40) +species(41);
    observables(16) = species(36);
    observables(17) = species(23) +species(24) +species(25) +species(30) +species(31) +species(32) +species(42) +species(43) +species(44) +species(45) +species(46) +species(47) +species(48) +species(49) +species(50) +species(51) +species(52) +species(53) +species(54) +species(55) +species(56) +species(73) +species(74) +species(75) +species(76) +species(77) +species(78) +species(79) +species(80) +species(81) +species(82) +species(83) +species(84) +species(85) +species(86) +species(87) +species(88) +species(89) +species(90) +species(91) +species(92) +species(93) +species(94) +species(97) +species(98) +species(99) +species(100) +species(101) +species(102) +species(113) +species(114) +species(115) +species(116) +species(117) +species(118) +species(119) +species(120) +species(121) +species(122) +species(123) +species(124) +species(127) +species(128) +species(129) +species(130) +species(131) +species(132);
    observables(18) = species(45) +species(46) +species(47) +species(48) +species(49) +species(50) +species(51) +species(52) +species(53) +species(54) +species(55) +species(56) +species(79) +species(80) +species(81) +species(82) +species(83) +species(84) +species(85) +species(86) +species(87) +species(88) +species(89) +species(90) +species(91) +species(92) +species(93) +species(94) +species(97) +species(98) +species(99) +species(100) +species(101) +species(102) +species(113) +species(114) +species(115) +species(116) +species(117) +species(118) +species(119) +species(120) +species(121) +species(122) +species(123) +species(124) +species(127) +species(128) +species(129) +species(130) +species(131) +species(132);
    observables(19) = species(42) +species(43) +species(44);
    observables(20) = species(65) +species(66) +species(103) +species(104) +species(113) +species(116) +species(118) +species(121);
    observables(21) = species(61) +species(62) +species(95) +species(96) +species(98) +species(99) +species(101) +species(102);
    observables(22) = species(33) +species(34) +species(73) +species(74) +species(75) +species(76) +species(77) +species(78);
    observables(23) = species(67) +species(68);
    observables(24) = species(105) +species(106) +species(107) +species(108) +species(110) +species(112) +species(123) +species(124) +species(125) +species(126) +species(127) +species(128) +species(129) +species(130) +species(131) +species(132);
    observables(25) = species(45) +species(46) +species(53) +species(54);
    observables(26) = species(47) +species(48) +species(55) +species(56);
    observables(27) = species(49) +species(50) +species(51) +species(52);
    observables(28) = species(91) +species(92) +species(93) +species(94);
    observables(29) = species(83) +species(84) +species(85) +species(86) +species(87) +species(88) +species(89) +species(90);
    observables(30) = species(47) +species(48) +species(55) +species(56) +species(79) +species(80) +species(81) +species(82) +species(97) +species(98) +species(99) +species(100) +species(101) +species(102) +species(113) +species(114) +species(115) +species(116) +species(117) +species(118) +species(119) +species(120) +species(121) +species(122) +species(123) +species(124) +species(127) +species(128) +species(129) +species(130) +species(131) +species(132);
    observables(31) = species(33) +species(34) +species(37) +species(38) +species(61) +species(62) +species(65) +species(66) +species(67) +species(68) +species(73) +species(74) +species(75) +species(76) +species(77) +species(78) +species(95) +species(96) +species(98) +species(99) +species(101) +species(102) +species(103) +species(104) +species(105) +species(106) +species(107) +species(108) +species(110) +species(112) +species(113) +species(116) +species(118) +species(121) +species(123) +species(124) +species(125) +species(126) +species(127) +species(128) +species(129) +species(130) +species(131) +species(132);
    observables(32) = species(45) +species(46) +species(47) +species(48) +species(49) +species(50) +species(51) +species(52) +species(53) +species(54) +species(55) +species(56) +species(83) +species(84) +species(85) +species(86) +species(87) +species(88) +species(89) +species(90) +species(91) +species(92) +species(93) +species(94);

end


% Calculate ratelaws
function [ ratelaws ] = calc_ratelaws ( species, expressions, observables )

    ratelaws = zeros(1,32);
    ratelaws(1) = ((((expressions(19)*expressions(73))/expressions(9))/expressions(12))*expressions(7))*species(2);
    ratelaws(2) = ((((expressions(19)*expressions(77))/expressions(11))/expressions(12))*expressions(8))*species(2);
    ratelaws(3) = ((((expressions(19)*expressions(76))/expressions(9))/expressions(12))*expressions(7))*species(2);
    ratelaws(4) = ((((expressions(19)*expressions(80))/expressions(11))/expressions(12))*expressions(8))*species(2);
    ratelaws(5) = ((((expressions(19)*expressions(74))/expressions(9))/expressions(12))*expressions(7))*species(2);
    ratelaws(6) = ((((expressions(19)*expressions(78))/expressions(11))/expressions(12))*expressions(8))*species(2);
    ratelaws(7) = ((((expressions(19)*expressions(75))/expressions(9))/expressions(12))*expressions(7))*species(2);
    ratelaws(8) = ((((expressions(19)*expressions(79))/expressions(11))/expressions(12))*expressions(8))*species(2);
    ratelaws(9) = ((((expressions(19)*expressions(81))/expressions(9))/expressions(12))*expressions(7))*species(2);
    ratelaws(10) = ((((expressions(19)*expressions(82))/expressions(9))/expressions(12))*expressions(7))*species(2);
    ratelaws(11) = ((((expressions(19)*expressions(83))/expressions(9))/expressions(12))*expressions(7))*species(2);
    ratelaws(12) = ((((expressions(19)*expressions(84))/expressions(9))/expressions(12))*expressions(7))*species(2);
    ratelaws(13) = ((((expressions(19)*expressions(85))/expressions(11))/expressions(12))*expressions(8))*species(2);
    ratelaws(14) = ((((expressions(19)*expressions(86))/expressions(11))/expressions(12))*expressions(8))*species(2);
    ratelaws(15) = ((((expressions(19)*expressions(87))/expressions(11))/expressions(12))*expressions(8))*species(2);
    ratelaws(16) = ((((expressions(19)*expressions(88))/expressions(11))/expressions(12))*expressions(8))*species(2);
    ratelaws(17) = expressions(40)*species(16);
    ratelaws(18) = expressions(40)*species(20);
    ratelaws(19) = expressions(40)*species(8);
    ratelaws(20) = expressions(40)*species(13);
    ratelaws(21) = expressions(40)*species(18);
    ratelaws(22) = expressions(40)*species(21);
    ratelaws(23) = expressions(40)*species(17);
    ratelaws(24) = expressions(40)*species(19);
    ratelaws(25) = expressions(40)*species(6);
    ratelaws(26) = expressions(40)*species(11);
    ratelaws(27) = expressions(40)*species(7);
    ratelaws(28) = expressions(40)*species(12);
    ratelaws(29) = expressions(40)*species(9);
    ratelaws(30) = expressions(40)*species(14);
    ratelaws(31) = expressions(40)*species(10);
    ratelaws(32) = expressions(40)*species(15);
    ratelaws(33) = (expressions(141)+expressions(135))*species(2);
    ratelaws(34) = expressions(65)*species(22);
    ratelaws(35) = expressions(67)*species(22);
    ratelaws(36) = expressions(65)*species(29);
    ratelaws(37) = (expressions(136)+expressions(138))*species(2);
    ratelaws(38) = expressions(64)*species(23);
    ratelaws(39) = (expressions(140)+expressions(137))*species(2);
    ratelaws(40) = expressions(64)*species(24);
    ratelaws(41) = expressions(64)*species(25);
    ratelaws(42) = ((expressions(23)*1000)/expressions(89))*species(22)*species(23);
    ratelaws(43) = ((expressions(23)*1000)/expressions(89))*species(22)*species(24);
    ratelaws(44) = ((expressions(23)*1000)/expressions(89))*species(22)*species(25);
    ratelaws(45) = ((expressions(25)*1000)/expressions(89))*species(22)*species(10);
    ratelaws(46) = ((expressions(25)*1000)/expressions(89))*species(22)*species(15);
    ratelaws(47) = ((expressions(27)*1000)/expressions(89))*species(26)*species(28);
    ratelaws(48) = ((expressions(28)*1000)/expressions(89))*species(28)*species(22);
    ratelaws(49) = (((expressions(16)/expressions(12))/expressions(11))*expressions(8))*species(2);
    ratelaws(50) = (((expressions(15)/expressions(12))/expressions(11))*expressions(8))*species(2);
    ratelaws(51) = expressions(66)*species(28);
    ratelaws(52) = expressions(66)*species(27);
    ratelaws(53) = expressions(66)*species(26);
    ratelaws(54) = ((expressions(72)*1000)/expressions(89))*species(28)*species(23);
    ratelaws(55) = ((expressions(72)*1000)/expressions(89))*species(27)*species(23);
    ratelaws(56) = (((expressions(35)/expressions(36))*1000)/expressions(89))*species(27)*species(10);
    ratelaws(57) = (((expressions(35)/expressions(36))*1000)/expressions(89))*species(27)*species(15);
    ratelaws(58) = ((expressions(21)*1000)/expressions(89))*species(22)*species(3);
    ratelaws(59) = ((expressions(21)*1000)/expressions(89))*species(22)*species(4);
    ratelaws(60) = ((expressions(21)*1000)/expressions(89))*species(22)*species(5);
    ratelaws(61) = ((expressions(24)*1000)/expressions(89))*species(23)*species(3);
    ratelaws(62) = ((expressions(24)*1000)/expressions(89))*species(23)*species(4);
    ratelaws(63) = ((expressions(24)*1000)/expressions(89))*species(23)*species(5);
    ratelaws(64) = ((expressions(31)*1000)/expressions(89))*species(23)*species(16);
    ratelaws(65) = ((expressions(31)*1000)/expressions(89))*species(23)*species(20);
    ratelaws(66) = ((expressions(32)*1000)/expressions(89))*species(23)*species(8);
    ratelaws(67) = ((expressions(32)*1000)/expressions(89))*species(23)*species(13);
    ratelaws(68) = ((expressions(29)*1000)/expressions(89))*species(23)*species(17);
    ratelaws(69) = ((expressions(29)*1000)/expressions(89))*species(23)*species(19);
    ratelaws(70) = ((expressions(30)*1000)/expressions(89))*species(23)*species(18);
    ratelaws(71) = ((expressions(30)*1000)/expressions(89))*species(23)*species(21);
    ratelaws(72) = ((expressions(33)*1000)/expressions(89))*species(24)*species(16);
    ratelaws(73) = ((expressions(33)*1000)/expressions(89))*species(24)*species(20);
    ratelaws(74) = ((expressions(34)*1000)/expressions(89))*species(24)*species(8);
    ratelaws(75) = ((expressions(34)*1000)/expressions(89))*species(24)*species(13);
    ratelaws(76) = (expressions(39)/expressions(7))*species(16)*species(18);
    ratelaws(77) = (expressions(39)/expressions(7))*species(16)*species(17);
    ratelaws(78) = (expressions(39)/expressions(8))*species(20)*species(19);
    ratelaws(79) = (expressions(39)/expressions(8))*species(20)*species(21);
    ratelaws(80) = ((expressions(22)*1000)/expressions(89))*species(22)*species(9);
    ratelaws(81) = ((expressions(22)*1000)/expressions(89))*species(22)*species(14);
    ratelaws(82) = (expressions(71)/expressions(7))*species(9)*species(8);
    ratelaws(83) = (expressions(71)/expressions(8))*species(14)*species(13);
    ratelaws(84) = ((expressions(20)*1000)/expressions(89))*species(22)*species(6);
    ratelaws(85) = ((expressions(20)*1000)/expressions(89))*species(22)*species(11);
    ratelaws(86) = ((expressions(26)*1000)/expressions(89))*species(22)*species(7);
    ratelaws(87) = ((expressions(26)*1000)/expressions(89))*species(22)*species(12);
    ratelaws(88) = (expressions(70)/expressions(7))*species(7)*species(6);
    ratelaws(89) = (expressions(70)/expressions(8))*species(12)*species(11);
    ratelaws(90) = (expressions(69)/expressions(7))*species(8)*species(6);
    ratelaws(91) = (expressions(69)/expressions(8))*species(13)*species(11);
    ratelaws(92) = expressions(40)*species(45);
    ratelaws(93) = expressions(40)*species(46);
    ratelaws(94) = expressions(40)*species(53);
    ratelaws(95) = expressions(40)*species(54);
    ratelaws(96) = expressions(40)*species(57);
    ratelaws(97) = expressions(40)*species(58);
    ratelaws(98) = expressions(40)*species(59);
    ratelaws(99) = expressions(40)*species(60);
    ratelaws(100) = expressions(40)*species(47);
    ratelaws(101) = expressions(40)*species(48);
    ratelaws(102) = expressions(40)*species(55);
    ratelaws(103) = expressions(40)*species(56);
    ratelaws(104) = expressions(40)*species(63);
    ratelaws(105) = expressions(40)*species(64);
    ratelaws(106) = expressions(40)*species(71);
    ratelaws(107) = expressions(40)*species(72);
    ratelaws(108) = expressions(40)*species(51);
    ratelaws(109) = expressions(40)*species(52);
    ratelaws(110) = expressions(40)*species(49);
    ratelaws(111) = expressions(40)*species(50);
    ratelaws(112) = expressions(40)*species(65);
    ratelaws(113) = expressions(40)*species(66);
    ratelaws(114) = expressions(40)*species(69);
    ratelaws(115) = expressions(40)*species(70);
    ratelaws(116) = expressions(40)*species(67);
    ratelaws(117) = expressions(40)*species(68);
    ratelaws(118) = expressions(40)*species(33);
    ratelaws(119) = expressions(40)*species(34);
    ratelaws(120) = expressions(40)*species(37);
    ratelaws(121) = expressions(40)*species(38);
    ratelaws(122) = expressions(45)*species(30);
    ratelaws(123) = expressions(45)*species(31);
    ratelaws(124) = expressions(45)*species(32);
    ratelaws(125) = expressions(64)*species(30);
    ratelaws(126) = expressions(64)*species(31);
    ratelaws(127) = expressions(64)*species(32);
    ratelaws(128) = ((expressions(25)*1000)/expressions(89))*species(30)*species(10);
    ratelaws(129) = ((expressions(25)*1000)/expressions(89))*species(30)*species(15);
    ratelaws(130) = ((expressions(25)*1000)/expressions(89))*species(31)*species(10);
    ratelaws(131) = ((expressions(25)*1000)/expressions(89))*species(31)*species(15);
    ratelaws(132) = ((expressions(25)*1000)/expressions(89))*species(32)*species(10);
    ratelaws(133) = ((expressions(25)*1000)/expressions(89))*species(32)*species(15);
    ratelaws(134) = expressions(46)*species(33);
    ratelaws(135) = expressions(46)*species(34);
    ratelaws(136) = ((expressions(23)*1000)/expressions(89))*species(23)*species(33);
    ratelaws(137) = ((expressions(23)*1000)/expressions(89))*species(23)*species(34);
    ratelaws(138) = ((expressions(23)*1000)/expressions(89))*species(24)*species(33);
    ratelaws(139) = ((expressions(23)*1000)/expressions(89))*species(24)*species(34);
    ratelaws(140) = expressions(48)*species(35);
    ratelaws(141) = expressions(68)*species(35);
    ratelaws(142) = expressions(49)*species(36);
    ratelaws(143) = expressions(66)*species(35);
    ratelaws(144) = expressions(65)*species(36);
    ratelaws(145) = expressions(35)*species(37);
    ratelaws(146) = expressions(35)*species(38);
    ratelaws(147) = expressions(43)*species(39);
    ratelaws(148) = expressions(43)*species(40);
    ratelaws(149) = expressions(43)*species(41);
    ratelaws(150) = expressions(42)*species(42);
    ratelaws(151) = expressions(42)*species(43);
    ratelaws(152) = expressions(42)*species(44);
    ratelaws(153) = ((expressions(72)*1000)/expressions(89))*species(42)*species(28);
    ratelaws(154) = ((expressions(72)*1000)/expressions(89))*species(42)*species(27);
    ratelaws(155) = expressions(58)*species(45);
    ratelaws(156) = expressions(58)*species(46);
    ratelaws(157) = ((expressions(32)*1000)/expressions(89))*species(23)*species(63);
    ratelaws(158) = ((expressions(32)*1000)/expressions(89))*species(23)*species(64);
    ratelaws(159) = ((expressions(32)*1000)/expressions(89))*species(23)*species(71);
    ratelaws(160) = ((expressions(32)*1000)/expressions(89))*species(23)*species(72);
    ratelaws(161) = expressions(59)*species(47);
    ratelaws(162) = expressions(59)*species(48);
    ratelaws(163) = expressions(56)*species(49);
    ratelaws(164) = expressions(56)*species(50);
    ratelaws(165) = expressions(57)*species(51);
    ratelaws(166) = expressions(57)*species(52);
    ratelaws(167) = expressions(51)*species(53);
    ratelaws(168) = expressions(51)*species(54);
    ratelaws(169) = expressions(53)*species(55);
    ratelaws(170) = expressions(53)*species(56);
    ratelaws(171) = (expressions(37)/expressions(7))*species(8)*species(49);
    ratelaws(172) = (expressions(38)/expressions(7))*species(17)*species(47);
    ratelaws(173) = (expressions(37)/expressions(7))*species(8)*species(51);
    ratelaws(174) = (expressions(38)/expressions(7))*species(18)*species(47);
    ratelaws(175) = (expressions(37)/expressions(8))*species(13)*species(50);
    ratelaws(176) = (expressions(38)/expressions(8))*species(48)*species(19);
    ratelaws(177) = (expressions(37)/expressions(8))*species(13)*species(52);
    ratelaws(178) = (expressions(38)/expressions(8))*species(48)*species(21);
    ratelaws(179) = expressions(63)*species(57);
    ratelaws(180) = expressions(63)*species(58);
    ratelaws(181) = (expressions(39)/expressions(7))*species(53)*species(18);
    ratelaws(182) = (expressions(39)/expressions(7))*species(53)*species(17);
    ratelaws(183) = expressions(63)*species(59);
    ratelaws(184) = expressions(63)*species(60);
    ratelaws(185) = (expressions(39)/expressions(8))*species(54)*species(19);
    ratelaws(186) = (expressions(39)/expressions(8))*species(54)*species(21);
    ratelaws(187) = ((expressions(33)*1000)/expressions(89))*species(24)*species(57);
    ratelaws(188) = ((expressions(33)*1000)/expressions(89))*species(24)*species(60);
    ratelaws(189) = ((expressions(33)*1000)/expressions(89))*species(24)*species(58);
    ratelaws(190) = ((expressions(33)*1000)/expressions(89))*species(24)*species(59);
    ratelaws(191) = ((expressions(22)*1000)/expressions(89))*species(22)*species(63);
    ratelaws(192) = ((expressions(22)*1000)/expressions(89))*species(22)*species(64);
    ratelaws(193) = expressions(44)*species(61);
    ratelaws(194) = expressions(44)*species(62);
    ratelaws(195) = (expressions(71)/expressions(7))*species(9)*species(47);
    ratelaws(196) = (expressions(71)/expressions(7))*species(9)*species(55);
    ratelaws(197) = (expressions(71)/expressions(7))*species(61)*species(8);
    ratelaws(198) = (expressions(71)/expressions(7))*species(61)*species(47);
    ratelaws(199) = (expressions(71)/expressions(7))*species(61)*species(55);
    ratelaws(200) = expressions(61)*species(63);
    ratelaws(201) = (expressions(71)/expressions(8))*species(14)*species(48);
    ratelaws(202) = (expressions(71)/expressions(8))*species(14)*species(56);
    ratelaws(203) = (expressions(71)/expressions(8))*species(62)*species(13);
    ratelaws(204) = (expressions(71)/expressions(8))*species(62)*species(48);
    ratelaws(205) = (expressions(71)/expressions(8))*species(62)*species(56);
    ratelaws(206) = expressions(61)*species(64);
    ratelaws(207) = ((expressions(20)*1000)/expressions(89))*species(22)*species(71);
    ratelaws(208) = ((expressions(20)*1000)/expressions(89))*species(22)*species(72);
    ratelaws(209) = expressions(41)*species(65);
    ratelaws(210) = expressions(41)*species(66);
    ratelaws(211) = expressions(47)*species(67);
    ratelaws(212) = expressions(47)*species(68);
    ratelaws(213) = ((expressions(20)*1000)/expressions(89))*species(22)*species(69);
    ratelaws(214) = ((expressions(20)*1000)/expressions(89))*species(22)*species(70);
    ratelaws(215) = ((expressions(26)*1000)/expressions(89))*species(22)*species(69);
    ratelaws(216) = ((expressions(26)*1000)/expressions(89))*species(22)*species(70);
    ratelaws(217) = (expressions(70)/expressions(7))*species(7)*species(65);
    ratelaws(218) = (expressions(70)/expressions(7))*species(7)*species(71);
    ratelaws(219) = expressions(60)*species(69);
    ratelaws(220) = (expressions(70)/expressions(7))*species(67)*species(6);
    ratelaws(221) = (expressions(70)/expressions(7))*species(67)*species(71);
    ratelaws(222) = (expressions(70)/expressions(8))*species(12)*species(66);
    ratelaws(223) = (expressions(70)/expressions(8))*species(12)*species(72);
    ratelaws(224) = expressions(60)*species(70);
    ratelaws(225) = (expressions(70)/expressions(8))*species(68)*species(11);
    ratelaws(226) = (expressions(70)/expressions(8))*species(68)*species(72);
    ratelaws(227) = (expressions(69)/expressions(7))*species(8)*species(65);
    ratelaws(228) = (expressions(69)/expressions(7))*species(8)*species(69);
    ratelaws(229) = (expressions(69)/expressions(7))*species(47)*species(6);
    ratelaws(230) = (expressions(69)/expressions(7))*species(47)*species(65);
    ratelaws(231) = (expressions(69)/expressions(7))*species(47)*species(69);
    ratelaws(232) = (expressions(69)/expressions(7))*species(55)*species(6);
    ratelaws(233) = (expressions(69)/expressions(7))*species(55)*species(65);
    ratelaws(234) = (expressions(69)/expressions(7))*species(55)*species(69);
    ratelaws(235) = expressions(62)*species(71);
    ratelaws(236) = (expressions(69)/expressions(8))*species(13)*species(66);
    ratelaws(237) = (expressions(69)/expressions(8))*species(13)*species(70);
    ratelaws(238) = (expressions(69)/expressions(8))*species(48)*species(11);
    ratelaws(239) = (expressions(69)/expressions(8))*species(48)*species(66);
    ratelaws(240) = (expressions(69)/expressions(8))*species(48)*species(70);
    ratelaws(241) = (expressions(69)/expressions(8))*species(56)*species(11);
    ratelaws(242) = (expressions(69)/expressions(8))*species(56)*species(66);
    ratelaws(243) = (expressions(69)/expressions(8))*species(56)*species(70);
    ratelaws(244) = expressions(62)*species(72);
    ratelaws(245) = expressions(40)*species(91);
    ratelaws(246) = expressions(40)*species(92);
    ratelaws(247) = expressions(40)*species(93);
    ratelaws(248) = expressions(40)*species(94);
    ratelaws(249) = expressions(40)*species(79);
    ratelaws(250) = expressions(40)*species(80);
    ratelaws(251) = expressions(40)*species(81);
    ratelaws(252) = expressions(40)*species(82);
    ratelaws(253) = expressions(40)*species(83);
    ratelaws(254) = expressions(40)*species(84);
    ratelaws(255) = expressions(40)*species(85);
    ratelaws(256) = expressions(40)*species(86);
    ratelaws(257) = expressions(40)*species(87);
    ratelaws(258) = expressions(40)*species(88);
    ratelaws(259) = expressions(40)*species(89);
    ratelaws(260) = expressions(40)*species(90);
    ratelaws(261) = expressions(40)*species(95);
    ratelaws(262) = expressions(40)*species(96);
    ratelaws(263) = expressions(40)*species(97);
    ratelaws(264) = expressions(40)*species(98);
    ratelaws(265) = expressions(40)*species(99);
    ratelaws(266) = expressions(40)*species(100);
    ratelaws(267) = expressions(40)*species(101);
    ratelaws(268) = expressions(40)*species(102);
    ratelaws(269) = expressions(40)*species(103);
    ratelaws(270) = expressions(40)*species(104);
    ratelaws(271) = expressions(40)*species(109);
    ratelaws(272) = expressions(40)*species(110);
    ratelaws(273) = expressions(40)*species(111);
    ratelaws(274) = expressions(40)*species(112);
    ratelaws(275) = expressions(40)*species(113);
    ratelaws(276) = expressions(40)*species(114);
    ratelaws(277) = expressions(40)*species(115);
    ratelaws(278) = expressions(40)*species(116);
    ratelaws(279) = expressions(40)*species(117);
    ratelaws(280) = expressions(40)*species(118);
    ratelaws(281) = expressions(40)*species(119);
    ratelaws(282) = expressions(40)*species(120);
    ratelaws(283) = expressions(40)*species(121);
    ratelaws(284) = expressions(40)*species(122);
    ratelaws(285) = expressions(40)*species(105);
    ratelaws(286) = expressions(40)*species(106);
    ratelaws(287) = expressions(40)*species(107);
    ratelaws(288) = expressions(40)*species(108);
    ratelaws(289) = expressions(40)*species(73);
    ratelaws(290) = expressions(40)*species(74);
    ratelaws(291) = expressions(40)*species(75);
    ratelaws(292) = expressions(40)*species(76);
    ratelaws(293) = expressions(40)*species(77);
    ratelaws(294) = expressions(40)*species(78);
    ratelaws(295) = expressions(46)*species(73);
    ratelaws(296) = expressions(46)*species(74);
    ratelaws(297) = expressions(46)*species(75);
    ratelaws(298) = expressions(46)*species(76);
    ratelaws(299) = expressions(46)*species(77);
    ratelaws(300) = expressions(46)*species(78);
    ratelaws(301) = expressions(45)*species(73);
    ratelaws(302) = expressions(45)*species(74);
    ratelaws(303) = expressions(45)*species(75);
    ratelaws(304) = expressions(45)*species(76);
    ratelaws(305) = ((expressions(32)*1000)/expressions(89))*species(23)*species(95);
    ratelaws(306) = ((expressions(32)*1000)/expressions(89))*species(23)*species(96);
    ratelaws(307) = ((expressions(32)*1000)/expressions(89))*species(23)*species(103);
    ratelaws(308) = ((expressions(32)*1000)/expressions(89))*species(23)*species(104);
    ratelaws(309) = ((expressions(32)*1000)/expressions(89))*species(23)*species(109);
    ratelaws(310) = ((expressions(32)*1000)/expressions(89))*species(23)*species(110);
    ratelaws(311) = ((expressions(32)*1000)/expressions(89))*species(23)*species(111);
    ratelaws(312) = ((expressions(32)*1000)/expressions(89))*species(23)*species(112);
    ratelaws(313) = expressions(59)*species(79);
    ratelaws(314) = expressions(59)*species(80);
    ratelaws(315) = expressions(59)*species(81);
    ratelaws(316) = expressions(59)*species(82);
    ratelaws(317) = expressions(59)*species(98);
    ratelaws(318) = expressions(59)*species(101);
    ratelaws(319) = expressions(59)*species(113);
    ratelaws(320) = expressions(59)*species(114);
    ratelaws(321) = expressions(59)*species(118);
    ratelaws(322) = expressions(59)*species(119);
    ratelaws(323) = expressions(50)*species(83);
    ratelaws(324) = expressions(52)*species(84);
    ratelaws(325) = expressions(54)*species(85);
    ratelaws(326) = expressions(55)*species(86);
    ratelaws(327) = expressions(50)*species(87);
    ratelaws(328) = expressions(52)*species(88);
    ratelaws(329) = expressions(54)*species(89);
    ratelaws(330) = expressions(55)*species(90);
    ratelaws(331) = expressions(63)*species(91);
    ratelaws(332) = expressions(63)*species(92);
    ratelaws(333) = expressions(63)*species(93);
    ratelaws(334) = expressions(63)*species(94);
    ratelaws(335) = expressions(51)*species(91);
    ratelaws(336) = expressions(51)*species(94);
    ratelaws(337) = expressions(51)*species(92);
    ratelaws(338) = expressions(51)*species(93);
    ratelaws(339) = ((expressions(22)*1000)/expressions(89))*species(22)*species(79);
    ratelaws(340) = ((expressions(22)*1000)/expressions(89))*species(22)*species(80);
    ratelaws(341) = ((expressions(22)*1000)/expressions(89))*species(22)*species(97);
    ratelaws(342) = ((expressions(22)*1000)/expressions(89))*species(22)*species(100);
    ratelaws(343) = expressions(44)*species(95);
    ratelaws(344) = expressions(44)*species(96);
    ratelaws(345) = expressions(44)*species(98);
    ratelaws(346) = expressions(44)*species(99);
    ratelaws(347) = expressions(44)*species(101);
    ratelaws(348) = expressions(44)*species(102);
    ratelaws(349) = expressions(61)*species(79);
    ratelaws(350) = expressions(61)*species(95);
    ratelaws(351) = expressions(61)*species(97);
    ratelaws(352) = expressions(61)*species(98);
    ratelaws(353) = expressions(61)*species(99);
    ratelaws(354) = expressions(61)*species(80);
    ratelaws(355) = expressions(61)*species(96);
    ratelaws(356) = expressions(61)*species(100);
    ratelaws(357) = expressions(61)*species(101);
    ratelaws(358) = expressions(61)*species(102);
    ratelaws(359) = ((expressions(20)*1000)/expressions(89))*species(22)*species(81);
    ratelaws(360) = ((expressions(20)*1000)/expressions(89))*species(22)*species(82);
    ratelaws(361) = ((expressions(20)*1000)/expressions(89))*species(22)*species(115);
    ratelaws(362) = ((expressions(20)*1000)/expressions(89))*species(22)*species(120);
    ratelaws(363) = expressions(41)*species(103);
    ratelaws(364) = expressions(41)*species(104);
    ratelaws(365) = expressions(41)*species(113);
    ratelaws(366) = expressions(41)*species(116);
    ratelaws(367) = expressions(41)*species(118);
    ratelaws(368) = expressions(41)*species(121);
    ratelaws(369) = ((expressions(20)*1000)/expressions(89))*species(22)*species(109);
    ratelaws(370) = ((expressions(20)*1000)/expressions(89))*species(22)*species(111);
    ratelaws(371) = ((expressions(20)*1000)/expressions(89))*species(22)*species(114);
    ratelaws(372) = ((expressions(20)*1000)/expressions(89))*species(22)*species(117);
    ratelaws(373) = ((expressions(20)*1000)/expressions(89))*species(22)*species(119);
    ratelaws(374) = ((expressions(20)*1000)/expressions(89))*species(22)*species(122);
    ratelaws(375) = expressions(41)*species(105);
    ratelaws(376) = expressions(41)*species(106);
    ratelaws(377) = ((expressions(26)*1000)/expressions(89))*species(22)*species(109);
    ratelaws(378) = ((expressions(26)*1000)/expressions(89))*species(22)*species(111);
    ratelaws(379) = ((expressions(26)*1000)/expressions(89))*species(22)*species(114);
    ratelaws(380) = ((expressions(26)*1000)/expressions(89))*species(22)*species(117);
    ratelaws(381) = ((expressions(26)*1000)/expressions(89))*species(22)*species(119);
    ratelaws(382) = ((expressions(26)*1000)/expressions(89))*species(22)*species(122);
    ratelaws(383) = expressions(47)*species(107);
    ratelaws(384) = expressions(47)*species(108);
    ratelaws(385) = expressions(47)*species(110);
    ratelaws(386) = expressions(47)*species(112);
    ratelaws(387) = (expressions(70)/expressions(7))*species(7)*species(81);
    ratelaws(388) = (expressions(70)/expressions(7))*species(7)*species(103);
    ratelaws(389) = (expressions(70)/expressions(7))*species(7)*species(113);
    ratelaws(390) = (expressions(70)/expressions(7))*species(7)*species(115);
    ratelaws(391) = (expressions(70)/expressions(7))*species(7)*species(116);
    ratelaws(392) = expressions(60)*species(105);
    ratelaws(393) = expressions(60)*species(109);
    ratelaws(394) = expressions(60)*species(114);
    ratelaws(395) = expressions(60)*species(117);
    ratelaws(396) = (expressions(70)/expressions(7))*species(67)*species(81);
    ratelaws(397) = (expressions(70)/expressions(7))*species(67)*species(115);
    ratelaws(398) = expressions(60)*species(107);
    ratelaws(399) = expressions(60)*species(110);
    ratelaws(400) = (expressions(70)/expressions(8))*species(12)*species(82);
    ratelaws(401) = (expressions(70)/expressions(8))*species(12)*species(104);
    ratelaws(402) = (expressions(70)/expressions(8))*species(12)*species(118);
    ratelaws(403) = (expressions(70)/expressions(8))*species(12)*species(120);
    ratelaws(404) = (expressions(70)/expressions(8))*species(12)*species(121);
    ratelaws(405) = expressions(60)*species(106);
    ratelaws(406) = expressions(60)*species(111);
    ratelaws(407) = expressions(60)*species(119);
    ratelaws(408) = expressions(60)*species(122);
    ratelaws(409) = (expressions(70)/expressions(8))*species(68)*species(82);
    ratelaws(410) = (expressions(70)/expressions(8))*species(68)*species(120);
    ratelaws(411) = expressions(60)*species(108);
    ratelaws(412) = expressions(60)*species(112);
    ratelaws(413) = (expressions(69)/expressions(7))*species(8)*species(105);
    ratelaws(414) = (expressions(69)/expressions(7))*species(8)*species(107);
    ratelaws(415) = (expressions(69)/expressions(7))*species(47)*species(105);
    ratelaws(416) = (expressions(69)/expressions(7))*species(47)*species(107);
    ratelaws(417) = (expressions(69)/expressions(7))*species(55)*species(105);
    ratelaws(418) = (expressions(69)/expressions(7))*species(55)*species(107);
    ratelaws(419) = expressions(62)*species(81);
    ratelaws(420) = expressions(62)*species(103);
    ratelaws(421) = expressions(62)*species(109);
    ratelaws(422) = expressions(62)*species(110);
    ratelaws(423) = expressions(62)*species(113);
    ratelaws(424) = expressions(62)*species(114);
    ratelaws(425) = expressions(62)*species(115);
    ratelaws(426) = expressions(62)*species(116);
    ratelaws(427) = expressions(62)*species(117);
    ratelaws(428) = (expressions(69)/expressions(8))*species(13)*species(106);
    ratelaws(429) = (expressions(69)/expressions(8))*species(13)*species(108);
    ratelaws(430) = (expressions(69)/expressions(8))*species(48)*species(106);
    ratelaws(431) = (expressions(69)/expressions(8))*species(48)*species(108);
    ratelaws(432) = (expressions(69)/expressions(8))*species(56)*species(106);
    ratelaws(433) = (expressions(69)/expressions(8))*species(56)*species(108);
    ratelaws(434) = expressions(62)*species(82);
    ratelaws(435) = expressions(62)*species(104);
    ratelaws(436) = expressions(62)*species(111);
    ratelaws(437) = expressions(62)*species(112);
    ratelaws(438) = expressions(62)*species(118);
    ratelaws(439) = expressions(62)*species(119);
    ratelaws(440) = expressions(62)*species(120);
    ratelaws(441) = expressions(62)*species(121);
    ratelaws(442) = expressions(62)*species(122);
    ratelaws(443) = expressions(40)*species(123);
    ratelaws(444) = expressions(40)*species(124);
    ratelaws(445) = expressions(40)*species(125);
    ratelaws(446) = expressions(40)*species(126);
    ratelaws(447) = expressions(40)*species(127);
    ratelaws(448) = expressions(40)*species(128);
    ratelaws(449) = expressions(40)*species(129);
    ratelaws(450) = expressions(40)*species(130);
    ratelaws(451) = expressions(40)*species(131);
    ratelaws(452) = expressions(40)*species(132);
    ratelaws(453) = ((expressions(32)*1000)/expressions(89))*species(23)*species(125);
    ratelaws(454) = ((expressions(32)*1000)/expressions(89))*species(23)*species(126);
    ratelaws(455) = expressions(59)*species(123);
    ratelaws(456) = expressions(59)*species(124);
    ratelaws(457) = expressions(59)*species(127);
    ratelaws(458) = expressions(59)*species(129);
    ratelaws(459) = expressions(41)*species(125);
    ratelaws(460) = expressions(41)*species(126);
    ratelaws(461) = expressions(41)*species(127);
    ratelaws(462) = expressions(41)*species(128);
    ratelaws(463) = expressions(41)*species(129);
    ratelaws(464) = expressions(41)*species(130);
    ratelaws(465) = expressions(47)*species(123);
    ratelaws(466) = expressions(47)*species(124);
    ratelaws(467) = expressions(47)*species(131);
    ratelaws(468) = expressions(47)*species(132);
    ratelaws(469) = expressions(60)*species(125);
    ratelaws(470) = expressions(60)*species(127);
    ratelaws(471) = expressions(60)*species(128);
    ratelaws(472) = expressions(60)*species(123);
    ratelaws(473) = expressions(60)*species(131);
    ratelaws(474) = expressions(60)*species(126);
    ratelaws(475) = expressions(60)*species(129);
    ratelaws(476) = expressions(60)*species(130);
    ratelaws(477) = expressions(60)*species(124);
    ratelaws(478) = expressions(60)*species(132);
    ratelaws(479) = expressions(62)*species(123);
    ratelaws(480) = expressions(62)*species(125);
    ratelaws(481) = expressions(62)*species(127);
    ratelaws(482) = expressions(62)*species(128);
    ratelaws(483) = expressions(62)*species(131);
    ratelaws(484) = expressions(62)*species(124);
    ratelaws(485) = expressions(62)*species(126);
    ratelaws(486) = expressions(62)*species(129);
    ratelaws(487) = expressions(62)*species(130);
    ratelaws(488) = expressions(62)*species(132);

end

% Calculate species derivatives
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros(132,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calc_ratelaws( species, expressions, observables );
                        
    % calculate derivatives
    Dspecies(1) = ratelaws(17) +ratelaws(18) +ratelaws(19) +ratelaws(20) +ratelaws(21) +ratelaws(22) +ratelaws(23) +ratelaws(24) +ratelaws(25) +ratelaws(26) +ratelaws(27) +ratelaws(28) +ratelaws(29) +ratelaws(30) +ratelaws(31) +ratelaws(32) +ratelaws(34) +ratelaws(36) +ratelaws(38) +ratelaws(40) +ratelaws(41) +ratelaws(51) +ratelaws(52) +ratelaws(53) +ratelaws(92) +ratelaws(93) +ratelaws(94) +ratelaws(95) +ratelaws(96) +ratelaws(97) +ratelaws(98) +ratelaws(99) +ratelaws(100) +ratelaws(101) +ratelaws(102) +ratelaws(103) +ratelaws(104) +ratelaws(105) +ratelaws(106) +ratelaws(107) +ratelaws(108) +ratelaws(109) +ratelaws(110) +ratelaws(111) +ratelaws(112) +ratelaws(113) +ratelaws(114) +ratelaws(115) +ratelaws(116) +ratelaws(117) +ratelaws(118) +ratelaws(119) +ratelaws(120) +ratelaws(121) +ratelaws(125) +ratelaws(126) +ratelaws(127) +ratelaws(143) +ratelaws(144) +ratelaws(245) +ratelaws(246) +ratelaws(247) +ratelaws(248) +ratelaws(249) +ratelaws(250) +ratelaws(251) +ratelaws(252) +ratelaws(253) +ratelaws(254) +ratelaws(255) +ratelaws(256) +ratelaws(257) +ratelaws(258) +ratelaws(259) +ratelaws(260) +ratelaws(261) +ratelaws(262) +ratelaws(263) +ratelaws(264) +ratelaws(265) +ratelaws(266) +ratelaws(267) +ratelaws(268) +ratelaws(269) +ratelaws(270) +ratelaws(271) +ratelaws(272) +ratelaws(273) +ratelaws(274) +ratelaws(275) +ratelaws(276) +ratelaws(277) +ratelaws(278) +ratelaws(279) +ratelaws(280) +ratelaws(281) +ratelaws(282) +ratelaws(283) +ratelaws(284) +ratelaws(285) +ratelaws(286) +ratelaws(287) +ratelaws(288) +ratelaws(289) +ratelaws(290) +ratelaws(291) +ratelaws(292) +ratelaws(293) +ratelaws(294) +ratelaws(443) +ratelaws(444) +ratelaws(445) +ratelaws(446) +ratelaws(447) +ratelaws(448) +ratelaws(449) +ratelaws(450) +ratelaws(451) +ratelaws(452);
    Dspecies(2) = 0.0;
    Dspecies(3) = -ratelaws(58) -ratelaws(61) +ratelaws(147) +ratelaws(150) +ratelaws(153) +ratelaws(154);
    Dspecies(4) = -ratelaws(59) -ratelaws(62) +ratelaws(148) +ratelaws(151);
    Dspecies(5) = -ratelaws(60) -ratelaws(63) +ratelaws(149) +ratelaws(152);
    Dspecies(6) = ratelaws(1) -ratelaws(25) -ratelaws(84) -ratelaws(88) -ratelaws(90) +ratelaws(209) +ratelaws(219) -ratelaws(220) -ratelaws(229) -ratelaws(232) +ratelaws(235) +ratelaws(398) +ratelaws(419) +ratelaws(425);
    Dspecies(7) = ratelaws(3) -ratelaws(27) -ratelaws(86) -ratelaws(88) +ratelaws(211) -ratelaws(217) -ratelaws(218) +ratelaws(219) -ratelaws(387) -ratelaws(388) -ratelaws(389) -ratelaws(390) -ratelaws(391) +ratelaws(392) +ratelaws(393) +ratelaws(394) +ratelaws(395) +ratelaws(469) +ratelaws(470) +ratelaws(471);
    Dspecies(8) = ratelaws(10) -ratelaws(19) -ratelaws(66) -ratelaws(74) -ratelaws(82) -ratelaws(90) +ratelaws(161) +ratelaws(169) -ratelaws(171) -ratelaws(173) -ratelaws(197) +ratelaws(200) -ratelaws(227) -ratelaws(228) +ratelaws(235) +ratelaws(323) +ratelaws(325) +ratelaws(350) -ratelaws(413) -ratelaws(414) +ratelaws(420) +ratelaws(421) +ratelaws(422) +ratelaws(480);
    Dspecies(9) = ratelaws(5) -ratelaws(29) -ratelaws(80) -ratelaws(82) +ratelaws(193) -ratelaws(195) -ratelaws(196) +ratelaws(200) +ratelaws(349) +ratelaws(351);
    Dspecies(10) = ratelaws(7) -ratelaws(31) -ratelaws(45) -ratelaws(56) -ratelaws(128) -ratelaws(130) -ratelaws(132) +ratelaws(134) +ratelaws(145) +ratelaws(295) +ratelaws(297) +ratelaws(299);
    Dspecies(11) = ratelaws(2) -ratelaws(26) -ratelaws(85) -ratelaws(89) -ratelaws(91) +ratelaws(210) +ratelaws(224) -ratelaws(225) -ratelaws(238) -ratelaws(241) +ratelaws(244) +ratelaws(411) +ratelaws(434) +ratelaws(440);
    Dspecies(12) = ratelaws(4) -ratelaws(28) -ratelaws(87) -ratelaws(89) +ratelaws(212) -ratelaws(222) -ratelaws(223) +ratelaws(224) -ratelaws(400) -ratelaws(401) -ratelaws(402) -ratelaws(403) -ratelaws(404) +ratelaws(405) +ratelaws(406) +ratelaws(407) +ratelaws(408) +ratelaws(474) +ratelaws(475) +ratelaws(476);
    Dspecies(13) = ratelaws(14) -ratelaws(20) -ratelaws(67) -ratelaws(75) -ratelaws(83) -ratelaws(91) +ratelaws(162) +ratelaws(170) -ratelaws(175) -ratelaws(177) -ratelaws(203) +ratelaws(206) -ratelaws(236) -ratelaws(237) +ratelaws(244) +ratelaws(327) +ratelaws(329) +ratelaws(355) -ratelaws(428) -ratelaws(429) +ratelaws(435) +ratelaws(436) +ratelaws(437) +ratelaws(485);
    Dspecies(14) = ratelaws(6) -ratelaws(30) -ratelaws(81) -ratelaws(83) +ratelaws(194) -ratelaws(201) -ratelaws(202) +ratelaws(206) +ratelaws(354) +ratelaws(356);
    Dspecies(15) = ratelaws(8) -ratelaws(32) -ratelaws(46) -ratelaws(57) -ratelaws(129) -ratelaws(131) -ratelaws(133) +ratelaws(135) +ratelaws(146) +ratelaws(296) +ratelaws(298) +ratelaws(300);
    Dspecies(16) = ratelaws(9) -ratelaws(17) -ratelaws(64) -ratelaws(72) -ratelaws(76) -ratelaws(77) +ratelaws(155) +ratelaws(167) +ratelaws(179) +ratelaws(180);
    Dspecies(17) = ratelaws(11) -ratelaws(23) -ratelaws(68) -ratelaws(77) +ratelaws(163) -ratelaws(172) +ratelaws(180) -ratelaws(182) +ratelaws(324) +ratelaws(332);
    Dspecies(18) = ratelaws(12) -ratelaws(21) -ratelaws(70) -ratelaws(76) +ratelaws(165) -ratelaws(174) +ratelaws(179) -ratelaws(181) +ratelaws(326) +ratelaws(331);
    Dspecies(19) = ratelaws(15) -ratelaws(24) -ratelaws(69) -ratelaws(78) +ratelaws(164) -ratelaws(176) +ratelaws(183) -ratelaws(185) +ratelaws(328) +ratelaws(333);
    Dspecies(20) = ratelaws(13) -ratelaws(18) -ratelaws(65) -ratelaws(73) -ratelaws(78) -ratelaws(79) +ratelaws(156) +ratelaws(168) +ratelaws(183) +ratelaws(184);
    Dspecies(21) = ratelaws(16) -ratelaws(22) -ratelaws(71) -ratelaws(79) +ratelaws(166) -ratelaws(178) +ratelaws(184) -ratelaws(186) +ratelaws(330) +ratelaws(334);
    Dspecies(22) = ratelaws(33) -ratelaws(34) -ratelaws(35) -ratelaws(42) -ratelaws(43) -ratelaws(44) -ratelaws(45) -ratelaws(46) -ratelaws(48) -ratelaws(58) -ratelaws(59) -ratelaws(60) -ratelaws(80) -ratelaws(81) -ratelaws(84) -ratelaws(85) -ratelaws(86) -ratelaws(87) +ratelaws(122) +ratelaws(123) +ratelaws(124) +ratelaws(134) +ratelaws(135) +ratelaws(142) +ratelaws(147) +ratelaws(148) +ratelaws(149) -ratelaws(191) -ratelaws(192) +ratelaws(193) +ratelaws(194) -ratelaws(207) -ratelaws(208) +ratelaws(209) +ratelaws(210) +ratelaws(211) +ratelaws(212) -ratelaws(213) -ratelaws(214) -ratelaws(215) -ratelaws(216) -ratelaws(339) -ratelaws(340) -ratelaws(341) -ratelaws(342) +ratelaws(343) +ratelaws(344) +ratelaws(345) +ratelaws(346) +ratelaws(347) +ratelaws(348) -ratelaws(359) -ratelaws(360) -ratelaws(361) -ratelaws(362) +ratelaws(363) +ratelaws(364) +ratelaws(365) +ratelaws(366) +ratelaws(367) +ratelaws(368) -ratelaws(369) -ratelaws(370) -ratelaws(371) -ratelaws(372) -ratelaws(373) -ratelaws(374) +ratelaws(375) +ratelaws(376) -ratelaws(377) -ratelaws(378) -ratelaws(379) -ratelaws(380) -ratelaws(381) -ratelaws(382) +ratelaws(383) +ratelaws(384) +ratelaws(385) +ratelaws(386) +ratelaws(459) +ratelaws(460) +ratelaws(461) +ratelaws(462) +ratelaws(463) +ratelaws(464) +ratelaws(465) +ratelaws(466) +ratelaws(467) +ratelaws(468);
    Dspecies(23) = ratelaws(37) -ratelaws(38) -ratelaws(42) -ratelaws(54) -ratelaws(55) -ratelaws(61) -ratelaws(62) -ratelaws(63) -ratelaws(64) -ratelaws(65) -ratelaws(66) -ratelaws(67) -ratelaws(68) -ratelaws(69) -ratelaws(70) -ratelaws(71) +ratelaws(122) -ratelaws(136) -ratelaws(137) +ratelaws(150) +ratelaws(151) +ratelaws(152) +ratelaws(155) +ratelaws(156) -ratelaws(157) -ratelaws(158) -ratelaws(159) -ratelaws(160) +ratelaws(161) +ratelaws(162) +ratelaws(163) +ratelaws(164) +ratelaws(165) +ratelaws(166) +ratelaws(301) +ratelaws(302) -ratelaws(305) -ratelaws(306) -ratelaws(307) -ratelaws(308) -ratelaws(309) -ratelaws(310) -ratelaws(311) -ratelaws(312) +ratelaws(313) +ratelaws(314) +ratelaws(315) +ratelaws(316) +ratelaws(317) +ratelaws(318) +ratelaws(319) +ratelaws(320) +ratelaws(321) +ratelaws(322) -ratelaws(453) -ratelaws(454) +ratelaws(455) +ratelaws(456) +ratelaws(457) +ratelaws(458);
    Dspecies(24) = ratelaws(39) -ratelaws(40) -ratelaws(43) -ratelaws(72) -ratelaws(73) -ratelaws(74) -ratelaws(75) +ratelaws(123) -ratelaws(138) -ratelaws(139) +ratelaws(167) +ratelaws(168) +ratelaws(169) +ratelaws(170) -ratelaws(187) -ratelaws(188) -ratelaws(189) -ratelaws(190) +ratelaws(303) +ratelaws(304) +ratelaws(335) +ratelaws(336) +ratelaws(337) +ratelaws(338);
    Dspecies(25) = -ratelaws(41) -ratelaws(44) +ratelaws(54) +ratelaws(55) +ratelaws(124) +ratelaws(153) +ratelaws(154);
    Dspecies(26) = -ratelaws(47) +ratelaws(49) -ratelaws(53) +ratelaws(140);
    Dspecies(27) = -ratelaws(52) -ratelaws(56) -ratelaws(57) +ratelaws(141) +ratelaws(145) +ratelaws(146);
    Dspecies(28) = -ratelaws(47) -ratelaws(48) +ratelaws(50) -ratelaws(51) +ratelaws(140) +ratelaws(141) +ratelaws(142);
    Dspecies(29) = ratelaws(35) -ratelaws(36);
    Dspecies(30) = ratelaws(42) -ratelaws(122) -ratelaws(125) -ratelaws(128) -ratelaws(129) +ratelaws(295) +ratelaws(296);
    Dspecies(31) = ratelaws(43) -ratelaws(123) -ratelaws(126) -ratelaws(130) -ratelaws(131) +ratelaws(297) +ratelaws(298);
    Dspecies(32) = ratelaws(44) -ratelaws(124) -ratelaws(127) -ratelaws(132) -ratelaws(133) +ratelaws(299) +ratelaws(300);
    Dspecies(33) = ratelaws(45) -ratelaws(118) -ratelaws(134) -ratelaws(136) -ratelaws(138) +ratelaws(301) +ratelaws(303);
    Dspecies(34) = ratelaws(46) -ratelaws(119) -ratelaws(135) -ratelaws(137) -ratelaws(139) +ratelaws(302) +ratelaws(304);
    Dspecies(35) = ratelaws(47) -ratelaws(140) -ratelaws(141) -ratelaws(143);
    Dspecies(36) = ratelaws(48) -ratelaws(142) -ratelaws(144);
    Dspecies(37) = ratelaws(56) -ratelaws(120) -ratelaws(145);
    Dspecies(38) = ratelaws(57) -ratelaws(121) -ratelaws(146);
    Dspecies(39) = ratelaws(58) -ratelaws(147);
    Dspecies(40) = ratelaws(59) -ratelaws(148);
    Dspecies(41) = ratelaws(60) -ratelaws(149);
    Dspecies(42) = ratelaws(61) -ratelaws(150) -ratelaws(153) -ratelaws(154);
    Dspecies(43) = ratelaws(62) -ratelaws(151);
    Dspecies(44) = ratelaws(63) -ratelaws(152);
    Dspecies(45) = ratelaws(64) -ratelaws(92) -ratelaws(155);
    Dspecies(46) = ratelaws(65) -ratelaws(93) -ratelaws(156);
    Dspecies(47) = ratelaws(66) -ratelaws(100) -ratelaws(161) -ratelaws(172) -ratelaws(174) -ratelaws(195) -ratelaws(198) -ratelaws(229) -ratelaws(230) -ratelaws(231) +ratelaws(324) +ratelaws(326) +ratelaws(349) +ratelaws(352) -ratelaws(415) -ratelaws(416) +ratelaws(419) +ratelaws(423) +ratelaws(424) +ratelaws(479) +ratelaws(481);
    Dspecies(48) = ratelaws(67) -ratelaws(101) -ratelaws(162) -ratelaws(176) -ratelaws(178) -ratelaws(201) -ratelaws(204) -ratelaws(238) -ratelaws(239) -ratelaws(240) +ratelaws(328) +ratelaws(330) +ratelaws(354) +ratelaws(357) -ratelaws(430) -ratelaws(431) +ratelaws(434) +ratelaws(438) +ratelaws(439) +ratelaws(484) +ratelaws(486);
    Dspecies(49) = ratelaws(68) -ratelaws(110) -ratelaws(163) -ratelaws(171) +ratelaws(323);
    Dspecies(50) = ratelaws(69) -ratelaws(111) -ratelaws(164) -ratelaws(175) +ratelaws(327);
    Dspecies(51) = ratelaws(70) -ratelaws(108) -ratelaws(165) -ratelaws(173) +ratelaws(325);
    Dspecies(52) = ratelaws(71) -ratelaws(109) -ratelaws(166) -ratelaws(177) +ratelaws(329);
    Dspecies(53) = ratelaws(72) -ratelaws(94) -ratelaws(167) -ratelaws(181) -ratelaws(182) +ratelaws(331) +ratelaws(332);
    Dspecies(54) = ratelaws(73) -ratelaws(95) -ratelaws(168) -ratelaws(185) -ratelaws(186) +ratelaws(333) +ratelaws(334);
    Dspecies(55) = ratelaws(74) -ratelaws(102) -ratelaws(169) -ratelaws(196) -ratelaws(199) -ratelaws(232) -ratelaws(233) -ratelaws(234) +ratelaws(351) +ratelaws(353) -ratelaws(417) -ratelaws(418) +ratelaws(425) +ratelaws(426) +ratelaws(427) +ratelaws(482) +ratelaws(483);
    Dspecies(56) = ratelaws(75) -ratelaws(103) -ratelaws(170) -ratelaws(202) -ratelaws(205) -ratelaws(241) -ratelaws(242) -ratelaws(243) +ratelaws(356) +ratelaws(358) -ratelaws(432) -ratelaws(433) +ratelaws(440) +ratelaws(441) +ratelaws(442) +ratelaws(487) +ratelaws(488);
    Dspecies(57) = ratelaws(76) -ratelaws(96) -ratelaws(179) -ratelaws(187) +ratelaws(335);
    Dspecies(58) = ratelaws(77) -ratelaws(97) -ratelaws(180) -ratelaws(189) +ratelaws(337);
    Dspecies(59) = ratelaws(78) -ratelaws(98) -ratelaws(183) -ratelaws(190) +ratelaws(338);
    Dspecies(60) = ratelaws(79) -ratelaws(99) -ratelaws(184) -ratelaws(188) +ratelaws(336);
    Dspecies(61) = ratelaws(80) -ratelaws(193) -ratelaws(197) -ratelaws(198) -ratelaws(199) +ratelaws(350) +ratelaws(352) +ratelaws(353);
    Dspecies(62) = ratelaws(81) -ratelaws(194) -ratelaws(203) -ratelaws(204) -ratelaws(205) +ratelaws(355) +ratelaws(357) +ratelaws(358);
    Dspecies(63) = ratelaws(82) -ratelaws(104) -ratelaws(157) -ratelaws(191) -ratelaws(200) +ratelaws(313) +ratelaws(343);
    Dspecies(64) = ratelaws(83) -ratelaws(105) -ratelaws(158) -ratelaws(192) -ratelaws(206) +ratelaws(314) +ratelaws(344);
    Dspecies(65) = ratelaws(84) -ratelaws(112) -ratelaws(209) -ratelaws(217) -ratelaws(227) -ratelaws(230) -ratelaws(233) +ratelaws(392) +ratelaws(420) +ratelaws(423) +ratelaws(426);
    Dspecies(66) = ratelaws(85) -ratelaws(113) -ratelaws(210) -ratelaws(222) -ratelaws(236) -ratelaws(239) -ratelaws(242) +ratelaws(405) +ratelaws(435) +ratelaws(438) +ratelaws(441);
    Dspecies(67) = ratelaws(86) -ratelaws(116) -ratelaws(211) -ratelaws(220) -ratelaws(221) -ratelaws(396) -ratelaws(397) +ratelaws(398) +ratelaws(399) +ratelaws(472) +ratelaws(473);
    Dspecies(68) = ratelaws(87) -ratelaws(117) -ratelaws(212) -ratelaws(225) -ratelaws(226) -ratelaws(409) -ratelaws(410) +ratelaws(411) +ratelaws(412) +ratelaws(477) +ratelaws(478);
    Dspecies(69) = ratelaws(88) -ratelaws(114) -ratelaws(213) -ratelaws(215) -ratelaws(219) -ratelaws(228) -ratelaws(231) -ratelaws(234) +ratelaws(375) +ratelaws(383) +ratelaws(421) +ratelaws(424) +ratelaws(427);
    Dspecies(70) = ratelaws(89) -ratelaws(115) -ratelaws(214) -ratelaws(216) -ratelaws(224) -ratelaws(237) -ratelaws(240) -ratelaws(243) +ratelaws(376) +ratelaws(384) +ratelaws(436) +ratelaws(439) +ratelaws(442);
    Dspecies(71) = ratelaws(90) -ratelaws(106) -ratelaws(159) -ratelaws(207) -ratelaws(218) -ratelaws(221) -ratelaws(235) +ratelaws(315) +ratelaws(363) +ratelaws(393) +ratelaws(399);
    Dspecies(72) = ratelaws(91) -ratelaws(107) -ratelaws(160) -ratelaws(208) -ratelaws(223) -ratelaws(226) -ratelaws(244) +ratelaws(316) +ratelaws(364) +ratelaws(406) +ratelaws(412);
    Dspecies(73) = ratelaws(128) +ratelaws(136) -ratelaws(289) -ratelaws(295) -ratelaws(301);
    Dspecies(74) = ratelaws(129) +ratelaws(137) -ratelaws(290) -ratelaws(296) -ratelaws(302);
    Dspecies(75) = ratelaws(130) +ratelaws(138) -ratelaws(291) -ratelaws(297) -ratelaws(303);
    Dspecies(76) = ratelaws(131) +ratelaws(139) -ratelaws(292) -ratelaws(298) -ratelaws(304);
    Dspecies(77) = ratelaws(132) -ratelaws(293) -ratelaws(299);
    Dspecies(78) = ratelaws(133) -ratelaws(294) -ratelaws(300);
    Dspecies(79) = ratelaws(157) +ratelaws(195) -ratelaws(249) -ratelaws(313) -ratelaws(339) +ratelaws(345) -ratelaws(349);
    Dspecies(80) = ratelaws(158) +ratelaws(201) -ratelaws(250) -ratelaws(314) -ratelaws(340) +ratelaws(347) -ratelaws(354);
    Dspecies(81) = ratelaws(159) +ratelaws(229) -ratelaws(251) -ratelaws(315) -ratelaws(359) +ratelaws(365) -ratelaws(387) +ratelaws(394) -ratelaws(396) -ratelaws(419) +ratelaws(472);
    Dspecies(82) = ratelaws(160) +ratelaws(238) -ratelaws(252) -ratelaws(316) -ratelaws(360) +ratelaws(367) -ratelaws(400) +ratelaws(407) -ratelaws(409) -ratelaws(434) +ratelaws(477);
    Dspecies(83) = ratelaws(171) -ratelaws(253) -ratelaws(323);
    Dspecies(84) = ratelaws(172) -ratelaws(254) -ratelaws(324);
    Dspecies(85) = ratelaws(173) -ratelaws(255) -ratelaws(325);
    Dspecies(86) = ratelaws(174) -ratelaws(256) -ratelaws(326);
    Dspecies(87) = ratelaws(175) -ratelaws(257) -ratelaws(327);
    Dspecies(88) = ratelaws(176) -ratelaws(258) -ratelaws(328);
    Dspecies(89) = ratelaws(177) -ratelaws(259) -ratelaws(329);
    Dspecies(90) = ratelaws(178) -ratelaws(260) -ratelaws(330);
    Dspecies(91) = ratelaws(181) +ratelaws(187) -ratelaws(245) -ratelaws(331) -ratelaws(335);
    Dspecies(92) = ratelaws(182) +ratelaws(189) -ratelaws(246) -ratelaws(332) -ratelaws(337);
    Dspecies(93) = ratelaws(185) +ratelaws(190) -ratelaws(247) -ratelaws(333) -ratelaws(338);
    Dspecies(94) = ratelaws(186) +ratelaws(188) -ratelaws(248) -ratelaws(334) -ratelaws(336);
    Dspecies(95) = ratelaws(191) +ratelaws(197) -ratelaws(261) -ratelaws(305) +ratelaws(317) -ratelaws(343) -ratelaws(350);
    Dspecies(96) = ratelaws(192) +ratelaws(203) -ratelaws(262) -ratelaws(306) +ratelaws(318) -ratelaws(344) -ratelaws(355);
    Dspecies(97) = ratelaws(196) -ratelaws(263) -ratelaws(341) +ratelaws(346) -ratelaws(351);
    Dspecies(98) = ratelaws(198) -ratelaws(264) +ratelaws(305) -ratelaws(317) +ratelaws(339) -ratelaws(345) -ratelaws(352);
    Dspecies(99) = ratelaws(199) -ratelaws(265) +ratelaws(341) -ratelaws(346) -ratelaws(353);
    Dspecies(100) = ratelaws(202) -ratelaws(266) -ratelaws(342) +ratelaws(348) -ratelaws(356);
    Dspecies(101) = ratelaws(204) -ratelaws(267) +ratelaws(306) -ratelaws(318) +ratelaws(340) -ratelaws(347) -ratelaws(357);
    Dspecies(102) = ratelaws(205) -ratelaws(268) +ratelaws(342) -ratelaws(348) -ratelaws(358);
    Dspecies(103) = ratelaws(207) +ratelaws(227) -ratelaws(269) -ratelaws(307) +ratelaws(319) -ratelaws(363) -ratelaws(388) -ratelaws(420) +ratelaws(469);
    Dspecies(104) = ratelaws(208) +ratelaws(236) -ratelaws(270) -ratelaws(308) +ratelaws(321) -ratelaws(364) -ratelaws(401) -ratelaws(435) +ratelaws(474);
    Dspecies(105) = ratelaws(213) +ratelaws(217) -ratelaws(285) -ratelaws(375) -ratelaws(392) -ratelaws(413) -ratelaws(415) -ratelaws(417) +ratelaws(480) +ratelaws(481) +ratelaws(482);
    Dspecies(106) = ratelaws(214) +ratelaws(222) -ratelaws(286) -ratelaws(376) -ratelaws(405) -ratelaws(428) -ratelaws(430) -ratelaws(432) +ratelaws(485) +ratelaws(486) +ratelaws(487);
    Dspecies(107) = ratelaws(215) +ratelaws(220) -ratelaws(287) -ratelaws(383) -ratelaws(398) -ratelaws(414) -ratelaws(416) -ratelaws(418) +ratelaws(422) +ratelaws(479) +ratelaws(483);
    Dspecies(108) = ratelaws(216) +ratelaws(225) -ratelaws(288) -ratelaws(384) -ratelaws(411) -ratelaws(429) -ratelaws(431) -ratelaws(433) +ratelaws(437) +ratelaws(484) +ratelaws(488);
    Dspecies(109) = ratelaws(218) +ratelaws(228) -ratelaws(271) -ratelaws(309) +ratelaws(320) -ratelaws(369) -ratelaws(377) +ratelaws(385) -ratelaws(393) -ratelaws(421) +ratelaws(459);
    Dspecies(110) = ratelaws(221) -ratelaws(272) -ratelaws(310) +ratelaws(377) -ratelaws(385) -ratelaws(399) +ratelaws(414) -ratelaws(422) +ratelaws(455);
    Dspecies(111) = ratelaws(223) +ratelaws(237) -ratelaws(273) -ratelaws(311) +ratelaws(322) -ratelaws(370) -ratelaws(378) +ratelaws(386) -ratelaws(406) -ratelaws(436) +ratelaws(460);
    Dspecies(112) = ratelaws(226) -ratelaws(274) -ratelaws(312) +ratelaws(378) -ratelaws(386) -ratelaws(412) +ratelaws(429) -ratelaws(437) +ratelaws(456);
    Dspecies(113) = ratelaws(230) -ratelaws(275) +ratelaws(307) -ratelaws(319) +ratelaws(359) -ratelaws(365) -ratelaws(389) -ratelaws(423) +ratelaws(470);
    Dspecies(114) = ratelaws(231) -ratelaws(276) +ratelaws(309) -ratelaws(320) -ratelaws(371) -ratelaws(379) +ratelaws(387) -ratelaws(394) -ratelaws(424) +ratelaws(461) +ratelaws(465);
    Dspecies(115) = ratelaws(232) -ratelaws(277) -ratelaws(361) +ratelaws(366) -ratelaws(390) +ratelaws(395) -ratelaws(397) -ratelaws(425) +ratelaws(473);
    Dspecies(116) = ratelaws(233) -ratelaws(278) +ratelaws(361) -ratelaws(366) -ratelaws(391) -ratelaws(426) +ratelaws(471);
    Dspecies(117) = ratelaws(234) -ratelaws(279) -ratelaws(372) -ratelaws(380) +ratelaws(390) -ratelaws(395) -ratelaws(427) +ratelaws(462) +ratelaws(467);
    Dspecies(118) = ratelaws(239) -ratelaws(280) +ratelaws(308) -ratelaws(321) +ratelaws(360) -ratelaws(367) -ratelaws(402) -ratelaws(438) +ratelaws(475);
    Dspecies(119) = ratelaws(240) -ratelaws(281) +ratelaws(311) -ratelaws(322) -ratelaws(373) -ratelaws(381) +ratelaws(400) -ratelaws(407) -ratelaws(439) +ratelaws(463) +ratelaws(466);
    Dspecies(120) = ratelaws(241) -ratelaws(282) -ratelaws(362) +ratelaws(368) -ratelaws(403) +ratelaws(408) -ratelaws(410) -ratelaws(440) +ratelaws(478);
    Dspecies(121) = ratelaws(242) -ratelaws(283) +ratelaws(362) -ratelaws(368) -ratelaws(404) -ratelaws(441) +ratelaws(476);
    Dspecies(122) = ratelaws(243) -ratelaws(284) -ratelaws(374) -ratelaws(382) +ratelaws(403) -ratelaws(408) -ratelaws(442) +ratelaws(464) +ratelaws(468);
    Dspecies(123) = ratelaws(310) +ratelaws(379) +ratelaws(396) +ratelaws(416) -ratelaws(443) -ratelaws(455) -ratelaws(465) -ratelaws(472) -ratelaws(479);
    Dspecies(124) = ratelaws(312) +ratelaws(381) +ratelaws(409) +ratelaws(431) -ratelaws(444) -ratelaws(456) -ratelaws(466) -ratelaws(477) -ratelaws(484);
    Dspecies(125) = ratelaws(369) +ratelaws(388) +ratelaws(413) -ratelaws(445) -ratelaws(453) +ratelaws(457) -ratelaws(459) -ratelaws(469) -ratelaws(480);
    Dspecies(126) = ratelaws(370) +ratelaws(401) +ratelaws(428) -ratelaws(446) -ratelaws(454) +ratelaws(458) -ratelaws(460) -ratelaws(474) -ratelaws(485);
    Dspecies(127) = ratelaws(371) +ratelaws(389) +ratelaws(415) -ratelaws(447) +ratelaws(453) -ratelaws(457) -ratelaws(461) -ratelaws(470) -ratelaws(481);
    Dspecies(128) = ratelaws(372) +ratelaws(391) +ratelaws(417) -ratelaws(448) -ratelaws(462) -ratelaws(471) -ratelaws(482);
    Dspecies(129) = ratelaws(373) +ratelaws(402) +ratelaws(430) -ratelaws(449) +ratelaws(454) -ratelaws(458) -ratelaws(463) -ratelaws(475) -ratelaws(486);
    Dspecies(130) = ratelaws(374) +ratelaws(404) +ratelaws(432) -ratelaws(450) -ratelaws(464) -ratelaws(476) -ratelaws(487);
    Dspecies(131) = ratelaws(380) +ratelaws(397) +ratelaws(418) -ratelaws(451) -ratelaws(467) -ratelaws(473) -ratelaws(483);
    Dspecies(132) = ratelaws(382) +ratelaws(410) +ratelaws(433) -ratelaws(452) -ratelaws(468) -ratelaws(478) -ratelaws(488);

end


end
