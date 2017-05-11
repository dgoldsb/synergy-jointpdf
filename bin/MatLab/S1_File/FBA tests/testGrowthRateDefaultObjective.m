function [results]=testGrowthRateDefaultObjective(model, limiting)

    % This function compares model predictions of biomass flux to growth
    % rates as in doi:10.1186/1752-0509-7-36. This simulation uses a
    % minimal medium with glucose as the carbon source, the model-supplied
    % objective, and predicts biomass yields under carbon or nitrogen
    % limitations aerobically (most models do not predict anaerobic growth
    % without media supplementation).
    %
    % Input:
    %   model     A COBRA Toolbox-format yeast model. Currently supports:
    %              iFF708 - doi: 10.1101/gr.234503
    %              iND750 - doi: 10.1101/gr.2250904
    %              iIN800 - doi: 10.1186/1752-0509-2-71
    %              iMM904 - doi: 10.1186/1752-0509-3-37
    %              Yeast 4 - doi: 10.1186/1752-0509-4-145
    %              iAZ900 - doi: 10.1186/1752-0509-4-178
    %              iMM904bz - doi: 10.1038/ng.846
    %              Yeast 5 - doi:10.1186/1752-0509-6-55
    %              iTO977 - doi:10.1186/1752-0509-7-36
    %              Yeast 6 - doi:10.1093/database/bat059
    %              Yeast 7 - doi:10.1089/ind.2013.0013
    %              biomodels - doi:10.1186/1752-0509-7-116, downloaded from
    %                https://www.ebi.ac.uk/biomodels-main/BMID000000141353
    %
    %   limiting  1 for a single limiting nutrient, 2 for multiple
    %              constraints
    %
    %
    % Output:
    %   results   correlation (R) between reported experimental growth
    %             rates and calculated maximum achievable biomass flux 
    %             values

    %% citation
    % please cite: Heavner, Benjamin D. and Nathan Price. �Comparing Yeast
    % Metabolic Network Reconstructions� NEED TO ADD CITATION DETAILS
  
    lower(model.description)
    
    switch lower(model.description)
        case {'iff708_valid_modified.xml'}          
            model.c(1368) = 0; % the SBML file has two objectives set
            model = minimal_iFF(model);
            Glu_EX = 'GLCxtI';
            O2_EX = 'O2xtI';
            N_EX = 'NH3xtI';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, false)
  %              pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, false)
  %              pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, false)
  %              pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, false)
  %              pause
            end
            
        case {'ind750.xml'}
            model = minimal_iND(model);
            Glu_EX = 'EX_glc(e)'; 
            O2_EX = 'EX_o2(e)';
            N_EX = 'EX_nh4(e)';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, true)
 %               pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, true)
%                pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, true)
 %               pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, true)
 %               pause
            end
            
        case {'iin800_cobra'}
            model.c(findRxnIDs(model,'GROWTH')) = 1; % set growth objective
            
            model = addExchangeRxn(model,{'m758'},-1000,0);
            model = addExchangeRxn(model,{'m928'},-1000,0);
            
            % also includes TAG and phosphatidate to enable growth
            model = minimal_iIN(model);
            Glu_EX = 'GLCxtI';
            O2_EX = 'O2xtI';
            N_EX = 'NH3xtI';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, false)
  %              pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, false)
 %               pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, false)
 %               pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, false)
 %               pause
            end
                
        case {'iin800.xml'}
            model = minimal_iIN2(model);
%            model = set_iIN_iFF_biomass(model);
            Glu_EX = 'GLCxtI';
            O2_EX = 'O2xtI';
            N_EX = 'NH3xtI';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, false)
  %              pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, false)
 %               pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, false)
 %               pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, false)
 %               pause
            end
            
        case {'imm904'} % should have - constraints
            model = minimal_iMM(model);
 %           model = set_iMM_iFF_biomass(model);
            Glu_EX = 'EX_glc(e)'; 
            O2_EX = 'EX_o2(e)';
            N_EX = 'EX_nh4(e)';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, true)
  %              pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, true)
 %               pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, true)
  %              pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, true)
   %             pause
            end
            
        case {'yeast_4.05'} % nonstandard exchange bounds
            model = minimal_Y4(model);
%            model = set_Y4_iFF_biomass(model);
            Glu_EX = 'r_1293';
            O2_EX = 'r_1435';
            N_EX = 'r_1157';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, false)
  %              pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, false)
 %               pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, false)
 %               pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, false)
  %              pause
            end
            
        case {'iaz900.xml'}
           model = minimal_iAZ(model);
%           model = set_iAZ_iFF_biomass(model);
            Glu_EX = 'EX_glc(e)'; 
            O2_EX = 'EX_o2(e)';
            N_EX = 'EX_nh4(e)';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, true)
 %               pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, true)
%                pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, true)
 %               pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, true)
 %               pause
            end
           
        case {'imm904_nadcorrected'} % should have - constraints
            model = minimal_iMMbz(model);
%            model = set_iMM_NAD_iFF_biomass(model);
            Glu_EX = 'EX_glc(e)'; 
            O2_EX = 'EX_o2(e)';
            N_EX = 'EX_nh4(e)';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, true)
   %             pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, true)
  %              pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, true)
 %               pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, true)
 %               pause
            end
            
        case {'yeast_5.01_model.xml'}
            model = minimal_Y5(model);
%            model = set_Y5_iFF_biomass(model);
            Glu_EX = 'r_1714';
            O2_EX = 'r_1992';
            N_EX = 'r_1654';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, true)
  %              pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, true)
  %              pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, true)
  %              pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, true)
 %               pause
            end
            
        case {'ito977_cobra_compatible'}
            model = minimal_iTO(model);
%            model = set_iTO_iFF_biomass(model);
            Glu_EX = 'GLCxtI';
            O2_EX = 'O2xtI';
            N_EX = 'NH3xtI';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, false)
   %             pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, false)
   %             pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, false)
   %             pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, false)
   %             pause
            end
             
        case {'yeast_6.06_cobra'}
            model = minimal_Y6(model);
%            model = set_Y6_iFF_biomass(model);
            Glu_EX = 'r_1714';
            O2_EX = 'r_1992';
            N_EX = 'r_1654';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, true)
   %             pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, true)
   %             pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, true)
   %             pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, true)
    %            pause
            end

        case {'yeast_7.00_cobra'}
            model = minimal_Y6(model);
%            model = set_Y7_iFF_biomass(model);
            Glu_EX = 'r_1714';
            O2_EX = 'r_1992';
            N_EX = 'r_1654';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            if limiting == 1
                N_correlation = get_N_limited_correlation(model, exchangeIDs, true)
    %            pause
                C_correlation = get_C_limited_correlation(model, exchangeIDs, true)
    %            pause
            elseif limiting == 2
                mult_N_correlation = get_multiple_N_limited_correlation(model, exchangeIDs, true)
    %            pause
                mult_C_correlation = get_multiple_C_limited_correlation(model, exchangeIDs, true)
    %            pause
            end
            
        case {'bmid000000141353'}
            
            % first for growth, leave defaults of all uptakes allowed. Add
            % amino acids - the supplemental data says the model can't make
            % the following: L-tyrosine (met bigg_tyr_L_bm), L-lysine (met
            % bigg_lys_L_bm), L-isoleucine (met bigg_ile_L_bm), L-arginine
            % (met bigg_arg_L_bm), L-histidine (met bigg_his_L_bm),
            % L-methionine (met bigg_met_L_bm), and L-tryptophan
            % (bigg_trp_L_bm). (note, BY4743 has mutants in his3, leu2, and
            % ura, and something with lys2 and met15)


            % add L-tyrosine by allowing the bigg_tyr_L_out rxn to be
            % reversible
            model.lb(findRxnIDs(model,'bigg_tyr_L_out')) = -1000;

            % add L-lysine by allowing the met bigg_lys_L_out rxn to be
            % reversible
            model.lb(findRxnIDs(model,'bigg_lys_L_out')) = -1000;

            % add L-isoleucine by allowing the met bigg_ile_L_out rxn to be
            % reversible
            model.lb(findRxnIDs(model,'bigg_ile_L_out')) = -1000;

            % add L-arginine by allowing the met bigg_arg_L_out rxn to be
            % reversible
            model.lb(findRxnIDs(model,'bigg_arg_L_out')) = -1000;

            % add L-histidine by allowing the met bigg_his_L_out rxn to be
            % reversible
            model.lb(findRxnIDs(model,'bigg_his_L_out')) = -1000;

            % add L-methionine by allowing the met bigg_met_L_out rxn to be
            % reversible
            model.lb(findRxnIDs(model,'bigg_met_L_out')) = -1000;

            % add L-tryptophan by allowing the met bigg_trp_L_out rxn to be
            % reversible
            model.lb(findRxnIDs(model,'bigg_trp_L_out')) = -1000;
            
            % modify model gene annotation
            model.genes  = regexprep(model.genes,'_i','');
            model.genes  = regexprep(model.genes,'_\d','');
            metaCyc_genes = ~cellfun('isempty',regexp(model.genes,'MetaCyc'));
            model.genes(metaCyc_genes) = {' '};
            
            model = minimal_bio(model);
%            model = set_bio_iFF_biomass(model);
            Glu_EX = 'MNXM105_in';
            O2_EX = 'MNXM4_in';
            N_EX = 'MNXM15_in';
            exchangeIDs = {Glu_EX, O2_EX, N_EX};
            N_correlation = get_N_limited_correlation(model, exchangeIDs, true)
     %       pause
            C_correlation = get_C_limited_correlation(model, exchangeIDs, true)
    %        pause
    
        otherwise
            error('Sorry, this model is not currently supported.')
            
    end
end

%% get correlations
function correlation = get_N_limited_correlation(model, exhangeIDs, signchange)
% signchange is a logical; if TRUE, use negative bounds for exchange
% constriants

N_limited_aerobic_constraints = ...
        [5.8 2.7 .4; 
        4.83 4.42 .42;
        3.50 7.80 .61;
        4.61 9.20 .74;
        5.67 8.70 .83;
        8 8.80 .90;
        9.45 9.30 1.09;
        12.68 8.2 1.33];
    
    if signchange
        N_limited_aerobic_constraints = N_limited_aerobic_constraints * -1;
    end
    
    N_limited_observed = [.1 .1 .15 .18 .2 .25 .28 .34];
    N_limited_simulated = [];

    for index=1:length(N_limited_aerobic_constraints)
%        model.ub(findRxnIDs(model, exhangeIDs(1))) = N_limited_aerobic_constraints(index,1);
%        model.lb(findRxnIDs(model, exhangeIDs(1))) = N_limited_aerobic_constraints(index,1);

%        model.ub(findRxnIDs(model, exhangeIDs(2))) = N_limited_aerobic_constraints(index,2);
%        model.lb(findRxnIDs(model, exhangeIDs(2))) = N_limited_aerobic_constraints(index,2);

        model.ub(findRxnIDs(model, exhangeIDs(3))) = N_limited_aerobic_constraints(index,3);
        model.lb(findRxnIDs(model, exhangeIDs(3))) = N_limited_aerobic_constraints(index,3);

        sln = optimizeCbModel(model);
        N_limited_simulated(index) = sln.f;
        
    end

    N_correlation = corrcoef(N_limited_observed, N_limited_simulated);
    scatter(N_limited_observed, N_limited_simulated);
    correlation = N_correlation(1,2);
end

function correlation = get_C_limited_correlation(model, exhangeIDs, signchange)
% signchange is a logical; if TRUE, use negative bounds for exchange
% constriants

C_limited_aerobic_constraints = ...
    [1.15 2.7;
    1.17 2.5;
    1.69 4;
    2.26 5;
    2.88 6.5;
    3.27 7.46;
    3.29 7.8;
    3.88 8;
    6.20 7;
    7.89 6.5;
    13.39 3];
    
    if signchange
        C_limited_aerobic_constraints = C_limited_aerobic_constraints * -1;
    end

C_limited_observed = [.1 .1 .15 .2 .25 .27 .28 .31 .33 .35 .38];
C_limited_simulated = [];

    for index=1:length(C_limited_aerobic_constraints)
        model.ub(findRxnIDs(model, exhangeIDs(1))) = C_limited_aerobic_constraints(index,1);
        model.lb(findRxnIDs(model, exhangeIDs(1))) = C_limited_aerobic_constraints(index,1);

%        model.ub(findRxnIDs(model, exhangeIDs(2))) = C_limited_aerobic_constraints(index,2);
%        model.lb(findRxnIDs(model, exhangeIDs(2))) = C_limited_aerobic_constraints(index,2);

        sln = optimizeCbModel(model);
        C_limited_simulated(index) = sln.f;
    end

    C_correlation = corrcoef(C_limited_observed, C_limited_simulated);
    scatter(C_limited_observed, C_limited_simulated);
    correlation = C_correlation(1,2);
end

function correlation = get_multiple_N_limited_correlation(model, exhangeIDs, signchange)
% signchange is a logical; if TRUE, use negative bounds for exchange
% constriants

N_limited_aerobic_constraints = ...
        [5.8 2.7 .4; 
        4.83 4.42 .42;
        3.50 7.80 .61;
        4.61 9.20 .74;
        5.67 8.70 .83;
        8 8.80 .90;
        9.45 9.30 1.09;
        12.68 8.2 1.33];
    
    if signchange
        N_limited_aerobic_constraints = N_limited_aerobic_constraints * -1;
    end
    
    N_limited_observed = [.1 .1 .15 .18 .2 .25 .28 .34];
    N_limited_simulated = [];

    for index=1:length(N_limited_aerobic_constraints)
       model.ub(findRxnIDs(model, exhangeIDs(1))) = N_limited_aerobic_constraints(index,1);
       model.lb(findRxnIDs(model, exhangeIDs(1))) = N_limited_aerobic_constraints(index,1);

       model.ub(findRxnIDs(model, exhangeIDs(2))) = N_limited_aerobic_constraints(index,2);
       model.lb(findRxnIDs(model, exhangeIDs(2))) = N_limited_aerobic_constraints(index,2);

        model.ub(findRxnIDs(model, exhangeIDs(3))) = N_limited_aerobic_constraints(index,3);
        model.lb(findRxnIDs(model, exhangeIDs(3))) = N_limited_aerobic_constraints(index,3);

        sln = optimizeCbModel(model);
        N_limited_simulated(index) = sln.f;
        
    end

    N_correlation = corrcoef(N_limited_observed, N_limited_simulated);
    scatter(N_limited_observed, N_limited_simulated);
    correlation = N_correlation(1,2);
end

function correlation = get_multiple_C_limited_correlation(model, exhangeIDs, signchange)
% signchange is a logical; if TRUE, use negative bounds for exchange
% constriants

C_limited_aerobic_constraints = ...
    [1.15 2.7;
    1.17 2.5;
    1.69 4;
    2.26 5;
    2.88 6.5;
    3.27 7.46;
    3.29 7.8;
    3.88 8;
    6.20 7;
    7.89 6.5;
    13.39 3];
    
    if signchange
        C_limited_aerobic_constraints = C_limited_aerobic_constraints * -1;
    end

C_limited_observed = [.1 .1 .15 .2 .25 .27 .28 .31 .33 .35 .38];
C_limited_simulated = [];

    for index=1:length(C_limited_aerobic_constraints)
        model.ub(findRxnIDs(model, exhangeIDs(1))) = C_limited_aerobic_constraints(index,1);
        model.lb(findRxnIDs(model, exhangeIDs(1))) = C_limited_aerobic_constraints(index,1);

       model.ub(findRxnIDs(model, exhangeIDs(2))) = C_limited_aerobic_constraints(index,2);
       model.lb(findRxnIDs(model, exhangeIDs(2))) = C_limited_aerobic_constraints(index,2);

        sln = optimizeCbModel(model);
        C_limited_simulated(index) = sln.f;
    end

    C_correlation = corrcoef(C_limited_observed, C_limited_simulated);
    scatter(C_limited_observed, C_limited_simulated);
    correlation = C_correlation(1,2);
end

%% functions to set minimal medium
function model = minimal_iFF(model)
    % change iFF model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
    
    % start with a clean slate: set to unconstrained excretion, no uptake)

    uptake_rxn_ids = {'ZYMSTxtI';'XANxtI';'URIxtI';'UREAxtI';'URAxtI';...
        'THYxtI';'DTxtI';'THMxtI';'SLFxtI';'SUCxtI';'SUCCxtI';...
        'C180xtI';'SPRMxtI';'SORxtI';'NAxtI';'MMETxtI';'SAMxtI';...
        'RFLAVxtI';'PYRxtI';'KxtI';'PIMExtI';'PIxtI';'PEPTxtI';...
        'C160xtI';'O2xtI';'OGTxtI';'OPEPxtI';'NMNxtI';'NH3xtI';...
        'C140xtI';'MIxtI';'MTHNxtI';'MELIxtI';'MLTxtI';'MALxtI';...
        'VALxtI';'TYRxtI';'TRPxtI';'THRxtI';'SERxtI';'PROxtI';...
        'PHExtI';'ORNxtI';'METxtI';'LYSxtI';'LEUxtI';'ILExtI';'HISxtI';...
        'GLNxtI';'GLUxtI';'GLTxtI';'CYSxtI';'ASPxtI';'ASNxtI';'ARGxtI';...
        'ALAxtI';'INSxtI';'HYXNxtI';'GSNxtI';'GNxtI';'GLALxtI';...
        'GLYxtI';'GLxtI';'GLAMxtI';'GABAxtI';'FUMxtI';'FORxtI';'ETHxtI';...
        'ERGOxtI';'DTTPxtI';'RIBxtI';'RMNxtI';'MNTxtI';'DIPEPxtI';...
        'GA6PxtI';'GLACxtI';'FRUxtI';'DUxtI';'DINxtI';'DGxtI';'DCxtI';...
        'DAxtI';'ARABxtI';'CYTSxtI';'CYTDxtI';'CO2xtI';'CITxtI';...
        'CHOxtI';'BTxtI';'BMxtI';'FUCxtI';'MANxtI';'GLCxtI';'TRExtI';...
        'ATNxtI';'PAPxtI';'ADNxtI';'ADxtI';'ACxtI';'ACALxtI';'AONAxtI';...
        'DANNAxtI';'AKGxtI';'PNTOxtI';'LACxtI';};

    excretion_rxn_ids = {'ZYMSTxtO';'XANxtO';'URIxtO';'UREAxtO';...
        'URAxtO';'THYxtO';'DTxtO';'THMxtO';'SLFxtO';'SUCxtO';'SUCCxtO';...
        'C180xtO';'SPRMxtO';'SORxtO';'NAxtO';'MMETxtO';'SAMxtO';...
        'RFLAVxtO';'PYRxtO';'KxtO';'PIMExtO';'PIxtO';'PEPTxtO';...
        'C160xtO';'O2xtO';'OGTxtO';'OPEPxtO';'NMNxtO';'NH3xtO';...
        'C140xtO';'MIxtO';'MTHNxtO';'MELIxtO';'MLTxtO';'MALxtO';...
        'VALxtO';'TYRxtO';'TRPxtO';'THRxtO';'SERxtO';'PROxtO';...
        'PHExtO';'ORNxtO';'METxtO';'LYSxtO';'LEUxtO';'ILExtO';'HISxtO';...
        'GLNxtO';'GLUxtO';'GLTxtO';'CYSxtO';'ASPxtO';'ASNxtO';'ARGxtO';...
        'ALAxtO';'INSxtO';'HYXNxtO';'GSNxtO';'GNxtO';'GLALxtO';...
        'GLYxtO';'GLxtO';'GLAMxtO';'GABAxtO';'FUMxtO';'FORxtO';'ETHxtO';...
        'ERGOxtO';'DTTPxtO';'RIBxtO';'RMNxtO';'MNTxtO';'DIPEPxtO';...
        'GA6PxtO';'GLACxtO';'FRUxtO';'DUxtO';'DINxtO';'DGxtO';'DCxtO';...
        'DAxtO';'ARABxtO';'CYTSxtO';'CYTDxtO';'CO2xtO';'CITxtO';...
        'CHOxtO';'BTxtO';'FUCxtO';'MANxtO';'GLCxtO';'TRExtO';'ATNxtO';...
        'PAPxtO';'ADNxtO';'ADxtO';'ACxtO';'ACALxtO';'AONAxtO';...
        'DANNAxtO';'AKGxtO';'PNTOxtO';'LACxtO';'BMxtO';};

    uptake_rxn_indexes = findRxnIDs(model, uptake_rxn_ids);
    excretion_rxn_indexes = findRxnIDs(model, ...
        excretion_rxn_ids);

    model.lb(uptake_rxn_indexes) = 0;
    model.ub(uptake_rxn_indexes) = 0;
    model.lb(excretion_rxn_indexes) = 0;
    model.ub(excretion_rxn_indexes) = 1000;

    desiredExchanges = {...
        'NH3xtI';  ... % 'ammonium exchange';
        'PIxtI'; ... % 'phosphate exchange';
        'SLFxtI'; ... % 'sulphate exchange'; 
        %'NAxtI'; 
        %'KxtI'; 
        'O2xtI';... % 'oxygen exchange'; 
        };
    glucoseExchange = {...
        'GLCxtI'; ... % D-glucose exchange'
        };
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(uptakeRxnIndexes) ~= 4;
        error('Not all exchange reactions were found.')
    end
    model.ub(uptakeRxnIndexes) = 1000;
    model.ub(glucoseExchangeIndex) = 10;
end

function model = minimal_iND(model)
    % change iND model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
    
    % start with a clean slate: set all exchange reactions to upper bound =
    % 1000 and lower bound = 0 (ie, unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'EX_nh4(e)'; ... % 'ammonium exchange';
        'EX_o2(e)'; ... % 'oxygen exchange'; 
        'EX_pi(e)'; ... % 'phosphate exchange';
        'EX_so4(e)'; ... % 'sulphate exchange'; 
        };
    glucoseExchange = {...
        'EX_glc(e)'; ... % D-glucose exchange'
        };
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(uptakeRxnIndexes) ~= 4;
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
    model.lb(glucoseExchangeIndex)=-10;
end

function model = minimal_iMM(model)
    % change iMM model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
    
    % start with a clean slate: set all exchange reactions to upper bound =
    % 1000 and lower bound = 0 (ie, unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'EX_nh4(e)'; ... % 'ammonium exchange';
        'EX_o2(e)'; ... % 'oxygen exchange'; 
        'EX_pi(e)'; ... % 'phosphate exchange';
        'EX_so4(e)'; ... % 'sulphate exchange';
        };
    glucoseExchange = {...
        'EX_glc(e)'; ... % D-glucose exchange'
        };
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(uptakeRxnIndexes) ~= 4;
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
    model.lb(glucoseExchangeIndex)=-10;
end
   
function model = minimal_iMMbz(model)
    % change iMM model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
    
    % iMMbz version also requires Fe2 to enable heme biosynthesis
    
    % start with a clean slate: set all exchange reactions to upper bound =
    % 1000 and lower bound = 0 (ie, unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'EX_nh4(e)'; ... % 'ammonium exchange';
        'EX_o2(e)'; ... % 'oxygen exchange'; 
        'EX_pi(e)'; ... % 'phosphate exchange';
        'EX_so4(e)'; ... % 'sulphate exchange'; 
        'EX_fe2(e)'; ... % 'Fe2 exchange';
        };
    glucoseExchange = {...
        'EX_glc(e)'; ... % D-glucose exchange'
        };
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(uptakeRxnIndexes) ~= 5;
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
    model.lb(glucoseExchangeIndex)=-10;
end

function model = minimal_iAZ(model)
    % change iAZ model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
                
    % start with a clean slate: set all exchange reactions to upper bound =
    % 1000 and lower bound = 0 (ie, unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'EX_nh4(e)'; ... % 'ammonium exchange';
        'EX_o2(e)'; ... % 'oxygen exchange'; 
        'EX_pi(e)'; ... % 'phosphate exchange';
        'EX_so4(e)'; ... % 'sulphate exchange'; 
        };
    glucoseExchange = {...
        'EX_glc(e)'; ... % D-glucose exchange'
        };
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(uptakeRxnIndexes) ~= 4;
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
    model.lb(glucoseExchangeIndex)=-10;
end

function model = minimal_Y4(model)
    % change Y4 model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
    
    % This one is hard b/c nonstandard exchange bounds, directions.

    % start with a clean slate (ie, unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    exchange_indexes=find(exchangeRxns);

    %thanks to Shuyi Ma (of the Price lab at ISB) for this!
    b2 = sum(model.S(:,exchange_indexes),1); %-1 for exchange as reactant +1 for exchange as product
    posexch = exchange_indexes(b2 == 1); % index of exchange rxns having product
    negexch = exchange_indexes(b2 == -1); % index of exchange rxns having substrate

    model.lb(posexch) = -1000;
    model.ub(posexch) = 0;

    model.lb(negexch) = 0;
    model.ub(negexch) = 1000;

    % now allow uptake of media components
    desiredPosExchanges = {...
        'r_1157'; ... % 'ammonia exchange';
        'r_1435'; ... % 'oxygen exchange';
        'r_1461'; ... % 'phosphate exchange';
        'r_1507'; ... % 'sulfate uniport ';
        };
    glucosePosExchange = {...
        'r_1293'; ... % 'glucose transport' - should be ub=10 (was 0 10)
        };

    desiredPosExchangesIndexes = findRxnIDs(model,desiredPosExchanges);
    glucosePosExchangeIndex = findRxnIDs(model,glucosePosExchange);
    if length(desiredPosExchangesIndexes) ~= 4;
        error('Not all exchange reactions were found.')
    end
    model.ub(desiredPosExchangesIndexes)=1000;
    model.ub(glucosePosExchangeIndex)=10;
end

function model = minimal_Y5(model)
    % change Y5 model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
                
    % start with a clean slate: set all exchange reactions to upper bound =
    % 1000 and lower bound = 0 (ie, unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'r_1654'; ... % 'ammonium exchange';
        'r_1992'; ... % 'oxygen exchange'; 
        'r_2005'; ... % 'phosphate exchange';
        'r_2060'; ... % 'sulphate exchange'; 
        };
    glucoseExchange = {...
        'r_1714'; ... % 'D-glucose exchange'
        };
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(uptakeRxnIndexes) ~= 4;
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
    model.lb(glucoseExchangeIndex)=-10;
end

function model = minimal_Y6(model)
    % change Y6 model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
    
    % start with a clean slate: set all exchange reactions to
    % upper bound = 1000 and lower bound = 0 (ie, unconstrained
    % excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'r_1654'; ... % 'ammonium exchange';
        'r_1992'; ... % 'oxygen exchange'; 
        'r_2005'; ... % 'phosphate exchange';
        'r_2060'; ... % 'sulphate exchange';
        %'r_1861'; ... % iron for test of expanded biomass def;
        };
    glucoseExchange = {...
        'r_1714'; ... % D-glucose exchange'
        };
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(uptakeRxnIndexes) ~= 4;%5;
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
    model.lb(glucoseExchangeIndex)=-10;
end
    
function model = minimal_iIN(model)
    % change iIN model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
    
    % start with a clean slate: set all exchange reactions to
    % unconstrained excretion, no uptake.
    
    % iIN also requires TAG and phosphatidate to enable growth
    
    uptake_rxn_ids = {'44DIMZYMSTxtI';'ACALxtI';'ACxtI';'ADNxtI';...
        'ADxtI';'AKGxtI';'ALAVxtI';'ALAxtI';'AMGxtI';'AONAxtI';...
        'ARABxtI';'ARGxtI';'ASNxtI';'ASPxtI';'ATNxtI';'ATTxtI';...
        'BMxtI';'BTxtI';'C100xtI';'C120xtI';'C140xtI';'C141xtI';...
        'C160xtI';'C161xtI';'C180xtI';'C181xtI';'C24xtI';'C26xtI';...
        'CARxtI';'CHOxtI';'CITxtI';'CO2xtI';'CYSxtI';'CYTDxtI';...
        'CYTSxtI';'DANNAxtI';'DAxtI';'DCxtI';'DGxtI';'DINxtI';...
        'DIPEPxtI';'DTTPxtI';'DTxtI';'DUxtI';'EPISTxtI';...
        'ERG572224xtI';'ERG722xtI';'ERGOSTxtI';'ETHxtI';...
        'FCOSTxtI';'FORxtI';'FRUxtI';'FUCxtI';'FUMxtI';'GA6PxtI';...
        'GABAxtI';'GLACxtI';'GLALxtI';'GLAMxtI';'GLCxtI';'GLNxtI';...
        'GLTxtI';'GLUxtI';'GLxtI';'GLYxtI';'GNxtI';'GROPCxtI';...
        'GROPIxtI';'GSNxtI';'HEXTI';'HISxtI';'HYXNxtI';'ILExtI';...
        'INSxtI';'KxtI';'LACxtI';'LANOSTxtI';'LEUxtI';'LYSxtI';...
        'MALxtI';'MANxtI';'MELIxtI';'METxtI';'MIxtI';'MLTxtI';...
        'MMETxtI';'MNTxtI';'MTHNxtI';'NAGxtI';'NAxtI';'NH3xtI';...
        'NMNxtI';'O2xtI';'OGTxtI';'OPEPxtI';'ORNxtI';'PAPxtI';...
        'PEPTxtI';'PHExtI';'PIMExtI';'PIxtI';'PNTOxtI';'PROxtI';...
        'PTRSCxtI';'PYRxtI';'RFLAVxtI';'RIBxtI';'RMNxtI';'SAMxtI';...
        'SERxtI';'SLFxtI';'SORxtI';'SPRMDxtI';'SPRMxtI';'SUCCxtI';...
        'SUCxtI';'THMxtI';'THRxtI';'THYxtI';'TRExtI';'TRPxtI';...
        'TYRxtI';'URAxtI';'UREAxtI';'URIxtI';'VALxtI';'VB6xtI';...
        'XANxtI';'XTSINExtI';'XYLxtI';'ZYMSTxtI';};

    excretion_rxn_ids = {'44DIMZYMSTxtO';'ACALxtO';'ACxtO';...
        'ADNxtO';'ADxtO';'AKGxtO';'ALAVxtO';'ALAxtO';'AMGxtO';...
        'AONAxtO';'ARABxtO';'ARGxtO';'ASNxtO';'ASPxtO';'ATNxtO';...
        'ATTxtO';'BMxtO';'BTxtO';'C100xtO';'C120xtO';'C140xtO';...
        'C141xtO';'C160xtO';'C161xtO';'C180xtO';'C181xtO';...
        'C24xtO';'C26xtO';'CARxtO';'CHOxtO';'CITxtO';'CO2xtO';...
        'CYSxtO';'CYTDxtO';'CYTSxtO';'DANNAxtO';'DAxtO';'DCxtO';...
        'DGxtO';'DINxtO';'DIPEPxtO';'DTTPxtO';'DTxtO';'DUxtO';...
        'EPISTxtO';'ERG572224xtO';'ERG722xtO';'ERGOSTxtO';...
        'ETHxtO';'FCOSTxtO';'FORxtO';'FRUxtO';'FUCxtO';'FUMxtO';...
        'GA6PxtO';'GABAxtO';'GLACxtO';'GLALxtO';'GLAMxtO';...
        'GLCxtO';'GLNxtO';'GLTxtO';'GLUxtO';'GLxtO';'GLYxtO';...
        'GNxtO';'GROPCxtO';'GROPIxtO';'GSNxtO';'HEXTO';'HISxtO';...
        'HYXNxtO';'ILExtO';'INSxtO';'KxtO';'LACxtO';'LANOSTxtO';...
        'LEUxtO';'LYSxtO';'MALxtO';'MANxtO';'MELIxtO';'METxtO';...
        'MIxtO';'MLTxtO';'MMETxtO';'MNTxtO';'MTHNxtO';'NAGxtO';...
        'NAxtO';'NH3xtO';'NMNxtO';'O2xtO';'OGTxtO';'OPEPxtO';...
        'ORNxtO';'PAPxtO';'PEPTxtO';'PHExtO';'PIMExtO';'PIxtO';...
        'PNTOxtO';'PROxtO';'PTRSCxtO';'PYRxtO';'RFLAVxtO';...
        'RIBxtO';'RMNxtO';'SAMxtO';'SERxtO';'SLFxtO';'SORxtO';...
        'SPRMDxtO';'SPRMxtO';'SUCCxtO';'SUCxtO';'THMxtO';'THRxtO';...
        'THYxtO';'TRExtO';'TRPxtO';'TYRxtO';'URAxtO';'UREAxtO';...
        'URIxtO';'VALxtO';'VB6xtO';'XANxtO';'XTSINExtO';'XYLxtO';...
        'ZYMSTxtO';};

        uptake_rxn_indexes = findRxnIDs(model,uptake_rxn_ids);
        excretion_rxn_indexes = findRxnIDs(model,excretion_rxn_ids);
        
        model.lb(uptake_rxn_indexes)=0;
        model.ub(uptake_rxn_indexes)=0;
        
        model.lb(excretion_rxn_indexes)=0;
        model.ub(excretion_rxn_indexes)=1000;
        
    % now allow uptake of media components
    desiredExchanges = {...
        'NH3xtI'; ... % 'ammonia exchange';
        'O2xtI'; ... % 'oxygen exchange';
        'PIxtI'; ... % 'phosphate exchange';
        'SLFxtI'; ... % 'sulfate uniport ';       
        };
    glucoseExchange = {...
        'GLCxtI'; ... % 'glucose transport' - should be ub=10 (was 0 10)
        };

    desiredExchangesIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(desiredExchangesIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.ub(desiredExchangesIndexes)=1000;
    model.ub(glucoseExchangeIndex)=10;
end

function model = minimal_iIN2(model)
    % change iIN model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
    
    % start with a clean slate: set all exchange reactions to
    % unconstrained excretion, no uptake.
       
    uptake_rxn_ids = {'44DIMZYMSTxtI';'ACALxtI';'ACxtI';'ADNxtI';...
        'ADxtI';'AKGxtI';'ALAVxtI';'ALAxtI';'AMGxtI';'AONAxtI';...
        'ARABxtI';'ARGxtI';'ASNxtI';'ASPxtI';'ATNxtI';'ATTxtI';...
        'BMxtI';'BTxtI';'C100xtI';'C120xtI';'C140xtI';'C141xtI';...
        'C160xtI';'C161xtI';'C180xtI';'C181xtI';'C24xtI';'C26xtI';...
        'CARxtI';'CHOxtI';'CITxtI';'CO2xtI';'CYSxtI';'CYTDxtI';...
        'CYTSxtI';'DANNAxtI';'DAxtI';'DCxtI';'DGxtI';'DINxtI';...
        'DIPEPxtI';'DTTPxtI';'DTxtI';'DUxtI';'EPISTxtI';...
        'ERG572224xtI';'ERG722xtI';'ERGOSTxtI';'ETHxtI';...
        'FCOSTxtI';'FORxtI';'FRUxtI';'FUCxtI';'FUMxtI';'GA6PxtI';...
        'GABAxtI';'GLACxtI';'GLALxtI';'GLAMxtI';'GLCxtI';'GLNxtI';...
        'GLTxtI';'GLUxtI';'GLxtI';'GLYxtI';'GNxtI';'GROPCxtI';...
        'GROPIxtI';'GSNxtI';'HEXTI';'HISxtI';'HYXNxtI';'ILExtI';...
        'INSxtI';'KxtI';'LACxtI';'LANOSTxtI';'LEUxtI';'LYSxtI';...
        'MALxtI';'MANxtI';'MELIxtI';'METxtI';'MIxtI';'MLTxtI';...
        'MMETxtI';'MNTxtI';'MTHNxtI';'NAGxtI';'NAxtI';'NH3xtI';...
        'NMNxtI';'O2xtI';'OGTxtI';'OPEPxtI';'ORNxtI';'PAPxtI';...
        'PEPTxtI';'PHExtI';'PIMExtI';'PIxtI';'PNTOxtI';'PROxtI';...
        'PTRSCxtI';'PYRxtI';'RFLAVxtI';'RIBxtI';'RMNxtI';'SAMxtI';...
        'SERxtI';'SLFxtI';'SORxtI';'SPRMDxtI';'SPRMxtI';'SUCCxtI';...
        'SUCxtI';'THMxtI';'THRxtI';'THYxtI';'TRExtI';'TRPxtI';...
        'TYRxtI';'URAxtI';'UREAxtI';'URIxtI';'VALxtI';'VB6xtI';...
        'XANxtI';'XTSINExtI';'XYLxtI';'ZYMSTxtI';};

    excretion_rxn_ids = {'44DIMZYMSTxtO';'ACALxtO';'ACxtO';...
        'ADNxtO';'ADxtO';'AKGxtO';'ALAVxtO';'ALAxtO';'AMGxtO';...
        'AONAxtO';'ARABxtO';'ARGxtO';'ASNxtO';'ASPxtO';'ATNxtO';...
        'ATTxtO';'BMxtO';'BTxtO';'C100xtO';'C120xtO';'C140xtO';...
        'C141xtO';'C160xtO';'C161xtO';'C180xtO';'C181xtO';...
        'C24xtO';'C26xtO';'CARxtO';'CHOxtO';'CITxtO';'CO2xtO';...
        'CYSxtO';'CYTDxtO';'CYTSxtO';'DANNAxtO';'DAxtO';'DCxtO';...
        'DGxtO';'DINxtO';'DIPEPxtO';'DTTPxtO';'DTxtO';'DUxtO';...
        'EPISTxtO';'ERG572224xtO';'ERG722xtO';'ERGOSTxtO';...
        'ETHxtO';'FCOSTxtO';'FORxtO';'FRUxtO';'FUCxtO';'FUMxtO';...
        'GA6PxtO';'GABAxtO';'GLACxtO';'GLALxtO';'GLAMxtO';...
        'GLCxtO';'GLNxtO';'GLTxtO';'GLUxtO';'GLxtO';'GLYxtO';...
        'GNxtO';'GROPCxtO';'GROPIxtO';'GSNxtO';'HEXTO';'HISxtO';...
        'HYXNxtO';'ILExtO';'INSxtO';'KxtO';'LACxtO';'LANOSTxtO';...
        'LEUxtO';'LYSxtO';'MALxtO';'MANxtO';'MELIxtO';'METxtO';...
        'MIxtO';'MLTxtO';'MMETxtO';'MNTxtO';'MTHNxtO';'NAGxtO';...
        'NAxtO';'NH3xtO';'NMNxtO';'O2xtO';'OGTxtO';'OPEPxtO';...
        'ORNxtO';'PAPxtO';'PEPTxtO';'PHExtO';'PIMExtO';'PIxtO';...
        'PNTOxtO';'PROxtO';'PTRSCxtO';'PYRxtO';'RFLAVxtO';...
        'RIBxtO';'RMNxtO';'SAMxtO';'SERxtO';'SLFxtO';'SORxtO';...
        'SPRMDxtO';'SPRMxtO';'SUCCxtO';'SUCxtO';'THMxtO';'THRxtO';...
        'THYxtO';'TRExtO';'TRPxtO';'TYRxtO';'URAxtO';'UREAxtO';...
        'URIxtO';'VALxtO';'VB6xtO';'XANxtO';'XTSINExtO';'XYLxtO';...
        'ZYMSTxtO';};

        uptake_rxn_indexes = findRxnIDs(model,uptake_rxn_ids);
        excretion_rxn_indexes = findRxnIDs(model,excretion_rxn_ids);
        
        model.lb(uptake_rxn_indexes)=0;
        model.ub(uptake_rxn_indexes)=0;
        
        model.lb(excretion_rxn_indexes)=0;
        model.ub(excretion_rxn_indexes)=1000;
        
    % now allow uptake of media components
    desiredExchanges = {...
        'NH3xtI'; ... % 'ammonia exchange';
        'O2xtI'; ... % 'oxygen exchange';
        'PIxtI'; ... % 'phosphate exchange';
        'SLFxtI'; ... % 'sulfate uniport ';       
        };
    glucoseExchange = {...
        'GLCxtI'; ... % 'glucose transport' - should be ub=10 (was 0 10)
        };

    desiredExchangesIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(desiredExchangesIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.ub(desiredExchangesIndexes)=1000;
    model.ub(glucoseExchangeIndex)=10;
end

function model = minimal_iTO(model)
    % change iTO model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate
    
    % start with a clean slate: set all exchange reactions to
    % unconstrained excretion, no uptake.
    
    uptake_rxn_ids = {'2MBACxtI';'2MBALDxtI';'44DIMZYMSTxtI';...
        'ACALxtI';'ACxtI';'ADNxtI';'ADxtI';'AKGxtI';'ALAVxtI';'ALAxtI';...
        'AMGxtI';'AONAxtI';'ARABxtI';'ARGxtI';'ASNxtI';'ASPxtI';...
        'ATNxtI';'ATTxtI';'BMxtI';'BTxtI';'C100xtI';'C120xtI';...
        'C140xtI';'C141xtI';'C160xtI';'C161xtI';'C180xtI';'C181xtI';...
        'C24xtI';'C26xtI';'CARxtI';'CHOxtI';'CITxtI';'CO2xtI';'CYSxtI';...
        'CYTDxtI';'CYTSxtI';'DANNAxtI';'DAxtI';'DCxtI';'DGxtI';...
        'DINxtI';'DIPEPxtI';'DTTPxtI';'DTxtI';'DUxtI';'EPISTxtI';...
        'ERG572224xtI';'ERG722xtI';'ERGOSTxtI';'ETHAxtI';'ETHxtI';...
        'FCOSTxtI';'FORxtI';'FRUxtI';'FUCxtI';'FUMxtI';'GA6PxtI';...
        'GABAxtI';'GLACxtI';'GLALxtI';'GLAMxtI';'GLCxtI';'GLNxtI';...
        'GLTxtI';'GLUxtI';'GLxtI';'GLYxtI';'GNxtI';'GROPCxtI';...
        'GROPIxtI';'GSNxtI';'HEXTI';'HISxtI';'HYXNxtI';'ILExtI';...
        'INSxtI';'KxtI';'LACxtI';'LANOSTxtI';'LEUxtI';'LYSxtI';...
        'MALxtI';'MANxtI';'MELIxtI';'METxtI';'MIxtI';'MLTxtI';...
        'MMETxtI';'MNTxtI';'MTHNxtI';'NAGxtI';'NAxtI';'NH3xtI';...
        'NMNxtI';'O2xtI';'OGTxtI';'OPEPxtI';'ORNxtI';'PAPxtI';...
        'PEPTxtI';'PHExtI';'PIMExtI';'PIxtI';'PNTOxtI';'PROxtI';...
        'PTRSCxtI';'PYRxtI';'RFLAVxtI';'RIBxtI';'RMNxtI';'SAMxtI';...
        'SERxtI';'SLFxtI';'SORxtI';'SPRMDxtI';'SPRMxtI';'SUCCxtI';...
        'SUCxtI';'THMxtI';'THRxtI';'THYxtI';'TRExtI';'TRPxtI';'TYRxtI';...
        'URAxtI';'UREAxtI';'URIxtI';'VALxtI';'VB6xtI';'XANxtI';...
        'XTSINExtI';'XYLxtI';'ZYMSTxtI';'DSERxtI'};
    
    excretion_rxn_ids = {'2MBACxtO';'44DIMZYMSTxtO';'ACALxtO';'ACxtO';...
        'ADNxtO';'ADxtO';'AKGxtO';'ALAVxtO';'ALAxtO';'AMGxtO';...
        'AONAxtO';'ARABxtO';'ARGxtO';'ASNxtO';'ASPxtO';'ATNxtO';...
        'ATTxtO';'BMxtO';'BTxtO';'C100xtO';'C120xtO';'C140xtO';...
        'C141xtO';'C160xtO';'C161xtO';'C180xtO';'C181xtO';'C24xtO';...
        'C26xtO';'CARxtO';'CHOxtO';'CITxtO';'CO2xtO';'CYSxtO';...
        'CYTDxtO';'CYTSxtO';'DANNAxtO';'DAxtO';'DCxtO';'DGxtO';...
        'DINxtO';'DIPEPxtO';'DTTPxtO';'DTxtO';'DUxtO';'EPISTxtO';...
        'ERG572224xtO';'ERG722xtO';'ERGOSTxtO';'ETHxtO';'FCOSTxtO';...
        'FORxtO';'FRUxtO';'FUCxtO';'FUMxtO';'GA6PxtO';'GABAxtO';...
        'GLACxtO';'GLALxtO';'GLAMxtO';'GLCxtO';'GLNxtO';'GLTxtO';...
        'GLUxtO';'GLxtO';'GLYxtO';'GNxtO';'GROPCxtO';'GROPIxtO';...
        'GSNxtO';'HEXTO';'HISxtO';'HYXNxtO';'ILExtO';'INSxtO';'KxtO';...
        'LACxtO';'LANOSTxtO';'LEUxtO';'LYSxtO';'MALxtO';'MANxtO';...
        'MELIxtO';'METxtO';'MIxtO';'MLTxtO';'MMETxtO';'MNTxtO';...
        'MTHNxtO';'NAGxtO';'NAxtO';'NH3xtO';'NMNxtO';'O2xtO';'OGTxtO';...
        'OPEPxtO';'ORNxtO';'PAPxtO';'PEPTxtO';'PHExtO';'PIMExtO';...
        'PIxtO';'PNTOxtO';'PROxtO';'PTRSCxtO';'PYRxtO';'RFLAVxtO';...
        'RIBxtO';'RMNxtO';'SAMxtO';'SERxtO';'SLFxtO';'SORxtO';...
        'SPRMDxtO';'SPRMxtO';'SUCCxtO';'SUCxtO';'THMxtO';'THRxtO';...
        'THYxtO';'TRExtO';'TRPxtO';'TYRxtO';'URAxtO';'UREAxtO';...
        'URIxtO';'VALxtO';'VB6xtO';'XANxtO';'XTSINExtO';'XYLxtO';...
        'ZYMSTxtO';'DSERxtO'};

        uptake_rxn_indexes = findRxnIDs(model,uptake_rxn_ids);
        excretion_rxn_indexes = findRxnIDs(model,excretion_rxn_ids);
        
        model.lb(uptake_rxn_indexes)=0;
        model.ub(uptake_rxn_indexes)=0;
        
        model.lb(excretion_rxn_indexes)=0;
        model.ub(excretion_rxn_indexes)=1000;
        
    % now allow uptake of media components
    desiredExchanges = {...
        'NH3xtI'; ... % 'ammonia exchange';
        'O2xtI'; ... % 'oxygen exchange';
        'PIxtI'; ... % 'phosphate exchange';
        'SLFxtI'; ... % 'sulfate uniport ';       
        };
    glucoseExchange = {...
        'GLCxtI'; ... % 'glucose transport' - should be ub=10 (was 0 10)
        };

    desiredExchangesIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    if length(desiredExchangesIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.ub(desiredExchangesIndexes)=1000;
    model.ub(glucoseExchangeIndex)=10;
end

function model = minimal_bio(model)
    % change bio model media to minimal - ammonium, glucose, oxygen,
    % phosphate, sulphate. Add L-tyrosine (met bigg_tyr_L_bm), L-lysine
    % (met bigg_lys_L_bm), L-isoleucine (met bigg_ile_L_bm), L-arginine
    % (met bigg_arg_L_bm), L-histidine (met bigg_his_L_bm), L-methionine
    % (met bigg_met_L_bm), and L-tryptophan (bigg_trp_L_bm) to enable
    % growth.
    
    % start with a clean slate: set uptake exchange reactions to
    % upper bound = 0 and lower bound = 0 (ie, unconstrained
    % excretion, no uptake)
 
    uptakeRxns = findRxnIDs(model, ...
        {'MNXR6704_i'; 'MNXM99_in'; 'MNXM105_in'; 'MNXM15_in'; ...
        'MNXM27_in'; 'MNXM95_in'; 'MNXM653_in'; 'MNXM128_in'; ...
        'MNXM58_in'; 'MNXM43_in'; 'MNXM9_in'; 'MNXM1_in'; 'MNXM2_in'; ...
        'MNXM13_in'; 'MNXM4_in';});

    otherRxns = setdiff(find(findExcRxns(model)), uptakeRxns);

    model.ub(uptakeRxns) = 0;
    model.ub(otherRxns) = 1000;
    model.lb(otherRxns) = 0;
    
    desiredExchanges = {...
        'MNXM15_in'; ... % 'ammonium exchange';
        'MNXM4_in'; ... % 'oxygen exchange'; 
        'MNXM9_in'; ... % 'phosphate exchange';
        'MNXM58_in'; ... % 'sulphate exchange';
        };
    glucoseExchange = {...
        'MNXM105_in'; ... % D-glucose exchange'
        };
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    
    if length(uptakeRxnIndexes) ~= 4;%5;
        error('Not all exchange reactions were found.')
    end
    
    model.ub(uptakeRxnIndexes)=1000;
    model.ub(glucoseExchangeIndex)=10;
    
    % add L-tyrosine by allowing the bigg_tyr_L_out rxn to be reversible
    model.lb(findRxnIDs(model,'bigg_tyr_L_out')) = -1000;
    
    % add L-lysine by allowing the met bigg_lys_L_out rxn to be reversible
    model.lb(findRxnIDs(model,'bigg_lys_L_out')) = -1000;

    % add L-isoleucine by allowing the met bigg_ile_L_out rxn to be
    % reversible
    model.lb(findRxnIDs(model,'bigg_ile_L_out')) = -1000;

    % add L-arginine by allowing the met bigg_arg_L_out rxn to be
    % reversible
    model.lb(findRxnIDs(model,'bigg_arg_L_out')) = -1000;

    % add L-histidine by allowing the met bigg_his_L_out rxn to be
    % reversible
    model.lb(findRxnIDs(model,'bigg_his_L_out')) = -1000;
    
    % add L-methionine by allowing the met bigg_met_L_out rxn to be
    % reversible
    model.lb(findRxnIDs(model,'bigg_met_L_out')) = -1000;
    
    % add L-tryptophan by allowing the met bigg_trp_L_out rxn to be
    % reversible
    model.lb(findRxnIDs(model,'bigg_trp_L_out')) = -1000;
end

