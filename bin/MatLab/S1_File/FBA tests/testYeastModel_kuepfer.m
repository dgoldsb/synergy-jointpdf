function [results]=testYeastModel_kuepfer(model,varargin)

    % testYeastModel_kuepfer tests model predictions against the genes
    % screened by Keupfer et al, using the media they used and glucose as
    % the carbon source. It is reliant on the COBRA Toolbox.
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
    %   biomass (optional)
    %             0 = use model-defined biomass def (default)
    %             1 = use iFF708 biomass definition
    %
    %   carbon source (optional)
    %             0 = glucose (default)
    %             1 = galactose
    %             2 = glycerol
    %             3 = ethanol
    %
    %   output (optional)
    %             0 = silent
    %             1 = default screen output
    %             2 = verbose screen output
    %
    % Output:
    %   screen output: description of model statistics and predictive
    %             accuracy compared to lists of essential genes and
    %             auxotrophs reported by Kuepfer et al.
    %
    %   results   a structure containing lists of true positive, true
    %             negative, false positive, and false negative results

    %% citation
    % based on testYeastmodel code by kieran smallbone and ben heavner,
    % doi: 10.1186/1752-0509-6-55
    %
    % please cite: Heavner, Benjamin D. and Nathan Price. �Comparative
    % Analysis of Yeast Metabolic Network Models Highlights Progress,
    % Opportunity� NEED TO ADD CITATION DETAILS

    %% process input arguments
    % 3 optional inputs (4 total) at most
    numvarargs = length(varargin);
    if numvarargs > 3
        error('myfuns:testYeast:TooManyInputs', ...
            'requires at most 3 optional inputs');
    end

    % set defaults for optional inputs
    optargs = {0 0 1};

    % put defaults values into into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;

    % Place optional args in memorable variable names
    [biomass carbon output] = optargs{:};
    
    switch lower(model.description)
        case {'iff708_valid_modified.xml'} 
            if output
                fprintf('\nSetting iFF708 model objective function.\n');
            end
            
            model.c(1368) = 0; % the SBML file has two objectives set
            
        case {'ind750.xml'}
            
        case {'iin800_cobra'}
            if output
                fprintf('\nSetting iIN800 model objective function.\n\n');
            end
            
            model.c(findRxnIDs(model,'GROWTH')) = 1; % set growth objective
            
            if output
                fprintf('\nAdding TAG and phosphatidate exchanges.\n\n');
            end
            
            model = addExchangeRxn(model,{'m758'},-1000,0);
            model = addExchangeRxn(model,{'m928'},-1000,0);
            
        case {'iin800.xml'}
                        
        case {'imm904'}
            
        case {'yeast_4.05'}            
            
        case {'iaz900.xml'}
           
        case {'imm904_nadcorrected'}
            
        case {'yeast_5.01_model.xml'}
            
        case {'ito977_cobra_compatible'}
             
        case {'yeast_6.06_cobra'}           

        case {'yeast_7.00_cobra'}
            
        case {'bmid000000141353'}
            if output
                fprintf('\nModifying biodb model to enable FBA.\n');
            end
            
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
    
        otherwise
            error('Sorry, this model is not currently supported.')
            
    end

    %% set tolerance for KO growth
    ko_tol = 1e-6;

    %% simple model statistics

    if output
        % special handling for biomodelsdb, which has different gene
        % formatting
        if strcmpi(model.description, 'bmid000000141353')
            gene_number = length(model.genes) - ...
                sum(strncmp(' ',model.genes,1));
            
            fprintf('\n%s\n\t%s\t%s\n\n\t%g\t%s\n\t%g\t%s\n\t%4g\t%s\n\n',...
                'Model Description:', 'name:', model.description,...
                gene_number, 'genes', ...
                length(model.mets), 'metabolites',...
                length(model.rxns), 'reactions');
        else
        fprintf('\n%s\n\t%s\t%s\n\n\t%g\t%s\n\t%g\t%s\n\t%4g\t%s\n\n',...
            'Model Description:', 'name:', model.description,...
            length(model.genes), 'genes', ...
            length(model.mets), 'metabolites',...
            length(model.rxns), 'reactions');
        end
    end

    %% dubious ORFs

    dubiousORFs = setdiff(model.genes,verifiedORFs);

    if output
        fprintf('\t%4g\t(%.2f%%)\t%s\n',length(dubiousORFs),...
            (100*length(dubiousORFs)/length(model.genes)),'dubious ORFs');
    end

    if output == 2
        fprintf('\n\t%s\n','list ORFs included in the model but annotated as dubious by SGD:');
        for k = 1:length(dubiousORFs)
            fprintf('\t\t%s\n',dubiousORFs{k});
        end
    end

    %% create model with kuepfer medium
    %     
    % kuepfer media: "minimal medium was, per liter (Verduyn et al.
    % 1992),5g of (NH4)2 SO4,3gofKH2PO4, 0.5 g of MgSO4�7H2O, 4.5 mg of
    % ZnSO4�7H2O, 0.3 mg of CoCl2�6H2O, 1.0 mg of MnCl2�4H2O, 0.3 mg of
    % CuSO4�5H2O, 4.5 mg of CaCl2�2H2O, 3.0 mg of FeSO4�7H2O, 0.4 mg of
    % NaMoO4�2H2O, 1.0 mg of H3BO3, 0.1mg of KI, 15 mg of EDTA, 0.05 mg of
    % biotin, 1.0 mg of Ca-pantothenate, 1.0 mg of nicotinic acid, 25 mg of
    % inositol, 1.0 mg of pyridoxine, 0.2 mg of p-amino-benzoic acid, and
    % 1.0 mg ofthiamine. The carbon sources (ethanol, galactose, glucose,
    % and glycerol) were added to a final concentration of 20 g/L. Strain
    % auxotrophies were complemented with 20 mg/L histidine, uracil,
    % methionine, and 60 mg/L leucine. About 50 strains of the yeast
    % collection are lysine auxotroph and were independently tested for
    % growth on plates supplemented with 20 mg/L lysine."

    if output
        fprintf('\nMaking Kuepfer medium...\n');
    end

    switch lower(model.description)
        case {'iff708_valid_modified.xml'} 
            model = kuepfer_iFF(model);

        case {'ind750.xml'}
            model = kuepfer_iND(model);

        case {'iin800_cobra', 'iin800.xml'}
            % also includes TAG and phosphatidate to enable growth
            model = kuepfer_iIN(model);
        
        case {'imm904'}
            model = kuepfer_iMM(model);

       case {'yeast_4.05'}
            model = kuepfer_Y4(model);

       case {'iaz900.xml'}
            model = kuepfer_iAZ(model);

        case {'imm904_nadcorrected'}
            model = kuepfer_iMMbs(model);

        case {'yeast_5.01_model.xml'}
            model = kuepfer_Y5(model);

        case {'ito977_cobra_compatible'}
            model = kuepfer_iTO(model);

        case {'yeast_6.06_cobra'}
            model = kuepfer_Y6(model);
            
        case {'yeast_7.00_cobra'}
            model = kuepfer_Y7(model);

        case {'bmid000000141353'}
            model = kuepfer_bio(model);

        otherwise
            error('Sorry, this model is not currently supported.')
    end

    %% add carbon source
    
    if output
        fprintf('\nSetting carbon source...\n');
    end
        
    if carbon == 0 % glucose
        switch lower(model.description)
            case {'iff708_valid_modified.xml'} 
                carbonExchange = {...
                    'GLCxtI'; ... % D-glucose exchange'
                };
            carbonExchangeIndex = findRxnIDs(model,carbonExchange);
            model.ub(carbonExchangeIndex) = 10;

            case {'ind750.xml'} 
                carbonExchange = {...
                    'EX_glc(e)'; ... % glucose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'imm904'}
                carbonExchange = {...
                    'EX_glc(e)'; ... % glucose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'iaz900.xml'}
                carbonExchange = {...
                    'EX_glc(e)'; ... % glucose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'yeast_4.05'}
                carbonExchange = {...
                    'r_1293'; ... % glucose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.ub(carbonExchangeIndex) = 10; % checked - this is right

            case {'yeast_5.01_model.xml'}
                carbonExchange = {...
                    'r_1714'; ... % glucose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'yeast_6.06_cobra'}
                carbonExchange = {...
                    'r_1714'; ... % glucose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;
                           
            case {'imm904_nadcorrected'}
                carbonExchange = {...
                    'EX_glc(e)'; ... % glucose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'iin800_cobra', 'iin800.xml', 'ito977_cobra_compatible'}
                carbonExchange = {...
                    'GLCxtI'; ... % glucose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.ub(carbonExchangeIndex) = 10;
                
            case {'yeast_7.00_cobra'}
                carbonExchange = {...
                    'r_1714'; ... % glucose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;
                
            case {'bmid000000141353'}
                carbonExchange = {...
                    'MNXM99_in'; ... % glucose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.ub(carbonExchangeIndex) = 10;

            otherwise
                error('Sorry, this model is not currently supported.')
        end
        
    elseif carbon == 1 % galactose
        switch lower(model.description)
            case {'iff708_valid_modified.xml'} 
                carbonExchange = {...
                    'GLACxtI'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.ub(carbonExchangeIndex) = 10;

            case {'ind750.xml'} 
                carbonExchange = {...
                    'EX_gal(e)'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'imm904'}
                carbonExchange = {...
                    'EX_gal(e)'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'iaz900.xml'}
                carbonExchange = {...
                    'EX_gal(e)'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'yeast_4.05'}
                carbonExchange = {...
                    'r_1205'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.ub(carbonExchangeIndex) = 10;  % checked - this is right

            case {'yeast_5.01_model.xml'}
                carbonExchange = {...
                    'r_1710'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'yeast_6.06_cobra'}
                carbonExchange = {...
                    'r_1710'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;
            
            case {'imm904_nadcorrected'}
                carbonExchange = {...
                    'EX_gal(e)'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'iin800_cobra', 'iin800.xml', 'ito977_cobra_compatible'}
                carbonExchange = {...
                    'GLACxtI'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.ub(carbonExchangeIndex) = 10;
                
            case {'yeast_7.00_cobra'}
                carbonExchange = {...
                    'r_1710'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;
                
            case {'bmid000000141353'}
                carbonExchange = {...
                    'MNXM390_out'; ... % galactose exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;
                
            otherwise
                error('Sorry, this model is not currently supported.')
        end
        
    elseif carbon == 2 % glycerol
        switch lower(model.description)
            case {'iff708_valid_modified.xml'} 
                carbonExchange = {...
                    'GLxtI'; ... % glycerol exchange'
                };
            carbonExchangeIndex = findRxnIDs(model,carbonExchange);
            model.ub(carbonExchangeIndex) = 10;

            case {'ind750.xml'} 
                carbonExchange = {...
                    'EX_glyc(e)'; ... %glycerol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'imm904'}
                carbonExchange = {...
                    'EX_glyc(e)'; ... %glycerol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'iaz900.xml'}
                carbonExchange = {...
                    'EX_glyc(e)'; ... %glycerol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'yeast_4.05'}
                carbonExchange = {...
                    'r_1299'; ... %glycerol exchange - could also use r_1300
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;  % checked - this is right

            case {'yeast_5.01_model.xml'}
                carbonExchange = {...
                    'r_1808'; ... %glycerol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'yeast_6.06_cobra'}
            % TO CONSIDER: TURN OFF ATP SYNTHASE IN nonfermentable carbon
            % source (ethanol/glycerol)?
                carbonExchange = {...
                    'r_1808'; ... %glycerol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;
                %model.ub(findRxnIDs(model,'r_0226'))=0; % was 1000 - turn off ATP synthase
            
            case {'imm904_nadcorrected'}
                carbonExchange = {...
                    'EX_glyc(e)'; ... %glycerol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'iin800_cobra', 'iin800.xml', 'ito977_cobra_compatible'}
                carbonExchange = {...
                    'GLxtI'; ... % glycerol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.ub(carbonExchangeIndex) = 10;
                
            case {'yeast_7.00_cobra'}
                carbonExchange = {...
                    'r_1808'; ... %glycerol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;
                
            case {'bmid000000141353'}
                error('Sorry, the biomodels model does not have a glycerol exchange reaction.')
                
            otherwise
                error('Sorry, this model is not currently supported.')
        end
        
    elseif carbon == 3 % ethanol
        switch lower(model.description)
            case {'iff708_valid_modified.xml'} 
                carbonExchange = {...
                    'ETHxtI'; ... % ethanol exchange'
                };
            carbonExchangeIndex = findRxnIDs(model,carbonExchange);
            model.ub(carbonExchangeIndex) = 10;

            case {'ind750.xml'} 
                carbonExchange = {...
                    'EX_etoh(e)'; ... % ethanol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'imm904'}
                carbonExchange = {...
                    'EX_etoh(e)'; ... % ethanol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'iaz900.xml'}
                carbonExchange = {...
                    'EX_etoh(e)'; ... % ethanol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'yeast_4.05'}
                carbonExchange = {...
                    'r_1247'; ... % ethanol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.ub(carbonExchangeIndex) = 10; % checked - this is right

            case {'yeast_5.01_model.xml'}
                carbonExchange = {...
                    'r_1761'; ... % ethanol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'yeast_6.06_cobra'}
            % TO CONSIDER: TURN OFF ATP SYNTHASE IN nonfermentable carbon
            % source (ethanol/glycerol)?
                carbonExchange = {...
                    'r_1761'; ... % ethanol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;
                
                %model.ub(findRxnIDs(model,'r_0226'))=0; % was 1000 - turn off ATP synthase

            case {'imm904_nadcorrected'}
                carbonExchange = {...
                    'EX_etoh(e)'; ... % ethanol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;

            case {'iin800_cobra', 'iin800.xml', 'ito977_cobra_compatible'}
                carbonExchange = {...
                    'ETHxtI'; ... % ethanol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.ub(carbonExchangeIndex) = 10;
                
            case {'yeast_7.00_cobra'}
                carbonExchange = {...
                    'r_1761'; ... % ethanol exchange
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;
                
            case {'bmid000000141353'}
                carbonExchange = {...
                    'bigg_etoh_out'; ... % ethanol exchange'
                };
                carbonExchangeIndex = findRxnIDs(model,carbonExchange);
                model.lb(carbonExchangeIndex) = -10;
           
            otherwise
                error('Sorry, this model is not currently supported.')
        end
    end
         
   %% select biomass def
    % all models but the biomodels.db model can grow with the iFF biomass
    % definition, so that's an option
    
    if biomass
        if output
            fprintf('Switching biomass def...\n');
        end
        
        switch lower(model.description)
            case {'iff708_valid_modified.xml'}
                if output
                    fprintf('\tiFF708 biomass already set: no change...\n\n');
                end

            case {'ind750.xml'}
                rxn_index = findRxnIDs(model,'biomass_SC4_bal');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'13BDglcn[c]';'ala-L[c]';'amp[c]';'arg-L[c]';...
                    'asn-L[c]';'asp-L[c]';'atp[c]';'cmp[c]';'cys-L[c]';...
                    'damp[c]';'dcmp[c]';'dgmp[c]';'dtmp[c]';'ergst[c]';...
                    'gln-L[c]';'glu-L[c]';'gly[c]';'glycogen[c]';'gmp[c]';...
                    'h2o[c]';'his-L[c]';'ile-L[c]';'leu-L[c]';'lys-L[c]';...
                    'mannan[c]';'met-L[c]';'pa_SC[c]';'pc_SC[c]';...
                    'pe_SC[c]';'phe-L[c]';'pro-L[c]';'ps_SC[c]';...
                    'ptd1ino_SC[c]';'ser-L[c]';'so4[c]';...
                    'thr-L[c]';'tre[c]';'triglyc_SC[c]';'trp-L[c]';...
                    'tyr-L[c]';'ump[c]';'val-L[c]';'zymst[c]';'adp[c]';...
                    'h[c]';'pi[c]';}, ... %mets to change
                    {'13BDglcn[c]';'ptd1ino_SC[c]';'tre[c]';'amp[c]';...
                    'atp[c]';'cmp[c]';'damp[c]';'dcmp[c]';'dgmp[c]';...
                    'dtmp[c]';'ergst[c]';'gly[c]';'glycogen[c]';'gmp[c]';...
                    'ala-L[c]';'arg-L[c]';'asn-L[c]';'asp-L[c]';'cys-L[c]';...
                    'glu-L[c]';'gln-L[c]';'his-L[c]';'ile-L[c]';'leu-L[c]';...
                    'lys-L[c]';'met-L[c]';'phe-L[c]';'pro-L[c]';'ser-L[c]';...
                    'thr-L[c]';'trp-L[c]';'tyr-L[c]';'val-L[c]';...
                    'mannan[c]';'pa_SC[c]';'pc_SC[c]';'pe_SC[c]';...
                    'ps_SC[c]';'so4[c]';'triglyc_SC[c]';'ump[c]';...
                    'zymst[c]';'adp[c]';'pi[c]';}, ... %new mets
                    {'biomass_SC4_bal';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;-0.006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;59.305;]); %stoichiometric coefficients
                if output
                    fprintf('\tBiomass changed to iFF708 definition...\n\n');
                end

            case {'iin800_cobra'}
                rxn_index = findRxnIDs(model,'CBIOMASS');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'m29';'m241';'m253';'m262';'m316';'m365';'m370';...
                    'm391';'m413';'m483';'m485';'m490';'m554';'m557';...
                    'm559';'m562';'m569';'m575';'m578';'m583';'m588';...
                    'm593';'m597';'m601';'m606';'m609';'m613';'m616';...
                    'm620';'m624';'m627';'m633';'m648';'m864';'m966';...
                    'm230';'m272';'m733';}, ... %mets to change
                    {'m29';'m38';'m241';'m253';'m262';'m316';'m365';...
                    'm370';'m391';'m413';'m429';'m483';'m485';'m490';...
                    'm554';'m557';'m559';'m562';'m569';'m575';'m578';...
                    'm583';'m588';'m593';'m597';'m601';'m606';'m609';...
                    'm613';'m616';'m620';'m624';'m627';'m648';'m757';...
                    'm761';'m762';'m766';'m864';'m928';'m966';'m984';...
                    'm230';'m272';'m733';}, ... %new mets
                    {'CBIOMASS';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;-0.006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;1;59.305;]); %stoichiometric coefficients
                
            case {'iin800.xml'}
                % note: iIN800 doesn't include the species "phosphatidate",
                % so not an exact iFF biomass deff
                rxn_index = findRxnIDs(model,'CBIOMASS');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model = changeRxnMets(model, ...
                    {'1,3-beta-D-glucan [c]';...
                    'alpha,alpha-trehalose [c]';'AMP [c]';'ATP [c]';...
                    'CMP [c]';'dAMP [c]';'dCMP [c]';'dGMP [c]';...
                    'dTMP [c]';'glycine [c]';'glycogen [c]';...
                    'GMP [c]';'L-alanine [c]';'L-arginine [c]';...
                    'L-asparagine [c]';'L-aspartate [c]';...
                    'L-cysteine [c]';'L-glutamate [c]';...
                    'L-glutamine [c]';'L-histidine [c]';...
                    'L-isoleucine [c]';'L-leucine [c]';...
                    'L-lysine [c]';'L-methionine [c]';...
                    'L-phenylalanine [c]';'L-proline [c]';...
                    'L-serine [c]';'L-threonine [c]';...
                    'L-tryptophan [c]';'L-tyrosine [c]';...
                    'L-valine [c]';'lipids [c]';'mannan [c]';...
                    'sulfate [c]';'UMP [c]';'ADP [c]';'biomass [c]';...
                    'phosphate [c]'}, ... %mets to change
                    {'1,3-beta-D-glucan [c]';...
                    '1-phosphatidyl-1D-myo-inositols [c]';...
                    'alpha,alpha-trehalose [c]';'AMP [c]';'ATP [c]';...
                    'CMP [c]';'dAMP [c]';'dCMP [c]';'dGMP [c]';...
                    'dTMP [c]';'ergosterol [c]';'glycine [c]';...
                    'glycogen [c]';'GMP [c]';'L-alanine [c]';...
                    'L-arginine [c]';'L-asparagine [c]';...
                    'L-aspartate [c]';'L-cysteine [c]';...
                    'L-glutamate [c]';'L-glutamine [c]';...
                    'L-histidine [c]';'L-isoleucine [c]';...
                    'L-leucine [c]';'L-lysine [c]';'L-methionine [c]';...
                    'L-phenylalanine [c]';'L-proline [c]';...
                    'L-serine [c]';'L-threonine [c]';...
                    'L-tryptophan [c]';'L-tyrosine [c]';...
                    'L-valine [c]';'mannan [c]';...
                    'phosphatidylcholines [c]';...
                    'phosphatidylethanolamines [c]';...
                    'phosphatidylserine [c]';'sulfate [c]';...
                    'triglycerides [c]';'UMP [c]';'zymosterol [c]';...
                    'ADP [c]';'biomass [c]';'phosphate [c]'}, ... %new mets
                    {'CBIOMASS';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;1;59.305;]); %stoichiometric coefficients
                
                if output
                    fprintf('\tBiomass changed to iFF708 definition...\n\n');
                end

            case {'imm904'}
                rxn_index = findRxnIDs(model,'biomass_SC5_notrace');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'13BDglcn[c]';'ala-L[c]';'amp[c]';'arg-L[c]';...
                    'asn-L[c]';'asp-L[c]';'atp[c]';'cmp[c]';'cys-L[c]';...
                    'damp[c]';'dcmp[c]';'dgmp[c]';'dtmp[c]';'ergst[c]';...
                    'gln-L[c]';'glu-L[c]';'gly[c]';'glycogen[c]';'gmp[c]';...
                    'h2o[c]';'his-L[c]';'ile-L[c]';'leu-L[c]';'lys-L[c]';...
                    'mannan[c]';'met-L[c]';'pa_SC[c]';'pc_SC[c]';...
                    'pe_SC[c]';'phe-L[c]';'pro-L[c]';'ps_SC[c]';...
                    'ptd1ino_SC[c]';'ribflv[c]';'ser-L[c]';'so4[c]';...
                    'thr-L[c]';'tre[c]';'triglyc_SC[c]';'trp-L[c]';...
                    'tyr-L[c]';'ump[c]';'val-L[c]';'zymst[c]';'adp[c]';...
                    'h[c]';'pi[c]';}, ... %mets to change
                    {'13BDglcn[c]';'ptd1ino_SC[c]';'tre[c]';'amp[c]';...
                    'atp[c]';'cmp[c]';'damp[c]';'dcmp[c]';'dgmp[c]';...
                    'dtmp[c]';'ergst[c]';'gly[c]';'glycogen[c]';'gmp[c]';...
                    'ala-L[c]';'arg-L[c]';'asn-L[c]';'asp-L[c]';'cys-L[c]';...
                    'glu-L[c]';'gln-L[c]';'his-L[c]';'ile-L[c]';'leu-L[c]';...
                    'lys-L[c]';'met-L[c]';'phe-L[c]';'pro-L[c]';'ser-L[c]';...
                    'thr-L[c]';'trp-L[c]';'tyr-L[c]';'val-L[c]';...
                    'mannan[c]';'pa_SC[c]';'pc_SC[c]';'pe_SC[c]';...
                    'ps_SC[c]';'so4[c]';'triglyc_SC[c]';'ump[c]';...
                    'zymst[c]';'adp[c]';'pi[c]';}, ... %new mets 
                    {'biomass_SC5_notrace';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;-0.006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;59.305;]); %stoichiometric coefficients
                if output
                    fprintf('\tBiomass changed to iFF708 definition...\n\n');
                end

            case {'iaz900.xml'}
                rxn_index = findRxnIDs(model,'biomass_wild');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'13BDglcn[c]';'16BDglcn[c]';'alatrna[c]';'amp[c]';...
                    'argtrna[c]';'asntrna[c]';'asptrna[c]';'atp[c]';...
                    'cmp[c]';'cystrna[c]';'damp[c]';'dcmp[c]';'dgmp[c]';...
                    'dtmp[c]';'glntrna[c]';'glutrna[c]';'glycogen[c]';...
                    'glytrna[c]';'gmp[c]';'h2o[c]';'histrna[c]';...
                    'iletrna[c]';'leutrna[c]';'lystrna[c]';'mannan[c]';...
                    'mettrna[c]';'pa_SC[c]';'pe_SC[c]';'phetrna[c]';...
                    'phospholipid[c]';'protrna[c]';'ps_SC[c]';...
                    'ptd1ino_SC[c]';'ribflv[c]';'sertrna[c]';'so4[c]';...
                    'sterol[c]';'thrtrna[c]';'tre[c]';'triglyc_SC[c]';...
                    'trptrna[c]';'tyrtrna[c]';'ump[c]';'valtrna[c]';...
                    'zymst[c]';'adp[c]';'h[c]';'pi[c]';'trnaala[c]';...
                    'trnaarg[c]';'trnaasn[c]';'trnaasp[c]';'trnacys[c]';...
                    'trnagln[c]';'trnaglu[c]';'trnagly[c]';'trnahis[c]';...
                    'trnaile[c]';'trnaleu[c]';'trnalys[c]';'trnamet[c]';...
                    'trnaphe[c]';'trnapro[c]';'trnaser[c]';'trnathr[c]';...
                    'trnatrp[c]';'trnatyr[c]';'trnaval[c]';}, ... %mets to change
                    {'13BDglcn[c]';'ptd1ino_SC[c]';'tre[c]';'amp[c]';...
                    'atp[c]';'cmp[c]';'damp[c]';'dcmp[c]';'dgmp[c]';...
                    'dtmp[c]';'ergst[c]';'gly[c]';'glycogen[c]';'gmp[c]';...
                    'ala_L[c]';'arg_L[c]';'asn_L[c]';'asp_L[c]';'cys_L[c]';...
                    'glu_L[c]';'gln_L[c]';'his_L[c]';'ile_L[c]';'leu_L[c]';...
                    'lys_L[c]';'met_L[c]';'phe_L[c]';'pro_L[c]';'ser_L[c]';...
                    'thr_L[c]';'trp_L[c]';'tyr_L[c]';'val_L[c]';...
                    'mannan[c]';'pa_SC[c]';'pc_SC[c]';'pe_SC[c]';...
                    'ps_SC[c]';'so4[c]';'triglyc_SC[c]';'ump[c]';...
                    'zymst[c]';'adp[c]';'pi[c]';}, ... %new mets
                    {'biomass_wild';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;-0.006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;59.305;]); %stoichiometric coefficients
                if output
                    fprintf('\tBiomass changed to iFF708 definition...\n\n');
                end

            case {'yeast_4.05'}
                rxn_index = findRxnIDs(model,'r_1812');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'s_0002';'s_0416';'s_0434';'s_0446';'s_0511';...
                    's_0564';'s_0569';'s_0593';'s_0619';'s_0740';'s_0743';...
                    's_0752';'s_0863';'s_0873';'s_0877';'s_0881';'s_0889';...
                    's_0899';'s_0907';'s_0911';'s_0920';'s_0925';'s_0929';...
                    's_0933';'s_0936';'s_0939';'s_0943';'s_0949';'s_0952';...
                    's_0955';'s_0960';'s_1000';'s_1011';'s_1283';'s_1347';...
                    's_1417';'s_0400';'s_0463';'s_1207';}, ... %mets to change
                    {'s_0002';'s_0091';'s_0416';'s_0434';'s_0446';...
                    's_0511';'s_0564';'s_0569';'s_0593';'s_0619';'s_0636';...
                    's_0740';'s_0743';'s_0752';'s_0863';'s_0873';'s_0877';...
                    's_0881';'s_0889';'s_0899';'s_0907';'s_0911';'s_0920';...
                    's_0925';'s_0929';'s_0933';'s_0936';'s_0939';'s_0943';...
                    's_0949';'s_0952';'s_0955';'s_0960';'s_1011';'s_1216';...
                    's_1229';'s_1234';'s_1220';'s_1347';'s_1399';'s_1417';...
                    's_1448';'s_0400';'s_0463';'s_1207';}, ... %new mets
                    {'r_1812';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;-0.006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;1;59.305;]); %stoichiometric coefficients
                if output
                    fprintf('\tBiomass changed to iFF708 definition...\n\n');
                end

            case {'yeast_5.01_model.xml'}
                rxn_index = findRxnIDs(model,'r_2110');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'s_0002';'s_0423';'s_0434';'s_0526';'s_0584';...
                    's_0589';'s_0615';'s_0649';'s_0773';'s_0782';'s_0803';...
                    's_0955';'s_0965';'s_0969';'s_0973';'s_0981';'s_0991';...
                    's_0999';'s_1003';'s_1006';'s_1016';'s_1021';'s_1025';...
                    's_1029';'s_1032';'s_1035';'s_1039';'s_1045';'s_1048';...
                    's_1051';'s_1056';'s_1096';'s_1107';'s_1405';'s_1467';...
                    's_1520';'s_1545';'s_0394';'s_0450';'s_0794';...
                    's_1322';}, ... %mets to change
                    {'s_0002';'s_0089';'s_1520';'s_0423';'s_0434';...
                    's_0526';'s_0584';'s_0589';'s_0615';'s_0649';'s_0666';...
                    's_1003';'s_0773';'s_0782';'s_0955';'s_0965';'s_0969';...
                    's_0973';'s_0981';'s_0991';'s_0999';'s_1006';'s_1016';...
                    's_1021';'s_1025';'s_1029';'s_1032';'s_1035';'s_1039';...
                    's_1045';'s_1048';'s_1051';'s_1056';'s_1107';'s_1331';...
                    's_1346';'s_1351';'s_1337';'s_1467';'s_1524';'s_1545';...
                    's_1569';'s_0394';'s_0450';'s_1322';}, ... %new mets
                    {'r_2110';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;-0.006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;1;59.305;]); %stoichiometric coefficients
                if output
                    fprintf('\tBiomass changed to iFF708 definition...\n\n');
                end

            case {'yeast_6.06_cobra'}
                % clear the current biomass def
                rxn_index = findRxnIDs(model,'r_2133');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'s_0002';'s_0423';'s_0434';'s_0526';'s_0584';...
                    's_0589';'s_0615';'s_0649';'s_0773';'s_0782';'s_0803';...
                    's_0955';'s_0965';'s_0969';'s_0973';'s_0981';'s_0991';...
                    's_0999';'s_1003';'s_1006';'s_1016';'s_1021';'s_1025';...
                    's_1029';'s_1032';'s_1035';'s_1039';'s_1045';'s_1048';...
                    's_1051';'s_1056';'s_1096';'s_1107';'s_1405';'s_1467';...
                    's_1520';'s_1545';'s_0394';'s_0450';'s_0794';...
                    's_1322';}, ... %mets to change
                    {'s_0002';'s_0089';'s_1520';'s_0423';'s_0434';...
                    's_0526';'s_0584';'s_0589';'s_0615';'s_0649';'s_0666';...
                    's_1003';'s_0773';'s_0782';'s_0955';'s_0965';'s_0969';...
                    's_0973';'s_0981';'s_0991';'s_0999';'s_1006';'s_1016';...
                    's_1021';'s_1025';'s_1029';'s_1032';'s_1035';'s_1039';...
                    's_1045';'s_1048';'s_1051';'s_1056';'s_1107';'s_1331';...
                    's_1346';'s_1351';'s_1337';'s_1467';'s_1524';'s_1545';...
                    's_1569';'s_0394';'s_0450';'s_1322';}, ... %new mets
                    {'r_2133';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;-0.006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;1;59.305;]); %stoichiometric coefficients
                if output
                    fprintf('\tBiomass changed to iFF708 definition...\n\n');
                end

            case {'imm904_nadcorrected'}
                rxn_index = findRxnIDs(model,'biomass_SC5_notrace');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'13BDglcn[c]';'ala_L[c]';'amp[c]';'arg_L[c]';...
                    'asn_L[c]';'asp_L[c]';'atp[c]';'camp[c]';...
                    'chitin[c]';'cmp[c]';'coa[c]';'cys_L[c]';...
                    'damp[c]';'dcmp[c]';'dgmp[c]';'dtmp[c]';'ergst[c]';...
                    'fad[c]';'gln_L[c]';'glu_L[c]';'gly[c]';...
                    'glycogen[c]';'gmp[c]';'gthrd[c]';'h2o[c]';...
                    'his_L[c]';'ile_L[c]';'leu_L[c]';'lys_L[c]';...
                    'mannan[c]';'met_L[c]';'nad[c]';'pa_SC[c]';...
                    'pc_SC[c]';'pe_SC[c]';'phe_L[c]';'pheme[m]';...
                    'pro_L[c]';'ps_SC[c]';'ptd1ino_SC[c]';'q6[m]';...
                    'ribflv[c]';'ser_L[c]';'so4[c]';'thf[c]';...
                    'thmtp[c]';'thr_L[c]';'tre[c]';'triglyc_SC[c]';'trp_L[c]';...
                    'tyr_L[c]';'ump[c]';'val_L[c]';'zymst[c]';'adp[c]';...
                    'h[c]';'pi[c]';}, ... %mets to change
                    {'13BDglcn[c]';'ptd1ino_SC[c]';'tre[c]';'amp[c]';...
                    'atp[c]';'cmp[c]';'damp[c]';'dcmp[c]';'dgmp[c]';...
                    'dtmp[c]';'ergst[c]';'gly[c]';'glycogen[c]';'gmp[c]';...
                    'ala_L[c]';'arg_L[c]';'asn_L[c]';'asp_L[c]';'cys_L[c]';...
                    'glu_L[c]';'gln_L[c]';'his_L[c]';'ile_L[c]';'leu_L[c]';...
                    'lys_L[c]';'met_L[c]';'phe_L[c]';'pro_L[c]';'ser_L[c]';...
                    'thr_L[c]';'trp_L[c]';'tyr_L[c]';'val_L[c]';...
                    'mannan[c]';'pa_SC[c]';'pc_SC[c]';'pe_SC[c]';...
                    'ps_SC[c]';'so4[c]';'triglyc_SC[c]';'ump[c]';...
                    'zymst[c]';'adp[c]';'pi[c]';}, ... %new mets 
                    {'biomass_published';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;-0.006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;59.305;]); %stoichiometric coefficients
                if output
                    fprintf('\tBiomass changed to iFF708 definition...\n\n');
                end
                
            case {'ito977_cobra_compatible'}
                rxn_index = findRxnIDs(model,'CBIOMASS');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'m160';'m584';'m339';'m340';'m358';'m618';'m621';...
                    'm624';'m627';'m682';'m683';'m415';'m440';'m442';...
                    'm444';'m447';'m451';'m457';'m458';'m461';'m466';...
                    'm469';'m471';'m473';'m476';'m478';'m481';'m483';...
                    'm485';'m487';'m489';'m732';'m816';'m562';'m724';...
                    'm335';'m600';'m768';}, ... %mets to change
                    {'m160';'m174';'m584';'m339';'m340';'m358';'m618';...
                    'm621';'m624';'m627';'m662';'m682';'m683';'m415';...
                    'm440';'m442';'m444';'m447';'m451';'m457';'m458';...
                    'm461';'m466';'m469';'m471';'m473';'m476';'m478';...
                    'm481';'m483';'m485';'m487';'m489';'m732';'m769';...
                    'm772';'m773';'m774';'m816';'m884';'m562';'m894';...
                    'm335';'m600';'m768';}, ... %new mets
                    {'CBIOMASS';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.0460;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.0460;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.1020;-0.2646;-0.8079;-0.0006;-0.0060;...
                    -0.0045;-0.0017;-0.0200;-0.0066;-0.0599;-0.0015;...
                    59.276;1;59.3050;]); %stoichiometric coefficients
                
            case {'yeast_7.00_cobra'}
                %note: Y7 includes the Y6 biomass def (rxn r_2133, with
                %lb=0, ub=1000), and the Y5 biomass def (r_2110, with lb=0,
                %ub=0). Both include the same lipid species, s_1096, made
                %in reaction r_2108. Below, I modify it as for Y6.

                % clear the current biomass def
                rxn_index = findRxnIDs(model,'r_2133');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'s_0002';'s_0423';'s_0434';'s_0526';'s_0584';...
                    's_0589';'s_0615';'s_0649';'s_0773';'s_0782';'s_0803';...
                    's_0955';'s_0965';'s_0969';'s_0973';'s_0981';'s_0991';...
                    's_0999';'s_1003';'s_1006';'s_1016';'s_1021';'s_1025';...
                    's_1029';'s_1032';'s_1035';'s_1039';'s_1045';'s_1048';...
                    's_1051';'s_1056';'s_1096';'s_1107';'s_1405';'s_1467';...
                    's_1520';'s_1545';'s_0394';'s_0450';'s_0794';...
                    's_1322';}, ... %mets to change
                    {'s_0002';'s_0089';'s_1520';'s_0423';'s_0434';...
                    's_0526';'s_0584';'s_0589';'s_0615';'s_0649';'s_0666';...
                    's_1003';'s_0773';'s_0782';'s_0955';'s_0965';'s_0969';...
                    's_0973';'s_0981';'s_0991';'s_0999';'s_1006';'s_1016';...
                    's_1021';'s_1025';'s_1029';'s_1032';'s_1035';'s_1039';...
                    's_1045';'s_1048';'s_1051';'s_1056';'s_1107'; ...
                    's_2954';... %-Y7 doesn't have phosphatidate [cytoplasm] (met 1190, s_1331), use s_2954 instead
                    's_1346';'s_1351';'s_1337';'s_1467';'s_1524';'s_1545';...
                    's_1569';'s_0394';'s_0450';'s_1322';}, ... %new mets 
                    {'r_2133';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;-0.006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;1;59.305;]); %stoichiometric coefficients
                
                if output
                    fprintf('\tBiomass changed to iFF708 definition...\n\n');
                end

            case {'bmid000000141353'}
                rxn_index = findRxnIDs(model,'BIOMASS_REACTION');
                model.S(find(model.S(:,rxn_index)),rxn_index) = 0;
                model=changeRxnMets(model,...
                    {'bigg_gly_bm'; 'bigg_tyr_L_bm'; 'bigg_asp_L_bm'; ...
                    'bigg_lys_L_bm'; 'bigg_thr_L_bm'; 'bigg_ile_L_bm'; ...
                    'bigg_leu_L_bm'; 'bigg_val_L_bm'; 'bigg_arg_L_bm'; ...
                    'bigg_asn_L_bm'; 'MNXM18_bm'; 'bigg_gln_L_bm'; ...
                    'bigg_pro_L_bm'; 'bigg_his_L_bm'; 'bigg_ala_L_bm'; ...
                    'bigg_ser_L_bm'; 'bigg_cys_L_bm'; 'bigg_met_L_bm'; ...
                    'bigg_phe_L_bm'; 'bigg_trp_L_bm'; 'bigg_amp_bm'; ...
                    'bigg_gmp_bm'; 'bigg_cmp_bm'; 'bigg_ump_bm'; ...
                    'bigg_damp_bm'; 'bigg_dgmp_bm'; 'bigg_dcmp_bm'; ...
                    'bigg_dtmp_bm'; 'MNXM876_bm'; 'bigg_atp_bm'; ...
                    'bigg_adp_bm'; 'bigg_pi_bm';}, ... %mets to change
                    {'MNXM2056[i]'; 'MNXM62[i]'; 'bigg_tre[i]'; ...
                    'bigg_amp_bm'; 'bigg_atp_bm'; 'bigg_cmp_bm'; ...
                    'bigg_damp_bm'; 'bigg_dcmp_bm'; 'bigg_dgmp_bm'; ...
                    'bigg_dtmp_bm'; 'bigg_ergst[i]'; 'bigg_gly_bm'; ...
                    'MNXM876_bm'; 'bigg_gmp_bm'; 'bigg_ala_L_bm'; ...
                    'bigg_arg_L_bm'; 'bigg_asn_L_bm'; 'bigg_asp_L_bm'; ...
                    'bigg_cys_L_bm'; 'MNXM18_bm'; 'bigg_gln_L_bm'; ...
                    'bigg_his_L_bm'; 'bigg_ile_L_bm'; 'bigg_leu_L_bm'; ...
                    'bigg_lys_L_bm'; 'bigg_met_L_bm'; 'bigg_phe_L_bm'; ...
                    'bigg_pro_L_bm'; 'bigg_ser_L_bm'; 'bigg_thr_L_bm'; ...
                    'bigg_trp_L_bm'; 'bigg_tyr_L_bm'; 'bigg_val_L_bm'; ...
                    'MNXM30402[i]'; 'MNXM17312[i]'; 'MNXM96041[i]'; ...
                    'bigg_ethamp[i]'; 'MNXM36806[i]'; 'bigg_so4[i]'; ...
                    'MNXM248[i]'; 'bigg_ump_bm'; 'bigg_zymst[i]'; ...
                    'bigg_adp_bm'; 'bigg_pi_bm';}, ... %new mets
                    {'BIOMASS_REACTION';}, ... %rxn id
                    [-1.1348;-0.0053;-0.0234;-0.046;-59.276;-0.0447;...
                    -0.0036;-0.0024;-0.0024;-0.0036;-0.0007;-0.2904;...
                    -0.5185;-0.046;-0.4588;-0.1607;-0.1017;-0.2975;...
                    -0.0066;-0.3018;-0.1054;-0.0663;-0.1927;-0.2964;...
                    -0.2862;-0.0507;-0.1339;-0.1647;-0.1854;-0.1914;...
                    -0.0284;-0.102;-0.2646;-0.8079;-0.0006;-0.006;...
                    -0.0045;-0.0017;-0.02;-0.0066;-0.0599;-0.0015;...
                    59.276;59.305;]); %stoichiometric coefficients
                
                if output
                    fprintf('\tBiomass changed to iFF708 definition...\n\n');
                end
                
            otherwise
                error('Sorry, this model is not currently supported.')
        end
    else
        if output
            fprintf('Using model-default biomass def...\n');
        end
    end
        
    %% describe medium
    if output
        fprintf('\nCurrent growth medium composition:\n');
    end
    
    switch lower(model.description)
        case {'iff708_valid_modified.xml'} 
            uptake_rxn_ids = {'ZYMSTxtI';'XANxtI';'URIxtI';'UREAxtI';...
                'URAxtI';'THYxtI';'DTxtI';'THMxtI';'SLFxtI';'SUCxtI';...
                'SUCCxtI';'C180xtI';'SPRMxtI';'SORxtI';'NAxtI';...
                'MMETxtI';'SAMxtI';'RFLAVxtI';'PYRxtI';'KxtI';...
                'PIMExtI';'PIxtI';'PEPTxtI';'C160xtI';'O2xtI';'OGTxtI';...
                'OPEPxtI';'NMNxtI';'NH3xtI';'C140xtI';'MIxtI';...
                'MTHNxtI';'MELIxtI';'MLTxtI';'MALxtI';'VALxtI';...
                'TYRxtI';'TRPxtI';'THRxtI';'SERxtI';'PROxtI';...
                'PHExtI';'ORNxtI';'METxtI';'LYSxtI';'LEUxtI';...
                'ILExtI';'HISxtI';'GLNxtI';'GLUxtI';'GLTxtI';'CYSxtI';...
                'ASPxtI';'ASNxtI';'ARGxtI';'ALAxtI';'INSxtI';'HYXNxtI';...
                'GSNxtI';'GNxtI';'GLALxtI';'GLYxtI';'GLxtI';'GLAMxtI';...
                'GABAxtI';'FUMxtI';'FORxtI';'ETHxtI';'ERGOxtI';...
                'DTTPxtI';'RIBxtI';'RMNxtI';'MNTxtI';'DIPEPxtI';...
                'GA6PxtI';'GLACxtI';'FRUxtI';'DUxtI';'DINxtI';'DGxtI';...
                'DCxtI';'DAxtI';'ARABxtI';'CYTSxtI';'CYTDxtI';...
                'CO2xtI';'CITxtI';'CHOxtI';'BTxtI';'BMxtI';'FUCxtI';...
                'MANxtI';'GLCxtI';'TRExtI';'ATNxtI';'PAPxtI';'ADNxtI';...
                'ADxtI';'ACxtI';'ACALxtI';'AONAxtI';'DANNAxtI';...
                'AKGxtI';'PNTOxtI';'LACxtI';};
            uptake_rxn_indexes = findRxnIDs(model, uptake_rxn_ids);

            media_indexes = uptake_rxn_indexes(find(model.ub(uptake_rxn_indexes)));

            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);

        case {'ind750.xml'}
            exchangeRxns = findExcRxns(model);
            media_indexes = intersect(find(exchangeRxns), find(model.lb ~=0));
            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);

        case {'iin800_cobra'}
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
            'XANxtI';'XTSINExtI';'XYLxtI';'ZYMSTxtI';...
            'Ex_m758';'Ex_m928'; ... % the TAG and phosphatidate reactions added
            };
        
            uptake_rxn_indexes = findRxnIDs(model, uptake_rxn_ids);

            media_indexes = uptake_rxn_indexes(find(model.ub(uptake_rxn_indexes)));
            media_indexes = [media_indexes; ...
                uptake_rxn_indexes(find(model.lb(uptake_rxn_indexes)))];

            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);

        case {'iin800.xml'}
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
        
            uptake_rxn_indexes = findRxnIDs(model, uptake_rxn_ids);

            media_indexes = uptake_rxn_indexes(find(model.ub(uptake_rxn_indexes)));
            media_indexes = [media_indexes; ...
                uptake_rxn_indexes(find(model.lb(uptake_rxn_indexes)))];

            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);
        
        case {'imm904'}
            exchangeRxns = findExcRxns(model);
            media_indexes = intersect(find(exchangeRxns), find(model.lb ~=0));
            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);

        case {'iaz900.xml'}
            exchangeRxns = findExcRxns(model);
            media_indexes = intersect(find(exchangeRxns), find(model.lb ~=0));
            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);

        case {'yeast_4.05'}
            %first, refind exchange reactions in case we're not doing minimal
            %media
            exchangeRxns = findExcRxns(model);
            exchange_indexes=find(exchangeRxns);

            %thanks to Shuyi Ma (of the Price lab at ISB) for this!
            b2 = sum(model.S(:,exchange_indexes),1); %-1 for exchange as reactant +1 for exchange as product
            posexch = exchange_indexes(b2 == 1); % index of exchange rxns having product
            negexch = exchange_indexes(b2 == -1); % index of exchange rxns having substrate
            
            source = intersect(posexch, find(model.ub > 0));
            source = [source; intersect(negexch, find(model.lb < 0))];
                
            % sink = intersect(negexch, find(model.ub > 0));
            % sink = [sink; intersect(negexch, find(model.lb < 0))];
                   
            disp([model.rxns(source) model.rxnNames(source) ...
            num2cell(model.lb(source)) num2cell(model.ub(source))]);

        case {'yeast_5.01_model.xml'}
            exchangeRxns = findExcRxns(model);
            media_indexes = intersect(find(exchangeRxns), find(model.lb ~=0));
            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);

        case {'yeast_6.06_cobra'}
            exchangeRxns = findExcRxns(model);
            media_indexes = intersect(find(exchangeRxns), find(model.lb ~=0));
            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);
 
        case {'yeast_7.00_cobra'}
            exchangeRxns = findExcRxns(model);
            media_indexes = intersect(find(exchangeRxns), find(model.lb ~=0));
            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);
 
        case {'imm904_nadcorrected'}
            exchangeRxns = findExcRxns(model);
            media_indexes = intersect(find(exchangeRxns), find(model.lb ~=0));
            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);
        
        case {'ito977_cobra_compatible'}
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
    
            uptake_rxn_indexes = findRxnIDs(model, uptake_rxn_ids);

            media_indexes = uptake_rxn_indexes(find(model.ub(uptake_rxn_indexes)));
            media_indexes = [media_indexes; ...
                uptake_rxn_indexes(find(model.lb(uptake_rxn_indexes)))];

            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);
        
       case {'bmid000000141353'} 
            uptakeRxns = findRxnIDs(model, ...
                {'MNXR6704_i'; 'MNXM99_in'; 'MNXM105_in'; ...
                'MNXM15_in'; 'MNXM27_in'; 'MNXM95_in'; 'MNXM653_in'; ...
                'MNXM128_in'; 'MNXM58_in'; 'MNXM43_in'; 'MNXM9_in'; ...
                'MNXM1_in'; 'MNXM2_in'; 'MNXM13_in'; 'MNXM4_in';});

            otherRxns = setdiff(find(findExcRxns(model)), uptakeRxns);
                        
            media_indexes = [uptakeRxns(model.ub(uptakeRxns)>0); ...
                otherRxns(model.lb(otherRxns)<0)];
            
            disp([model.rxns(media_indexes) model.rxnNames(media_indexes) ...
            num2cell(model.lb(media_indexes)) num2cell(model.ub(media_indexes))]);

        
        otherwise
            error('Sorry, this model is not currently supported.')
    end
        
    %% add auxotrophic markers
    
    if output
        fprintf('\nAdding strain auxotrophic markers...\n');
    end
    
    KO_genes = {...
        'YOR202W', ... % his3d1
        'YCL018W', ... % leu2d0
        'YLR303W', ... % met15d0 - now called met17
        'YEL021W'... % ura3d0
        };
    
    [model,~,~,~] = deleteModelGenes(model, KO_genes);
    
    %% prepare KO analysis
    
    if output
        fprintf('\nKnockout analysis...\n');
    end
    
    analyzed = intersect(model.genes, kuepfer_all_genes);
    dubious = setdiff(intersect(analyzed, verifiedORFs), analyzed);
    
    if output
        fprintf('\tNumber of genes analzyed: %i\n', length(analyzed));
    end
    
    if output
        fprintf('\tDubious ORFs analzyed: %i\n', length(dubious));
    end
    
    %% knockout analysis
    
    if carbon == 0
        inviableORFsAll = kuepfer_glucose_aerobic;
    elseif carbon == 1
        inviableORFsAll = kuepfer_galactose_aerobic;
    elseif carbon == 2
        inviableORFsAll = kuepfer_glycerol_aerobic;
    elseif carbon == 3
        inviableORFsAll = kuepfer_ethanol_aerobic;
    end
        
    % restrict essentlial gene list to those analyzed by Kuepfer et al.
    exp_inviable = intersect(analyzed,inviableORFsAll);
    exp_inviable = intersect(exp_inviable,verifiedORFs);

    exp_viable = setdiff(analyzed,inviableORFsAll);
    % exp_viable = intersect(exp_viable,verifiedORFs);

    grRatio = singleGeneDeletion(model, '', analyzed);

    % mod_viable  = model.genes(grRatio >= ko_tol);% WRONG b/c only doing
    % a subset
    mod_viable  = analyzed(grRatio >= ko_tol);
    % mod_viable = intersect(mod_viable,verifiedORFs);
    % mod_inviable = model.genes(grRatio < ko_tol); % WRONG b/c only doing
    % a subset
    mod_inviable = analyzed(grRatio < ko_tol);
    % mod_inviable = intersect(mod_inviable,verifiedORFs);

    tp = intersect(exp_viable,mod_viable); n_tp = length(tp);
    tn = intersect(exp_inviable,mod_inviable); n_tn = length(tn);
    fp = intersect(exp_inviable,mod_viable); n_fp = length(fp);
    fn = intersect(exp_viable,mod_inviable); n_fn = length(fn);

    results.tp = tp;
    results.tn = tn;
    results.fp = fp;
    results.fn = fn;

    if output
        fprintf('\nknockout analysis (positive = viable):');
        fprintf('\n\ttp: %i\n',n_tp);
        fprintf('\ttn: %i\n',n_tn);
        fprintf('\tfp: %i\n',n_fp);
        fprintf('\tfn: %i\n\n',n_fn);
    end

    % n_genes = length(intersect(model.genes,verifiedORFs));

    sensitivity = (100*n_tp/(n_tp+n_fn));
    specificity = (100*n_tn/(n_tn+n_fp));
    positivePredictive = (100*n_tp/(n_tp+n_fp));
    negativePredictive = (100*n_tn/(n_fn+n_tn));
    mcc = (n_tp * n_tn - n_fp * n_fn)/...
        sqrt((n_tp + n_fp)*(n_tp + n_fn)*(n_tn + n_fp)*(n_tn + n_fn));
    geoMean = (sensitivity * specificity)^.5;

    if output
        fprintf('\t%.2f%%\t%s\n\t%.2f%%\t%s\n\t%.2f%%\t%s\n\t%.2f%%\t%s\n\n\t%.2f%%\t%s\t%s\n\t%.2f\t%s\n\t\t\t%s\n\n',...
            sensitivity,'sensitivity = recall = tp/(tp+fn)',...
            specificity,'specificity = tn/(tn+fp)',...
            positivePredictive,'positive predictive value = precision = tp/(tp+fp)',...
            negativePredictive,'negative predictive value = tn/(fn+tn)',...
            geoMean, 'Geometric mean accuracy', ...
            '(sensitivity * specificity)^.5', ...
            mcc, 'Matthews correlation coefficient',...
            '(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn_fp)*(tn+fn))');
    end
    
    if output == 2
        fprintf('\tlist of ko false positives:\n');
        for k = 1:n_fp
            fprintf('\t\t%s\n',fp{k});
        end

        fprintf('\n\tlist of ko false negatives:\n');
        for k = 1:n_fn
            fprintf('\t\t%s\n',fn{k});
        end
        disp(' ');
    end

end

%% required functions

    % kuepfer media: "minimal medium was, per liter (Verduyn et al.
    % 1992),5g of (NH4)2SO4, 3g of KH2PO4, 0.5 g of MgSO4�7H2O, 4.5 mg of
    % ZnSO4�7H2O, 0.3 mg of CoCl2�6H2O, 1.0 mg of MnCl2�4H2O, 0.3 mg of
    % CuSO4�5H2O, 4.5 mg of CaCl2�2H2O, 3.0 mg of FeSO4�7H2O, 0.4 mg of
    % NaMoO4�2H2O, 1.0 mg of H3BO3, 0.1 mg of KI, 15 mg of EDTA, 0.05 mg of
    % biotin, 1.0 mg of Ca-pantothenate, 1.0 mg of nicotinic acid, 25 mg of
    % inositol, 1.0 mg of pyridoxine, 0.2 mg of p-amino-benzoic acid, and
    % 1.0 mg of thiamine. The carbon sources(ethanol, galactose, glucose,
    % and glycerol) were added to a final concentration of 20 g/L. Strain
    % auxotrophies were complemented with 20 mg/L histidine, uracil,
    % methionine, and 60 mg/L leucine. About 50 strains of the yeast
    % collection are lysine auxotroph and were independently tested for
    % growth on plates supplemented with 20 mg/L lysine."
    %
    % I will assume that this media is carbon-limited, so allow
    % unconstrained uptake of all media components I can identify in each
    % model.
    
function model = kuepfer_iFF(model)
    % change iFF model media to kuepfer
    
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
        'DANNAxtO';'AKGxtO';'PNTOxtO';'LACxtO';};

    uptake_rxn_indexes = findRxnIDs(model, uptake_rxn_ids);
    excretion_rxn_indexes = findRxnIDs(model, ...
        excretion_rxn_ids);

    model.lb(uptake_rxn_indexes) = 0;
    model.ub(uptake_rxn_indexes) = 0;
    model.lb(excretion_rxn_indexes) = 0;
    model.ub(excretion_rxn_indexes) = 1000;

    desiredExchanges = {...
        'BTxtI'; ... % biotin
        'PNTOxtI'; ... % Pantothenate
        'MIxtI'; ... % inositol
        'THMxtI'; ... % thiamine
        'HISxtI'; ... % histidine
        'URAxtI'; ... % uracil
        'METxtI'; ... % methionine
        'LEUxtI'; ... % leucine
        'SLFxtI'; ... % sufate
        'PIxtI'; ... % phosphate
        'NAxtI'; ... % sodium
        'KxtI'; ... % potassium
        'O2xtI'; ... % oxygen
        'NH3xtI'; ... % ammonia
        'NMNxtI'; ... % nicotinamide mononucleotide
        };

    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);

    if length(uptakeRxnIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.ub(uptakeRxnIndexes) = 1000;
end

function model = kuepfer_iND(model)
    % change iND model media to kuepfer
    
    % start with a clean slate: set to unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'EX_btn(e)'; ... % biotin
        'EX_pnto_R(e)'; ... % pantothenate
        'EX_inost(e)'; ... % inositol
        'EX_thm(e)'; ... % thiamine
        'EX_his_L(e)'; ... % histidine
        'EX_ura(e)'; ... % uracil
        'EX_met_L(e)'; ... % methionine
        'EX_leu_L(e)'; ... % leucine
        'EX_so4(e)'; ... % sulfate
        'EX_pi(e)'; ... % phosphate
        'EX_na1(e)'; ... % sodium
        'EX_k(e)'; ... % potassium
        'EX_o2(e)'; ... % oxygen
        'EX_h2o(e)'; ... % water
        'EX_h(e)'; ... % protons
        'EX_nh4(e)'; ... % ammonium
        'EX_thmpp(e)'; ... % thiamine diphosphate
        };

    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    if length(uptakeRxnIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
end

function model = kuepfer_iMM(model)
    % change iMM model media to kuepfer
    
    % start with a clean slate: set to unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'EX_btn(e)'; ... % biotin
        'EX_pnto-R(e)'; ... % pantothenate
        'EX_inost(e)'; ... % inositol
        'EX_thm(e)'; ... % thiamine
        'EX_his-L(e)'; ... % histidine
        'EX_ura(e)'; ... % uracil
        'EX_met-L(e)'; ... % methionine
        'EX_leu-L(e)'; ... % leucine
        'EX_so4(e)'; ... % sulfate
        'EX_pi(e)'; ... % phosphate
        'EX_na1(e)'; ... % sodium
        'EX_k(e)'; ... % potassium
        'EX_o2(e)'; ... % oxygen
        'EX_h2o(e)'; ... % water
        'EX_h(e)'; ... % protons
        'EX_nh4(e)'; ... % ammonium
        'EX_thmpp(e)'; ... % thiamine diphosphate
        'EX_fe2(e)'; ... % 'Fe2 exchange';
        };

    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    if length(uptakeRxnIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
end

function model = kuepfer_iMMbs(model)
    % change iMM model media to kuepfer
    
    % start with a clean slate: set to unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'EX_btn(e)'; ... % biotin
        'EX_pnto_R(e)'; ... % pantothenate
        'EX_inost(e)'; ... % inositol
        'EX_thm(e)'; ... % thiamine
        'EX_his_L(e)'; ... % histidine
        'EX_ura(e)'; ... % uracil
        'EX_met_L(e)'; ... % methionine
        'EX_leu_L(e)'; ... % leucine
        'EX_so4(e)'; ... % sulfate
        'EX_pi(e)'; ... % phosphate
        'EX_na1(e)'; ... % sodium
        'EX_k(e)'; ... % potassium
        'EX_o2(e)'; ... % oxygen
        'EX_h2o(e)'; ... % water
        'EX_h(e)'; ... % protons
        'EX_nh4(e)'; ... % ammonium
        'EX_thmpp(e)'; ... % thiamine diphosphate
        'EX_fe2(e)'; ... % 'Fe2 exchange';
        };

    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    if length(uptakeRxnIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
end
                
function model = kuepfer_iAZ(model)
    % change iAZ model media to kuepfer
    
    % start with a clean slate: set to unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'EX_fe2(e)'; ... % iron
        'EX_btn(e)'; ... % biotin
        'EX_pnto_R(e)'; ... % pantothenate
        'EX_inost(e)'; ... % inositol
        'EX_thm(e)'; ... % thiamine
        'EX_his_L(e)'; ... % histidine
        'EX_ura(e)'; ... % uracil
        'EX_met_L(e)'; ... % methionine
        'EX_leu_L(e)'; ... % leucine
        'EX_so4(e)'; ... % sulfate
        'EX_pi(e)'; ... % phosphate
        'EX_na1(e)'; ... % sodium
        'EX_k(e)'; ... % potassium
        'EX_o2(e)'; ... % oxygen
        'EX_h2o(e)'; ... % water
        'EX_h(e)'; ... % protons
        'EX_nh4(e)'; ... % ammonium
        'EX_thmpp(e)'; ... % thiamine diphosphate
        % 'EX_lys_L(e)'; ... % lysine - to test
        };

    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    if length(uptakeRxnIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
end

function model = kuepfer_Y4(model)
    % change Y4 model media to kuepfer
        
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
        'r_1157'; ... % ammonia want ub=1000
        'r_1170'; ... % biotin ub
        'r_1335'; ... % inositol ub
        'r_1336'; ... % iron ub
        'r_1373'; ... % histidine ub
        'r_1384'; ... % leucine ub
        'r_1388'; ... % methionine ub
        'r_1427'; ... % nicotinate ub
        'r_1435'; ... % oxygen ub
        'r_1451'; ... % pantothenate ub
        'r_1461'; ... % phosphate ub
        'r_1476'; ... % potassium ub
        'r_1481'; ... % pyridoxine ub
        'r_1507'; ... % sulfate  ub
        'r_1512'; ... % thiamine ub
        'r_1527'; ... % uracil ub
        };
        
    desiredNegExchanges = {...
        'r_1497'; ... % sodium lb
        };
    
    desiredPosExchangesIndexes = findRxnIDs(model,desiredPosExchanges);
    desiredNegExchangesIndexes = findRxnIDs(model,desiredNegExchanges);
    
    if length(desiredPosExchangesIndexes) ~= length(desiredPosExchanges);
        error('Not all exchange reactions were found.')
    end
    
    if length(desiredNegExchangesIndexes) ~= length(desiredNegExchanges);
        error('Not all exchange reactions were found.')
    end
    
    model.ub(desiredPosExchangesIndexes)=1000;
    model.lb(desiredNegExchangesIndexes)=-1000;
end

function model = kuepfer_Y5(model)
    % change Y5 model media to kuepfer
    
    % start with a clean slate: set to unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'r_1861'; ... % iron
        'r_1671'; ... % biotin
        'r_1548'; ... % pantothenate
        'r_1967'; ... % nicotinic acid
        'r_1947'; ... % inositol
        'r_2028'; ... % pyridoxine
        'r_2067'; ... % thiamine
        'r_1893'; ... % histidine
        'r_2090'; ... % uracil
        'r_1902'; ... % methionine
        'r_1899'; ... % leucine
        'r_2060'; ... % sulfate
        'r_2005'; ... % phosphate
        'r_2049'; ... % sodium
        'r_2020'; ... % potassium
        'r_1992'; ... % oxygen
        'r_2100'; ... % water
        'r_1832'; ... % protons
        'r_1654'; ... % ammonium
        };

    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    if length(uptakeRxnIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
end

function model = kuepfer_Y6(model)
    % change Y6 model media to kuepfer
    
    % start with a clean slate: set to unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'r_1861'; ... % iron
        'r_1671'; ... % biotin
        'r_1548'; ... % pantothenate
        'r_1967'; ... % nicotinic acid
        'r_1947'; ... % inositol
        'r_2028'; ... % pyridoxine
        'r_2067'; ... % thiamine
        'r_1893'; ... % histidine
        'r_2090'; ... % uracil
        'r_1902'; ... % methionine
        'r_1899'; ... % leucine
        'r_2060'; ... % sulfate
        'r_2005'; ... % phosphate
        'r_2049'; ... % sodium
        %'r_2020'; ... % potassium -- removed from Y6 b/c a dead end
        'r_1992'; ... % oxygen
        'r_2100'; ... % water
        'r_1832'; ... % protons
        'r_1654'; ... % ammonium
        % 'r_1900'; ... % lysine - to test
        };
    
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    if length(uptakeRxnIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
end
   
function model = kuepfer_iIN(model)
    % change iIN model to kuepfer
    
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
        % iron
        'BTxtI'; ... % biotin
        'PNTOxtI'; ... % pantothenate
        'NAxtI'; ... % nicotinic acid
        'INSxtI'; ... % inositol
        % pyridoxine
        'THMxtI'; ... % thiamine
        'HISxtI'; ... % histidine
        'URAxtI'; ... % uracil
        'METxtI'; ... % methionine
        'LEUxtI'; ... % leucine
        'SLFxtI'; ... % sulfate
        'PIxtI'; ... % phosphate
        'NAxtI'; ... % sodium
        % potassium
        'O2xtI'; ... % oxygen
        % water
        % protons
        'NH3xtI'; ... % ammonium
        };

    desiredExchangesIndexes = findRxnIDs(model,desiredExchanges);
    if length(desiredExchangesIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.ub(desiredExchangesIndexes)=1000;
end

function model = kuepfer_iTO(model)
    % change iIN model media to minimal - ammonium, glucose, oxygen,
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
        % iron
        'BTxtI'; ... % biotin
        'PNTOxtI'; ... % pantothenate
        'NAxtI'; ... % nicotinic acid
        'INSxtI'; ... % inositol
        'VB6xtI'; ... % pyridoxine
        'THMxtI'; ... % thiamine
        'HISxtI'; ... % histidine
        'URAxtI'; ... % uracil
        'METxtI'; ... % methionine
        'LEUxtI'; ... % leucine
        'SLFxtI'; ... % sulfate
        'PIxtI'; ... % phosphate
        'NAxtI'; ... % sodium
        'KxtI'; ... % potassium
        'O2xtI'; ... % oxygen
        % water
        % protons
        'NH3xtI'; ... % ammonium
        };

    desiredExchangesIndexes = findRxnIDs(model,desiredExchanges);
    if length(desiredExchangesIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.ub(desiredExchangesIndexes)=1000;
end

function model = kuepfer_Y7(model)
    % change Y7 model media to kuepfer
    
    % start with a clean slate: set to unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
    desiredExchanges = {...
        'r_1861'; ... % iron
        'r_1671'; ... % biotin
        'r_1548'; ... % pantothenate
        'r_1967'; ... % nicotinic acid
        'r_1947'; ... % inositol
        'r_2028'; ... % pyridoxine
        'r_2067'; ... % thiamine
        'r_1893'; ... % histidine
        'r_2090'; ... % uracil
        'r_1902'; ... % methionine
        'r_1899'; ... % leucine
        'r_2060'; ... % sulfate
        'r_2005'; ... % phosphate
        'r_2049'; ... % sodium
        'r_2020'; ... % potassium
        'r_1992'; ... % oxygen
        'r_2100'; ... % water
        'r_1832'; ... % protons
        'r_1654'; ... % ammonium
        };
    
    uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
    if length(uptakeRxnIndexes) ~= length(desiredExchanges);
        error('Not all exchange reactions were found.')
    end
    model.lb(uptakeRxnIndexes)=-1000;
end

function model = kuepfer_bio(model)
    % change bio model media to kupefer medium
    
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
    
    desiredExchanges_in = {...
        'MNXM58_in'; ... % sulfate
        'MNXM9_in'; ... % phosphate
        'MNXM43_in'; ... % chloride
        'MNXM27_in'; ... % sodium
        'MNXM95_in'; ... % potassium
        'MNXM4_in'; ... % oxygen
        'MNXM2_in'; ... % water
        'MNXM1_in'; ... % protons
        'MNXM15_in'; ... % NH4
        };
    
    desiredExchanges_out = {...
        'bigg_btn_out'; ... % biotin
        'bigg_pnto_R_out'; ... % pantothenate
        'bigg_nac_out'; ... % nicotinic acid
        'bigg_inost_out'; ... % inositol
        'bigg_pydxn_out'; ... % pyridoxine
        'MNXM322_out'; ... % thiamine
        'bigg_his_L_out'; ... % histidine
        'bigg_ura_out'; ... % uracil
        'bigg_met_L_out'; ... % methionine
        'bigg_leu_L_out'; ... % leucine
        'bigg_ncam_out'; ... % nicotinamide
        'bigg_thmpp_out'; ... thiamine diphosphate
        };
    
    uptakeRxnIndexes_in = findRxnIDs(model,desiredExchanges_in);
    uptakeRxnIndexes_out = findRxnIDs(model,desiredExchanges_out);
    
    if length(uptakeRxnIndexes_in) ~= length(desiredExchanges_in);
        error('Not all exchange reactions were found.')
    end
    if length(uptakeRxnIndexes_out) ~= length(desiredExchanges_out);
        error('Not all exchange reactions were found.')
    end
    
    model.ub(uptakeRxnIndexes_in)=1000;
    model.lb(uptakeRxnIndexes_out)=-1000;
    
%     % add L-tyrosine by allowing the bigg_tyr_L_out rxn to be reversible
%     model.lb(findRxnIDs(model,'bigg_tyr_L_out')) = -1000;
%     
%     % add L-lysine by allowing the met bigg_lys_L_out rxn to be reversible
%     model.lb(findRxnIDs(model,'bigg_lys_L_out')) = -1000;
% 
%     % add L-isoleucine by allowing the met bigg_ile_L_out rxn to be
%     % reversible
%     model.lb(findRxnIDs(model,'bigg_ile_L_out')) = -1000;
% 
%     % add L-arginine by allowing the met bigg_arg_L_out rxn to be
%     % reversible
%     model.lb(findRxnIDs(model,'bigg_arg_L_out')) = -1000;
% 
%     % add L-histidine by allowing the met bigg_his_L_out rxn to be
%     % reversible
%     model.lb(findRxnIDs(model,'bigg_his_L_out')) = -1000;
%     
%     % add L-methionine by allowing the met bigg_met_L_out rxn to be
%     % reversible
%     model.lb(findRxnIDs(model,'bigg_met_L_out')) = -1000;
%     
%     % add L-tryptophan by allowing the met bigg_trp_L_out rxn to be
%     % reversible
%     model.lb(findRxnIDs(model,'bigg_trp_L_out')) = -1000;
end

%% gene lists
% BH is Ben Heavner (bheavner@gmail.com), HA is Hnin Aung
% (hnin.w.aung@gmail.com)

function genes = verifiedORFs
    % list of 4992 verified ORFs taken from SGD (6 Feb 2013)
    % http://www.yeastgenome.org/cgi-bin/search/featureSearch?featuretype=ORF&qualifier=Verified
    % Note: this still includes some ORFs with "dubious" function

    genes = {'Q0045';'Q0050';'Q0055';'Q0060';'Q0065';'Q0070';'Q0080';'Q0085';'Q0105';'Q0110';'Q0115';'Q0120';'Q0130';'Q0140';'Q0160';'Q0250';'Q0275';'R0010W';'R0020C';'R0030W';'R0040C';'YAL001C';'YAL002W';'YAL003W';'YAL005C';'YAL007C';'YAL008W';'YAL009W';'YAL010C';'YAL011W';'YAL012W';'YAL013W';'YAL014C';'YAL015C';'YAL016W';'YAL017W';'YAL019W';'YAL020C';'YAL021C';'YAL022C';'YAL023C';'YAL024C';'YAL025C';'YAL026C';'YAL027W';'YAL028W';'YAL029C';'YAL030W';'YAL031C';'YAL032C';'YAL033W';'YAL034C';'YAL034W-A';'YAL035W';'YAL036C';'YAL038W';'YAL039C';'YAL040C';'YAL041W';'YAL042W';'YAL043C';'YAL044C';'YAL046C';'YAL047C';'YAL048C';'YAL049C';'YAL051W';'YAL053W';'YAL054C';'YAL055W';'YAL056W';'YAL058W';'YAL059W';'YAL060W';'YAL062W';'YAL063C';'YAL064W';'YAL067C';'YAL068C';'YAR002C-A';'YAR002W';'YAR003W';'YAR007C';'YAR008W';'YAR014C';'YAR015W';'YAR018C';'YAR019C';'YAR020C';'YAR027W';'YAR031W';'YAR033W';'YAR035W';'YAR042W';'YAR050W';'YAR071W';'YBL001C';'YBL002W';'YBL003C';'YBL004W';'YBL005W';'YBL006C';'YBL007C';'YBL008W';'YBL009W';'YBL011W';'YBL013W';'YBL014C';'YBL015W';'YBL016W';'YBL017C';'YBL018C';'YBL019W';'YBL020W';'YBL021C';'YBL022C';'YBL023C';'YBL024W';'YBL025W';'YBL026W';'YBL027W';'YBL028C';'YBL030C';'YBL031W';'YBL032W';'YBL033C';'YBL034C';'YBL035C';'YBL036C';'YBL037W';'YBL038W';'YBL039C';'YBL040C';'YBL041W';'YBL042C';'YBL043W';'YBL045C';'YBL046W';'YBL047C';'YBL049W';'YBL050W';'YBL051C';'YBL052C';'YBL054W';'YBL055C';'YBL056W';'YBL057C';'YBL058W';'YBL059C-A';'YBL060W';'YBL061C';'YBL063W';'YBL064C';'YBL066C';'YBL067C';'YBL068W';'YBL069W';'YBL071W-A';'YBL072C';'YBL074C';'YBL075C';'YBL076C';'YBL078C';'YBL079W';'YBL080C';'YBL082C';'YBL084C';'YBL085W';'YBL087C';'YBL088C';'YBL089W';'YBL090W';'YBL091C';'YBL091C-A';'YBL092W';'YBL093C';'YBL097W';'YBL098W';'YBL099W';'YBL101C';'YBL102W';'YBL103C';'YBL104C';'YBL105C';'YBL106C';'YBL107C';'YBL108C-A';'YBR001C';'YBR002C';'YBR003W';'YBR004C';'YBR005W';'YBR006W';'YBR008C';'YBR009C';'YBR010W';'YBR011C';'YBR014C';'YBR015C';'YBR016W';'YBR017C';'YBR018C';'YBR019C';'YBR020W';'YBR021W';'YBR022W';'YBR023C';'YBR024W';'YBR025C';'YBR026C';'YBR028C';'YBR029C';'YBR030W';'YBR031W';'YBR034C';'YBR035C';'YBR036C';'YBR037C';'YBR038W';'YBR039W';'YBR040W';'YBR041W';'YBR042C';'YBR043C';'YBR044C';'YBR045C';'YBR046C';'YBR048W';'YBR049C';'YBR050C';'YBR052C';'YBR054W';'YBR055C';'YBR056W';'YBR057C';'YBR058C';'YBR058C-A';'YBR059C';'YBR060C';'YBR061C';'YBR065C';'YBR066C';'YBR067C';'YBR068C';'YBR069C';'YBR070C';'YBR071W';'YBR072W';'YBR073W';'YBR076W';'YBR077C';'YBR078W';'YBR079C';'YBR080C';'YBR081C';'YBR082C';'YBR083W';'YBR084C-A';'YBR084W';'YBR085W';'YBR086C';'YBR087W';'YBR088C';'YBR089C-A';'YBR091C';'YBR092C';'YBR093C';'YBR094W';'YBR095C';'YBR097W';'YBR098W';'YBR101C';'YBR102C';'YBR103W';'YBR104W';'YBR105C';'YBR106W';'YBR107C';'YBR108W';'YBR109C';'YBR110W';'YBR111C';'YBR111W-A';'YBR112C';'YBR114W';'YBR115C';'YBR117C';'YBR118W';'YBR119W';'YBR120C';'YBR121C';'YBR122C';'YBR123C';'YBR125C';'YBR126C';'YBR127C';'YBR128C';'YBR129C';'YBR130C';'YBR131W';'YBR132C';'YBR133C';'YBR135W';'YBR136W';'YBR137W';'YBR139W';'YBR140C';'YBR142W';'YBR143C';'YBR145W';'YBR146W';'YBR147W';'YBR148W';'YBR149W';'YBR150C';'YBR151W';'YBR152W';'YBR153W';'YBR154C';'YBR155W';'YBR156C';'YBR157C';'YBR158W';'YBR159W';'YBR160W';'YBR161W';'YBR162C';'YBR162W-A';'YBR163W';'YBR164C';'YBR165W';'YBR166C';'YBR167C';'YBR168W';'YBR169C';'YBR170C';'YBR171W';'YBR172C';'YBR173C';'YBR175W';'YBR176W';'YBR177C';'YBR179C';'YBR180W';'YBR181C';'YBR182C';'YBR183W';'YBR185C';'YBR186W';'YBR188C';'YBR189W';'YBR191W';'YBR192W';'YBR193C';'YBR194W';'YBR195C';'YBR196C';'YBR198C';'YBR199W';'YBR200W';'YBR201W';'YBR202W';'YBR203W';'YBR204C';'YBR205W';'YBR207W';'YBR208C';'YBR210W';'YBR211C';'YBR212W';'YBR213W';'YBR214W';'YBR215W';'YBR216C';'YBR217W';'YBR218C';'YBR221C';'YBR222C';'YBR223C';'YBR227C';'YBR228W';'YBR229C';'YBR230C';'YBR231C';'YBR233W';'YBR233W-A';'YBR234C';'YBR235W';'YBR236C';'YBR237W';'YBR238C';'YBR240C';'YBR243C';'YBR244W';'YBR245C';'YBR246W';'YBR247C';'YBR248C';'YBR249C';'YBR250W';'YBR251W';'YBR252W';'YBR253W';'YBR254C';'YBR255W';'YBR256C';'YBR257W';'YBR258C';'YBR260C';'YBR261C';'YBR262C';'YBR263W';'YBR264C';'YBR265W';'YBR267W';'YBR268W';'YBR271W';'YBR272C';'YBR273C';'YBR274W';'YBR275C';'YBR276C';'YBR278W';'YBR279W';'YBR280C';'YBR281C';'YBR282W';'YBR283C';'YBR286W';'YBR288C';'YBR289W';'YBR290W';'YBR291C';'YBR293W';'YBR294W';'YBR295W';'YBR296C';'YBR297W';'YBR298C';'YBR299W';'YBR301W';'YBR302C';'YCL001W';'YCL004W';'YCL005W';'YCL005W-A';'YCL008C';'YCL009C';'YCL010C';'YCL011C';'YCL012C';'YCL014W';'YCL016C';'YCL017C';'YCL018W';'YCL024W';'YCL025C';'YCL026C-A';'YCL027W';'YCL028W';'YCL029C';'YCL030C';'YCL031C';'YCL032W';'YCL033C';'YCL034W';'YCL035C';'YCL036W';'YCL037C';'YCL038C';'YCL039W';'YCL040W';'YCL043C';'YCL044C';'YCL045C';'YCL047C';'YCL048W';'YCL050C';'YCL051W';'YCL052C';'YCL054W';'YCL055W';'YCL056C';'YCL057C-A';'YCL057W';'YCL058C';'YCL058W-A';'YCL059C';'YCL061C';'YCL063W';'YCL064C';'YCL066W';'YCL067C';'YCL069W';'YCL073C';'YCR002C';'YCR003W';'YCR004C';'YCR005C';'YCR008W';'YCR009C';'YCR010C';'YCR011C';'YCR012W';'YCR014C';'YCR017C';'YCR018C';'YCR019W';'YCR020C';'YCR020C-A';'YCR020W-B';'YCR021C';'YCR023C';'YCR024C';'YCR024C-A';'YCR026C';'YCR027C';'YCR028C';'YCR028C-A';'YCR030C';'YCR031C';'YCR032W';'YCR033W';'YCR034W';'YCR035C';'YCR036W';'YCR037C';'YCR038C';'YCR039C';'YCR040W';'YCR042C';'YCR044C';'YCR045C';'YCR046C';'YCR047C';'YCR048W';'YCR052W';'YCR053W';'YCR054C';'YCR057C';'YCR059C';'YCR060W';'YCR063W';'YCR065W';'YCR066W';'YCR067C';'YCR068W';'YCR069W';'YCR071C';'YCR072C';'YCR073C';'YCR073W-A';'YCR075C';'YCR076C';'YCR077C';'YCR079W';'YCR081W';'YCR082W';'YCR083W';'YCR084C';'YCR086W';'YCR088W';'YCR089W';'YCR091W';'YCR092C';'YCR093W';'YCR094W';'YCR096C';'YCR097W';'YCR098C';'YCR104W';'YCR105W';'YCR106W';'YCR107W';'YDL001W';'YDL002C';'YDL003W';'YDL004W';'YDL005C';'YDL006W';'YDL007W';'YDL008W';'YDL010W';'YDL012C';'YDL013W';'YDL014W';'YDL015C';'YDL017W';'YDL018C';'YDL019C';'YDL020C';'YDL021W';'YDL022W';'YDL024C';'YDL025C';'YDL028C';'YDL029W';'YDL030W';'YDL031W';'YDL033C';'YDL035C';'YDL036C';'YDL037C';'YDL039C';'YDL040C';'YDL042C';'YDL043C';'YDL044C';'YDL045C';'YDL045W-A';'YDL046W';'YDL047W';'YDL048C';'YDL049C';'YDL051W';'YDL052C';'YDL053C';'YDL054C';'YDL055C';'YDL056W';'YDL058W';'YDL059C';'YDL060W';'YDL061C';'YDL063C';'YDL064W';'YDL065C';'YDL066W';'YDL067C';'YDL069C';'YDL070W';'YDL072C';'YDL074C';'YDL075W';'YDL076C';'YDL077C';'YDL078C';'YDL079C';'YDL080C';'YDL081C';'YDL082W';'YDL083C';'YDL084W';'YDL085W';'YDL087C';'YDL088C';'YDL089W';'YDL090C';'YDL091C';'YDL092W';'YDL093W';'YDL095W';'YDL097C';'YDL098C';'YDL099W';'YDL100C';'YDL101C';'YDL102W';'YDL103C';'YDL104C';'YDL105W';'YDL106C';'YDL107W';'YDL108W';'YDL110C';'YDL111C';'YDL112W';'YDL113C';'YDL115C';'YDL116W';'YDL117W';'YDL120W';'YDL122W';'YDL123W';'YDL124W';'YDL125C';'YDL126C';'YDL127W';'YDL128W';'YDL130W';'YDL130W-A';'YDL131W';'YDL132W';'YDL133C-A';'YDL133W';'YDL134C';'YDL135C';'YDL136W';'YDL137W';'YDL138W';'YDL139C';'YDL140C';'YDL141W';'YDL142C';'YDL143W';'YDL145C';'YDL146W';'YDL147W';'YDL148C';'YDL149W';'YDL150W';'YDL153C';'YDL154W';'YDL155W';'YDL156W';'YDL159W';'YDL160C';'YDL160C-A';'YDL161W';'YDL164C';'YDL165W';'YDL166C';'YDL167C';'YDL168W';'YDL169C';'YDL170W';'YDL171C';'YDL173W';'YDL174C';'YDL175C';'YDL176W';'YDL178W';'YDL179W';'YDL181W';'YDL182W';'YDL183C';'YDL184C';'YDL185W';'YDL188C';'YDL189W';'YDL190C';'YDL191W';'YDL192W';'YDL193W';'YDL194W';'YDL195W';'YDL197C';'YDL198C';'YDL200C';'YDL201W';'YDL202W';'YDL203C';'YDL204W';'YDL205C';'YDL207W';'YDL208W';'YDL209C';'YDL210W';'YDL212W';'YDL213C';'YDL214C';'YDL215C';'YDL216C';'YDL217C';'YDL219W';'YDL220C';'YDL222C';'YDL223C';'YDL224C';'YDL225W';'YDL226C';'YDL227C';'YDL229W';'YDL230W';'YDL231C';'YDL232W';'YDL234C';'YDL235C';'YDL236W';'YDL237W';'YDL238C';'YDL239C';'YDL240W';'YDL243C';'YDL244W';'YDL245C';'YDL247W';'YDL248W';'YDR001C';'YDR002W';'YDR003W';'YDR004W';'YDR005C';'YDR006C';'YDR007W';'YDR009W';'YDR011W';'YDR012W';'YDR013W';'YDR014W';'YDR014W-A';'YDR016C';'YDR017C';'YDR019C';'YDR021W';'YDR022C';'YDR023W';'YDR025W';'YDR026C';'YDR027C';'YDR028C';'YDR030C';'YDR031W';'YDR032C';'YDR033W';'YDR034C';'YDR035W';'YDR036C';'YDR037W';'YDR038C';'YDR039C';'YDR040C';'YDR041W';'YDR043C';'YDR044W';'YDR045C';'YDR046C';'YDR047W';'YDR049W';'YDR050C';'YDR051C';'YDR052C';'YDR054C';'YDR055W';'YDR057W';'YDR058C';'YDR059C';'YDR060W';'YDR062W';'YDR063W';'YDR064W';'YDR065W';'YDR068W';'YDR069C';'YDR071C';'YDR072C';'YDR073W';'YDR074W';'YDR075W';'YDR076W';'YDR077W';'YDR078C';'YDR079C-A';'YDR079W';'YDR080W';'YDR081C';'YDR082W';'YDR083W';'YDR084C';'YDR085C';'YDR086C';'YDR087C';'YDR088C';'YDR091C';'YDR092W';'YDR093W';'YDR096W';'YDR097C';'YDR098C';'YDR099W';'YDR100W';'YDR101C';'YDR103W';'YDR104C';'YDR105C';'YDR106W';'YDR107C';'YDR108W';'YDR110W';'YDR113C';'YDR116C';'YDR117C';'YDR118W';'YDR119W-A';'YDR120C';'YDR121W';'YDR122W';'YDR123C';'YDR125C';'YDR126W';'YDR127W';'YDR128W';'YDR129C';'YDR130C';'YDR135C';'YDR137W';'YDR138W';'YDR139C';'YDR140W';'YDR141C';'YDR142C';'YDR143C';'YDR144C';'YDR145W';'YDR146C';'YDR147W';'YDR148C';'YDR150W';'YDR151C';'YDR152W';'YDR153C';'YDR155C';'YDR156W';'YDR158W';'YDR159W';'YDR160W';'YDR162C';'YDR163W';'YDR164C';'YDR165W';'YDR166C';'YDR167W';'YDR168W';'YDR169C';'YDR170C';'YDR171W';'YDR172W';'YDR173C';'YDR174W';'YDR175C';'YDR176W';'YDR177W';'YDR178W';'YDR179C';'YDR180W';'YDR181C';'YDR182W';'YDR183W';'YDR184C';'YDR185C';'YDR186C';'YDR188W';'YDR189W';'YDR190C';'YDR191W';'YDR192C';'YDR194C';'YDR195W';'YDR196C';'YDR197W';'YDR198C';'YDR200C';'YDR201W';'YDR202C';'YDR204W';'YDR205W';'YDR206W';'YDR207C';'YDR208W';'YDR211W';'YDR212W';'YDR213W';'YDR214W';'YDR216W';'YDR217C';'YDR218C';'YDR219C';'YDR221W';'YDR223W';'YDR224C';'YDR225W';'YDR226W';'YDR227W';'YDR228C';'YDR229W';'YDR231C';'YDR232W';'YDR233C';'YDR234W';'YDR235W';'YDR236C';'YDR237W';'YDR238C';'YDR239C';'YDR240C';'YDR242W';'YDR243C';'YDR244W';'YDR245W';'YDR246W';'YDR247W';'YDR251W';'YDR252W';'YDR253C';'YDR254W';'YDR255C';'YDR256C';'YDR257C';'YDR258C';'YDR259C';'YDR260C';'YDR261C';'YDR263C';'YDR264C';'YDR265W';'YDR266C';'YDR267C';'YDR268W';'YDR270W';'YDR272W';'YDR273W';'YDR275W';'YDR276C';'YDR277C';'YDR279W';'YDR280W';'YDR281C';'YDR283C';'YDR284C';'YDR285W';'YDR287W';'YDR288W';'YDR289C';'YDR291W';'YDR292C';'YDR293C';'YDR294C';'YDR295C';'YDR296W';'YDR297W';'YDR298C';'YDR299W';'YDR300C';'YDR301W';'YDR302W';'YDR303C';'YDR304C';'YDR305C';'YDR308C';'YDR309C';'YDR310C';'YDR311W';'YDR312W';'YDR313C';'YDR314C';'YDR315C';'YDR316W';'YDR317W';'YDR318W';'YDR320C';'YDR320C-A';'YDR321W';'YDR322C-A';'YDR322W';'YDR323C';'YDR324C';'YDR325W';'YDR326C';'YDR328C';'YDR329C';'YDR330W';'YDR331W';'YDR332W';'YDR333C';'YDR334W';'YDR335W';'YDR337W';'YDR339C';'YDR341C';'YDR342C';'YDR343C';'YDR345C';'YDR346C';'YDR347W';'YDR348C';'YDR349C';'YDR350C';'YDR351W';'YDR352W';'YDR353W';'YDR354W';'YDR356W';'YDR358W';'YDR359C';'YDR361C';'YDR362C';'YDR363W';'YDR363W-A';'YDR364C';'YDR365C';'YDR367W';'YDR368W';'YDR369C';'YDR372C';'YDR373W';'YDR374W-A';'YDR375C';'YDR376W';'YDR377W';'YDR378C';'YDR379C-A';'YDR379W';'YDR380W';'YDR381C-A';'YDR381W';'YDR382W';'YDR383C';'YDR384C';'YDR385W';'YDR386W';'YDR388W';'YDR389W';'YDR390C';'YDR392W';'YDR393W';'YDR394W';'YDR395W';'YDR397C';'YDR398W';'YDR399W';'YDR400W';'YDR402C';'YDR403W';'YDR404C';'YDR405W';'YDR406W';'YDR407C';'YDR408C';'YDR409W';'YDR410C';'YDR411C';'YDR412W';'YDR414C';'YDR416W';'YDR418W';'YDR419W';'YDR420W';'YDR421W';'YDR422C';'YDR423C';'YDR424C';'YDR425W';'YDR427W';'YDR428C';'YDR429C';'YDR430C';'YDR432W';'YDR434W';'YDR435C';'YDR436W';'YDR437W';'YDR438W';'YDR439W';'YDR440W';'YDR441C';'YDR443C';'YDR446W';'YDR447C';'YDR448W';'YDR449C';'YDR450W';'YDR451C';'YDR452W';'YDR453C';'YDR454C';'YDR456W';'YDR457W';'YDR458C';'YDR459C';'YDR460W';'YDR461W';'YDR462W';'YDR463W';'YDR464W';'YDR465C';'YDR466W';'YDR468C';'YDR469W';'YDR470C';'YDR471W';'YDR472W';'YDR473C';'YDR475C';'YDR477W';'YDR478W';'YDR479C';'YDR480W';'YDR481C';'YDR482C';'YDR483W';'YDR484W';'YDR485C';'YDR486C';'YDR487C';'YDR488C';'YDR489W';'YDR490C';'YDR492W';'YDR493W';'YDR494W';'YDR495C';'YDR496C';'YDR497C';'YDR498C';'YDR499W';'YDR500C';'YDR501W';'YDR502C';'YDR503C';'YDR504C';'YDR505C';'YDR506C';'YDR507C';'YDR508C';'YDR510W';'YDR511W';'YDR512C';'YDR513W';'YDR514C';'YDR515W';'YDR516C';'YDR517W';'YDR518W';'YDR519W';'YDR522C';'YDR523C';'YDR524C';'YDR525W-A';'YDR527W';'YDR528W';'YDR529C';'YDR530C';'YDR531W';'YDR532C';'YDR533C';'YDR534C';'YDR536W';'YDR538W';'YDR539W';'YDR540C';'YDR542W';'YDR545W';'YEL001C';'YEL002C';'YEL003W';'YEL004W';'YEL005C';'YEL006W';'YEL007W';'YEL009C';'YEL011W';'YEL012W';'YEL013W';'YEL015W';'YEL016C';'YEL017C-A';'YEL017W';'YEL018W';'YEL019C';'YEL020W-A';'YEL021W';'YEL022W';'YEL024W';'YEL026W';'YEL027W';'YEL029C';'YEL030W';'YEL031W';'YEL032W';'YEL034W';'YEL036C';'YEL037C';'YEL038W';'YEL039C';'YEL040W';'YEL041W';'YEL042W';'YEL043W';'YEL044W';'YEL046C';'YEL047C';'YEL048C';'YEL049W';'YEL050C';'YEL051W';'YEL052W';'YEL053C';'YEL054C';'YEL055C';'YEL056W';'YEL058W';'YEL059C-A';'YEL060C';'YEL061C';'YEL062W';'YEL063C';'YEL064C';'YEL065W';'YEL066W';'YEL069C';'YEL071W';'YEL072W';'YER001W';'YER002W';'YER003C';'YER004W';'YER005W';'YER006W';'YER007C-A';'YER007W';'YER008C';'YER009W';'YER010C';'YER011W';'YER012W';'YER013W';'YER014C-A';'YER014W';'YER015W';'YER016W';'YER017C';'YER018C';'YER019C-A';'YER019W';'YER020W';'YER021W';'YER022W';'YER023W';'YER024W';'YER025W';'YER026C';'YER027C';'YER028C';'YER029C';'YER030W';'YER031C';'YER032W';'YER033C';'YER035W';'YER036C';'YER037W';'YER038C';'YER039C';'YER040W';'YER041W';'YER042W';'YER043C';'YER044C';'YER044C-A';'YER045C';'YER046W';'YER047C';'YER048C';'YER048W-A';'YER049W';'YER050C';'YER051W';'YER052C';'YER053C';'YER054C';'YER055C';'YER056C';'YER056C-A';'YER057C';'YER058W';'YER059W';'YER060W';'YER060W-A';'YER061C';'YER062C';'YER063W';'YER065C';'YER067W';'YER068W';'YER069W';'YER070W';'YER072W';'YER073W';'YER074W';'YER074W-A';'YER075C';'YER078C';'YER080W';'YER081W';'YER082C';'YER083C';'YER086W';'YER087C-B';'YER087W';'YER088C';'YER089C';'YER090W';'YER091C';'YER092W';'YER093C';'YER093C-A';'YER094C';'YER095W';'YER096W';'YER098W';'YER099C';'YER100W';'YER101C';'YER102W';'YER103W';'YER104W';'YER105C';'YER106W';'YER107C';'YER109C';'YER110C';'YER111C';'YER112W';'YER113C';'YER114C';'YER115C';'YER116C';'YER117W';'YER118C';'YER119C';'YER120W';'YER122C';'YER123W';'YER124C';'YER125W';'YER126C';'YER127W';'YER128W';'YER129W';'YER131W';'YER132C';'YER133W';'YER134C';'YER136W';'YER139C';'YER140W';'YER141W';'YER142C';'YER143W';'YER144C';'YER145C';'YER146W';'YER147C';'YER148W';'YER149C';'YER150W';'YER151C';'YER152C';'YER153C';'YER154W';'YER155C';'YER157W';'YER159C';'YER161C';'YER162C';'YER163C';'YER164W';'YER165W';'YER166W';'YER167W';'YER168C';'YER169W';'YER170W';'YER171W';'YER172C';'YER173W';'YER174C';'YER175C';'YER176W';'YER177W';'YER178W';'YER179W';'YER180C';'YER180C-A';'YER183C';'YER185W';'YER190W';'YFL001W';'YFL002C';'YFL003C';'YFL004W';'YFL005W';'YFL007W';'YFL008W';'YFL009W';'YFL010C';'YFL010W-A';'YFL011W';'YFL013C';'YFL014W';'YFL016C';'YFL017C';'YFL017W-A';'YFL018C';'YFL020C';'YFL021W';'YFL022C';'YFL023W';'YFL024C';'YFL025C';'YFL026W';'YFL027C';'YFL028C';'YFL029C';'YFL030W';'YFL031W';'YFL033C';'YFL034C-A';'YFL034C-B';'YFL036W';'YFL037W';'YFL038C';'YFL039C';'YFL041W';'YFL044C';'YFL045C';'YFL047W';'YFL048C';'YFL049W';'YFL050C';'YFL053W';'YFL055W';'YFL056C';'YFL057C';'YFL058W';'YFL059W';'YFL060C';'YFL062W';'YFR001W';'YFR002W';'YFR003C';'YFR004W';'YFR005C';'YFR007W';'YFR008W';'YFR009W';'YFR010W';'YFR011C';'YFR012W';'YFR013W';'YFR014C';'YFR015C';'YFR016C';'YFR017C';'YFR019W';'YFR021W';'YFR022W';'YFR023W';'YFR024C-A';'YFR025C';'YFR026C';'YFR027W';'YFR028C';'YFR029W';'YFR030W';'YFR031C';'YFR031C-A';'YFR032C-A';'YFR033C';'YFR034C';'YFR036W';'YFR037C';'YFR038W';'YFR040W';'YFR041C';'YFR042W';'YFR043C';'YFR044C';'YFR046C';'YFR047C';'YFR048W';'YFR049W';'YFR050C';'YFR051C';'YFR052W';'YFR053C';'YGL001C';'YGL002W';'YGL003C';'YGL004C';'YGL005C';'YGL006W';'YGL008C';'YGL009C';'YGL011C';'YGL012W';'YGL013C';'YGL014W';'YGL016W';'YGL017W';'YGL018C';'YGL019W';'YGL020C';'YGL021W';'YGL022W';'YGL023C';'YGL025C';'YGL026C';'YGL027C';'YGL028C';'YGL029W';'YGL030W';'YGL031C';'YGL032C';'YGL033W';'YGL035C';'YGL037C';'YGL038C';'YGL039W';'YGL040C';'YGL043W';'YGL044C';'YGL045W';'YGL047W';'YGL048C';'YGL049C';'YGL050W';'YGL051W';'YGL053W';'YGL054C';'YGL055W';'YGL056C';'YGL057C';'YGL058W';'YGL059W';'YGL060W';'YGL061C';'YGL062W';'YGL063W';'YGL064C';'YGL065C';'YGL066W';'YGL067W';'YGL068W';'YGL070C';'YGL071W';'YGL073W';'YGL075C';'YGL076C';'YGL077C';'YGL078C';'YGL080W';'YGL083W';'YGL084C';'YGL086W';'YGL087C';'YGL089C';'YGL090W';'YGL091C';'YGL092W';'YGL093W';'YGL094C';'YGL095C';'YGL096W';'YGL097W';'YGL098W';'YGL099W';'YGL100W';'YGL103W';'YGL104C';'YGL105W';'YGL106W';'YGL107C';'YGL110C';'YGL111W';'YGL112C';'YGL113W';'YGL115W';'YGL116W';'YGL119W';'YGL120C';'YGL121C';'YGL122C';'YGL123W';'YGL124C';'YGL125W';'YGL126W';'YGL127C';'YGL128C';'YGL129C';'YGL130W';'YGL131C';'YGL133W';'YGL134W';'YGL135W';'YGL136C';'YGL137W';'YGL139W';'YGL141W';'YGL142C';'YGL143C';'YGL144C';'YGL145W';'YGL147C';'YGL148W';'YGL150C';'YGL151W';'YGL153W';'YGL154C';'YGL155W';'YGL156W';'YGL157W';'YGL158W';'YGL160W';'YGL161C';'YGL162W';'YGL163C';'YGL164C';'YGL166W';'YGL167C';'YGL168W';'YGL169W';'YGL170C';'YGL171W';'YGL172W';'YGL173C';'YGL174W';'YGL175C';'YGL178W';'YGL179C';'YGL180W';'YGL181W';'YGL183C';'YGL184C';'YGL186C';'YGL187C';'YGL189C';'YGL190C';'YGL191W';'YGL192W';'YGL194C';'YGL195W';'YGL196W';'YGL197W';'YGL198W';'YGL200C';'YGL201C';'YGL202W';'YGL203C';'YGL205W';'YGL206C';'YGL207W';'YGL208W';'YGL209W';'YGL210W';'YGL211W';'YGL212W';'YGL213C';'YGL215W';'YGL216W';'YGL219C';'YGL220W';'YGL221C';'YGL222C';'YGL223C';'YGL224C';'YGL225W';'YGL226C-A';'YGL226W';'YGL227W';'YGL228W';'YGL229C';'YGL231C';'YGL232W';'YGL233W';'YGL234W';'YGL236C';'YGL237C';'YGL238W';'YGL240W';'YGL241W';'YGL243W';'YGL244W';'YGL245W';'YGL246C';'YGL247W';'YGL248W';'YGL249W';'YGL250W';'YGL251C';'YGL252C';'YGL253W';'YGL254W';'YGL255W';'YGL256W';'YGL257C';'YGL258W';'YGL263W';'YGR002C';'YGR003W';'YGR004W';'YGR005C';'YGR006W';'YGR007W';'YGR008C';'YGR009C';'YGR010W';'YGR012W';'YGR013W';'YGR014W';'YGR019W';'YGR020C';'YGR023W';'YGR024C';'YGR027C';'YGR028W';'YGR029W';'YGR030C';'YGR031C-A';'YGR031W';'YGR032W';'YGR033C';'YGR034W';'YGR036C';'YGR037C';'YGR038W';'YGR040W';'YGR041W';'YGR043C';'YGR044C';'YGR046W';'YGR047C';'YGR048W';'YGR049W';'YGR054W';'YGR055W';'YGR056W';'YGR057C';'YGR058W';'YGR059W';'YGR060W';'YGR061C';'YGR062C';'YGR063C';'YGR065C';'YGR068C';'YGR070W';'YGR071C';'YGR072W';'YGR074W';'YGR075C';'YGR076C';'YGR077C';'YGR078C';'YGR080W';'YGR081C';'YGR082W';'YGR083C';'YGR084C';'YGR085C';'YGR086C';'YGR087C';'YGR088W';'YGR089W';'YGR090W';'YGR091W';'YGR092W';'YGR094W';'YGR095C';'YGR096W';'YGR097W';'YGR098C';'YGR099W';'YGR100W';'YGR101W';'YGR102C';'YGR103W';'YGR104C';'YGR105W';'YGR106C';'YGR108W';'YGR109C';'YGR110W';'YGR112W';'YGR113W';'YGR116W';'YGR118W';'YGR119C';'YGR120C';'YGR121C';'YGR122W';'YGR123C';'YGR124W';'YGR128C';'YGR129W';'YGR130C';'YGR131W';'YGR132C';'YGR133W';'YGR134W';'YGR135W';'YGR136W';'YGR138C';'YGR140W';'YGR141W';'YGR142W';'YGR143W';'YGR144W';'YGR145W';'YGR146C';'YGR147C';'YGR148C';'YGR150C';'YGR152C';'YGR154C';'YGR155W';'YGR156W';'YGR157W';'YGR158C';'YGR159C';'YGR162W';'YGR163W';'YGR165W';'YGR166W';'YGR167W';'YGR169C';'YGR170W';'YGR171C';'YGR172C';'YGR173W';'YGR174C';'YGR175C';'YGR177C';'YGR178C';'YGR179C';'YGR180C';'YGR181W';'YGR183C';'YGR184C';'YGR185C';'YGR186W';'YGR187C';'YGR188C';'YGR189C';'YGR191W';'YGR192C';'YGR193C';'YGR194C';'YGR195W';'YGR196C';'YGR197C';'YGR198W';'YGR199W';'YGR200C';'YGR202C';'YGR203W';'YGR204W';'YGR205W';'YGR206W';'YGR207C';'YGR208W';'YGR209C';'YGR211W';'YGR212W';'YGR213C';'YGR214W';'YGR215W';'YGR216C';'YGR217W';'YGR218W';'YGR220C';'YGR221C';'YGR222W';'YGR223C';'YGR224W';'YGR225W';'YGR227W';'YGR229C';'YGR230W';'YGR231C';'YGR232W';'YGR233C';'YGR234W';'YGR235C';'YGR236C';'YGR238C';'YGR239C';'YGR240C';'YGR241C';'YGR243W';'YGR244C';'YGR245C';'YGR246C';'YGR247W';'YGR248W';'YGR249W';'YGR250C';'YGR251W';'YGR252W';'YGR253C';'YGR254W';'YGR255C';'YGR256W';'YGR257C';'YGR258C';'YGR260W';'YGR261C';'YGR262C';'YGR263C';'YGR264C';'YGR266W';'YGR267C';'YGR268C';'YGR270W';'YGR271C-A';'YGR271W';'YGR274C';'YGR275W';'YGR276C';'YGR277C';'YGR278W';'YGR279C';'YGR280C';'YGR281W';'YGR282C';'YGR283C';'YGR284C';'YGR285C';'YGR286C';'YGR287C';'YGR288W';'YGR289C';'YGR292W';'YGR294W';'YGR295C';'YGR296W';'YHL001W';'YHL002W';'YHL003C';'YHL004W';'YHL006C';'YHL007C';'YHL009C';'YHL010C';'YHL011C';'YHL013C';'YHL014C';'YHL015W';'YHL016C';'YHL019C';'YHL020C';'YHL021C';'YHL022C';'YHL023C';'YHL024W';'YHL025W';'YHL027W';'YHL028W';'YHL030W';'YHL031C';'YHL032C';'YHL033C';'YHL034C';'YHL035C';'YHL036W';'YHL038C';'YHL039W';'YHL040C';'YHL043W';'YHL046C';'YHL047C';'YHL048W';'YHR001W';'YHR001W-A';'YHR002W';'YHR003C';'YHR004C';'YHR005C';'YHR005C-A';'YHR006W';'YHR007C';'YHR008C';'YHR010W';'YHR011W';'YHR012W';'YHR013C';'YHR014W';'YHR015W';'YHR016C';'YHR017W';'YHR018C';'YHR019C';'YHR020W';'YHR021C';'YHR023W';'YHR024C';'YHR025W';'YHR026W';'YHR027C';'YHR028C';'YHR029C';'YHR030C';'YHR031C';'YHR032W';'YHR034C';'YHR036W';'YHR037W';'YHR038W';'YHR039C';'YHR039C-A';'YHR040W';'YHR041C';'YHR042W';'YHR043C';'YHR044C';'YHR046C';'YHR047C';'YHR049W';'YHR050W';'YHR051W';'YHR052W';'YHR053C';'YHR055C';'YHR056C';'YHR057C';'YHR058C';'YHR059W';'YHR060W';'YHR061C';'YHR062C';'YHR063C';'YHR064C';'YHR065C';'YHR066W';'YHR067W';'YHR068W';'YHR069C';'YHR070W';'YHR071W';'YHR072W';'YHR072W-A';'YHR073W';'YHR074W';'YHR075C';'YHR076W';'YHR077C';'YHR079C';'YHR079C-A';'YHR080C';'YHR081W';'YHR082C';'YHR083W';'YHR084W';'YHR085W';'YHR086W';'YHR087W';'YHR088W';'YHR089C';'YHR090C';'YHR091C';'YHR092C';'YHR094C';'YHR096C';'YHR098C';'YHR099W';'YHR100C';'YHR101C';'YHR102W';'YHR103W';'YHR104W';'YHR105W';'YHR106W';'YHR107C';'YHR108W';'YHR109W';'YHR110W';'YHR111W';'YHR112C';'YHR113W';'YHR114W';'YHR115C';'YHR116W';'YHR117W';'YHR118C';'YHR119W';'YHR120W';'YHR121W';'YHR122W';'YHR123W';'YHR124W';'YHR127W';'YHR128W';'YHR129C';'YHR132C';'YHR132W-A';'YHR133C';'YHR134W';'YHR135C';'YHR136C';'YHR137W';'YHR139C';'YHR141C';'YHR142W';'YHR143W';'YHR143W-A';'YHR144C';'YHR146W';'YHR147C';'YHR148W';'YHR149C';'YHR150W';'YHR151C';'YHR152W';'YHR153C';'YHR154W';'YHR155W';'YHR156C';'YHR157W';'YHR158C';'YHR160C';'YHR161C';'YHR162W';'YHR163W';'YHR164C';'YHR165C';'YHR166C';'YHR167W';'YHR168W';'YHR169W';'YHR170W';'YHR171W';'YHR172W';'YHR174W';'YHR175W';'YHR176W';'YHR178W';'YHR179W';'YHR181W';'YHR183W';'YHR184W';'YHR185C';'YHR186C';'YHR187W';'YHR188C';'YHR189W';'YHR190W';'YHR191C';'YHR192W';'YHR193C';'YHR194W';'YHR195W';'YHR196W';'YHR197W';'YHR198C';'YHR199C';'YHR199C-A';'YHR200W';'YHR201C';'YHR203C';'YHR204W';'YHR205W';'YHR206W';'YHR207C';'YHR208W';'YHR209W';'YHR211W';'YHR215W';'YHR216W';'YIL002C';'YIL003W';'YIL004C';'YIL005W';'YIL006W';'YIL007C';'YIL008W';'YIL009C-A';'YIL009W';'YIL010W';'YIL011W';'YIL013C';'YIL014W';'YIL015W';'YIL016W';'YIL017C';'YIL018W';'YIL019W';'YIL020C';'YIL021W';'YIL022W';'YIL023C';'YIL026C';'YIL027C';'YIL030C';'YIL031W';'YIL033C';'YIL034C';'YIL035C';'YIL036W';'YIL037C';'YIL038C';'YIL039W';'YIL040W';'YIL041W';'YIL042C';'YIL043C';'YIL044C';'YIL045W';'YIL046W';'YIL047C';'YIL048W';'YIL049W';'YIL050W';'YIL051C';'YIL052C';'YIL053W';'YIL056W';'YIL057C';'YIL061C';'YIL062C';'YIL063C';'YIL064W';'YIL065C';'YIL066C';'YIL068C';'YIL069C';'YIL070C';'YIL071C';'YIL072W';'YIL073C';'YIL074C';'YIL075C';'YIL076W';'YIL078W';'YIL079C';'YIL083C';'YIL084C';'YIL085C';'YIL087C';'YIL088C';'YIL089W';'YIL090W';'YIL091C';'YIL093C';'YIL094C';'YIL095W';'YIL097W';'YIL098C';'YIL099W';'YIL101C';'YIL103W';'YIL104C';'YIL105C';'YIL106W';'YIL107C';'YIL109C';'YIL110W';'YIL111W';'YIL112W';'YIL113W';'YIL114C';'YIL115C';'YIL116W';'YIL117C';'YIL118W';'YIL119C';'YIL120W';'YIL121W';'YIL122W';'YIL123W';'YIL124W';'YIL125W';'YIL126W';'YIL128W';'YIL129C';'YIL130W';'YIL131C';'YIL132C';'YIL133C';'YIL134W';'YIL135C';'YIL136W';'YIL137C';'YIL138C';'YIL139C';'YIL140W';'YIL142W';'YIL143C';'YIL144W';'YIL145C';'YIL146C';'YIL147C';'YIL148W';'YIL149C';'YIL150C';'YIL153W';'YIL154C';'YIL155C';'YIL156W';'YIL157C';'YIL158W';'YIL159W';'YIL160C';'YIL162W';'YIL164C';'YIL172C';'YIL173W';'YIR001C';'YIR002C';'YIR003W';'YIR004W';'YIR005W';'YIR006C';'YIR008C';'YIR009W';'YIR010W';'YIR011C';'YIR012W';'YIR013C';'YIR015W';'YIR017C';'YIR018W';'YIR019C';'YIR021W';'YIR022W';'YIR023W';'YIR024C';'YIR025W';'YIR026C';'YIR027C';'YIR028W';'YIR029W';'YIR030C';'YIR031C';'YIR032C';'YIR033W';'YIR034C';'YIR037W';'YIR038C';'YIR039C';'YIR041W';'YJL001W';'YJL002C';'YJL003W';'YJL004C';'YJL005W';'YJL006C';'YJL008C';'YJL010C';'YJL011C';'YJL012C';'YJL013C';'YJL014W';'YJL019W';'YJL020C';'YJL023C';'YJL024C';'YJL025W';'YJL026W';'YJL028W';'YJL029C';'YJL030W';'YJL031C';'YJL033W';'YJL034W';'YJL035C';'YJL036W';'YJL037W';'YJL038C';'YJL039C';'YJL041W';'YJL042W';'YJL044C';'YJL045W';'YJL046W';'YJL047C';'YJL048C';'YJL050W';'YJL051W';'YJL052W';'YJL053W';'YJL054W';'YJL056C';'YJL057C';'YJL058C';'YJL059W';'YJL060W';'YJL061W';'YJL062W';'YJL062W-A';'YJL063C';'YJL065C';'YJL066C';'YJL068C';'YJL069C';'YJL071W';'YJL072C';'YJL073W';'YJL074C';'YJL076W';'YJL077C';'YJL078C';'YJL079C';'YJL080C';'YJL081C';'YJL082W';'YJL083W';'YJL084C';'YJL085W';'YJL087C';'YJL088W';'YJL089W';'YJL090C';'YJL091C';'YJL092W';'YJL093C';'YJL094C';'YJL095W';'YJL096W';'YJL097W';'YJL098W';'YJL099W';'YJL100W';'YJL101C';'YJL102W';'YJL103C';'YJL104W';'YJL105W';'YJL106W';'YJL108C';'YJL109C';'YJL110C';'YJL111W';'YJL112W';'YJL115W';'YJL116C';'YJL117W';'YJL118W';'YJL121C';'YJL122W';'YJL123C';'YJL124C';'YJL125C';'YJL126W';'YJL127C';'YJL128C';'YJL129C';'YJL130C';'YJL131C';'YJL133W';'YJL134W';'YJL136C';'YJL137C';'YJL138C';'YJL139C';'YJL140W';'YJL141C';'YJL143W';'YJL144W';'YJL145W';'YJL146W';'YJL148W';'YJL149W';'YJL151C';'YJL153C';'YJL154C';'YJL155C';'YJL156C';'YJL157C';'YJL158C';'YJL159W';'YJL162C';'YJL164C';'YJL165C';'YJL166W';'YJL167W';'YJL168C';'YJL170C';'YJL171C';'YJL172W';'YJL173C';'YJL174W';'YJL176C';'YJL177W';'YJL178C';'YJL179W';'YJL180C';'YJL183W';'YJL184W';'YJL185C';'YJL186W';'YJL187C';'YJL189W';'YJL190C';'YJL191W';'YJL192C';'YJL194W';'YJL196C';'YJL197W';'YJL198W';'YJL200C';'YJL201W';'YJL203W';'YJL204C';'YJL205C';'YJL207C';'YJL208C';'YJL209W';'YJL210W';'YJL212C';'YJL213W';'YJL214W';'YJL216C';'YJL217W';'YJL219W';'YJL221C';'YJL222W';'YJL223C';'YJR001W';'YJR002W';'YJR004C';'YJR005W';'YJR006W';'YJR007W';'YJR008W';'YJR009C';'YJR010C-A';'YJR010W';'YJR013W';'YJR014W';'YJR016C';'YJR017C';'YJR019C';'YJR021C';'YJR022W';'YJR024C';'YJR025C';'YJR031C';'YJR032W';'YJR033C';'YJR034W';'YJR035W';'YJR036C';'YJR040W';'YJR041C';'YJR042W';'YJR043C';'YJR044C';'YJR045C';'YJR046W';'YJR047C';'YJR048W';'YJR049C';'YJR050W';'YJR051W';'YJR052W';'YJR053W';'YJR054W';'YJR055W';'YJR057W';'YJR058C';'YJR059W';'YJR060W';'YJR062C';'YJR063W';'YJR064W';'YJR065C';'YJR066W';'YJR067C';'YJR068W';'YJR069C';'YJR070C';'YJR072C';'YJR073C';'YJR074W';'YJR075W';'YJR076C';'YJR077C';'YJR078W';'YJR080C';'YJR082C';'YJR083C';'YJR084W';'YJR086W';'YJR088C';'YJR089W';'YJR090C';'YJR091C';'YJR092W';'YJR093C';'YJR094C';'YJR094W-A';'YJR095W';'YJR096W';'YJR097W';'YJR099W';'YJR100C';'YJR101W';'YJR102C';'YJR103W';'YJR104C';'YJR105W';'YJR106W';'YJR108W';'YJR109C';'YJR110W';'YJR112W';'YJR113C';'YJR117W';'YJR118C';'YJR119C';'YJR120W';'YJR121W';'YJR122W';'YJR123W';'YJR125C';'YJR126C';'YJR127C';'YJR130C';'YJR131W';'YJR132W';'YJR133W';'YJR134C';'YJR135C';'YJR135W-A';'YJR136C';'YJR137C';'YJR138W';'YJR139C';'YJR140C';'YJR143C';'YJR144W';'YJR145C';'YJR147W';'YJR148W';'YJR150C';'YJR151C';'YJR152W';'YJR153W';'YJR155W';'YJR156C';'YJR158W';'YJR159W';'YJR160C';'YJR161C';'YKL001C';'YKL002W';'YKL003C';'YKL004W';'YKL005C';'YKL006C-A';'YKL006W';'YKL007W';'YKL008C';'YKL009W';'YKL010C';'YKL011C';'YKL012W';'YKL013C';'YKL014C';'YKL015W';'YKL016C';'YKL017C';'YKL018W';'YKL019W';'YKL020C';'YKL021C';'YKL022C';'YKL024C';'YKL025C';'YKL026C';'YKL027W';'YKL028W';'YKL029C';'YKL032C';'YKL033W';'YKL034W';'YKL035W';'YKL037W';'YKL038W';'YKL039W';'YKL040C';'YKL041W';'YKL042W';'YKL043W';'YKL045W';'YKL046C';'YKL048C';'YKL049C';'YKL050C';'YKL051W';'YKL052C';'YKL053C-A';'YKL054C';'YKL055C';'YKL056C';'YKL057C';'YKL058W';'YKL059C';'YKL060C';'YKL062W';'YKL064W';'YKL065C';'YKL067W';'YKL068W';'YKL069W';'YKL072W';'YKL073W';'YKL074C';'YKL078W';'YKL079W';'YKL080W';'YKL081W';'YKL082C';'YKL084W';'YKL085W';'YKL086W';'YKL087C';'YKL088W';'YKL089W';'YKL090W';'YKL091C';'YKL092C';'YKL093W';'YKL094W';'YKL095W';'YKL096W';'YKL096W-A';'YKL098W';'YKL099C';'YKL101W';'YKL103C';'YKL104C';'YKL105C';'YKL106W';'YKL108W';'YKL109W';'YKL110C';'YKL112W';'YKL113C';'YKL114C';'YKL116C';'YKL117W';'YKL119C';'YKL120W';'YKL122C';'YKL124W';'YKL125W';'YKL126W';'YKL127W';'YKL128C';'YKL129C';'YKL130C';'YKL132C';'YKL134C';'YKL135C';'YKL137W';'YKL138C';'YKL138C-A';'YKL139W';'YKL140W';'YKL141W';'YKL142W';'YKL143W';'YKL144C';'YKL145W';'YKL146W';'YKL148C';'YKL149C';'YKL150W';'YKL151C';'YKL152C';'YKL154W';'YKL155C';'YKL156W';'YKL157W';'YKL159C';'YKL160W';'YKL161C';'YKL163W';'YKL164C';'YKL165C';'YKL166C';'YKL167C';'YKL168C';'YKL170W';'YKL171W';'YKL172W';'YKL173W';'YKL174C';'YKL175W';'YKL176C';'YKL178C';'YKL179C';'YKL180W';'YKL181W';'YKL182W';'YKL183W';'YKL184W';'YKL185W';'YKL186C';'YKL188C';'YKL189W';'YKL190W';'YKL191W';'YKL192C';'YKL193C';'YKL194C';'YKL195W';'YKL196C';'YKL197C';'YKL198C';'YKL201C';'YKL203C';'YKL204W';'YKL205W';'YKL206C';'YKL207W';'YKL208W';'YKL209C';'YKL210W';'YKL211C';'YKL212W';'YKL213C';'YKL214C';'YKL215C';'YKL216W';'YKL217W';'YKL218C';'YKL219W';'YKL220C';'YKL221W';'YKL222C';'YKL224C';'YKR001C';'YKR002W';'YKR003W';'YKR004C';'YKR006C';'YKR007W';'YKR008W';'YKR009C';'YKR010C';'YKR013W';'YKR014C';'YKR016W';'YKR017C';'YKR019C';'YKR020W';'YKR021W';'YKR022C';'YKR024C';'YKR025W';'YKR026C';'YKR027W';'YKR028W';'YKR029C';'YKR030W';'YKR031C';'YKR034W';'YKR035W-A';'YKR036C';'YKR037C';'YKR038C';'YKR039W';'YKR041W';'YKR042W';'YKR043C';'YKR044W';'YKR046C';'YKR048C';'YKR049C';'YKR050W';'YKR052C';'YKR053C';'YKR054C';'YKR055W';'YKR056W';'YKR057W';'YKR058W';'YKR059W';'YKR060W';'YKR061W';'YKR062W';'YKR063C';'YKR064W';'YKR065C';'YKR066C';'YKR067W';'YKR068C';'YKR069W';'YKR071C';'YKR072C';'YKR074W';'YKR076W';'YKR077W';'YKR078W';'YKR079C';'YKR080W';'YKR081C';'YKR082W';'YKR083C';'YKR084C';'YKR085C';'YKR086W';'YKR087C';'YKR088C';'YKR089C';'YKR090W';'YKR091W';'YKR092C';'YKR093W';'YKR094C';'YKR095W';'YKR095W-A';'YKR096W';'YKR097W';'YKR098C';'YKR099W';'YKR100C';'YKR101W';'YKR102W';'YKR103W';'YKR104W';'YKR106W';'YLL001W';'YLL002W';'YLL003W';'YLL004W';'YLL005C';'YLL006W';'YLL008W';'YLL009C';'YLL010C';'YLL011W';'YLL012W';'YLL013C';'YLL014W';'YLL015W';'YLL018C';'YLL018C-A';'YLL019C';'YLL021W';'YLL022C';'YLL023C';'YLL024C';'YLL025W';'YLL026W';'YLL027W';'YLL028W';'YLL029W';'YLL031C';'YLL032C';'YLL033W';'YLL034C';'YLL035W';'YLL036C';'YLL038C';'YLL039C';'YLL040C';'YLL041C';'YLL042C';'YLL043W';'YLL045C';'YLL046C';'YLL048C';'YLL049W';'YLL050C';'YLL051C';'YLL052C';'YLL055W';'YLL057C';'YLL060C';'YLL061W';'YLL062C';'YLL063C';'YLL064C';'YLR002C';'YLR003C';'YLR004C';'YLR005W';'YLR006C';'YLR007W';'YLR008C';'YLR009W';'YLR010C';'YLR011W';'YLR013W';'YLR014C';'YLR015W';'YLR016C';'YLR017W';'YLR018C';'YLR019W';'YLR020C';'YLR021W';'YLR022C';'YLR023C';'YLR024C';'YLR025W';'YLR026C';'YLR027C';'YLR028C';'YLR029C';'YLR032W';'YLR033W';'YLR034C';'YLR035C';'YLR037C';'YLR038C';'YLR039C';'YLR043C';'YLR044C';'YLR045C';'YLR047C';'YLR048W';'YLR051C';'YLR052W';'YLR054C';'YLR055C';'YLR056W';'YLR057W';'YLR058C';'YLR059C';'YLR060W';'YLR061W';'YLR064W';'YLR065C';'YLR066W';'YLR067C';'YLR068W';'YLR069C';'YLR070C';'YLR071C';'YLR073C';'YLR074C';'YLR075W';'YLR077W';'YLR078C';'YLR079W';'YLR080W';'YLR081W';'YLR082C';'YLR083C';'YLR084C';'YLR085C';'YLR086W';'YLR087C';'YLR088W';'YLR089C';'YLR090W';'YLR091W';'YLR092W';'YLR093C';'YLR094C';'YLR095C';'YLR096W';'YLR097C';'YLR098C';'YLR099C';'YLR099W-A';'YLR100W';'YLR102C';'YLR103C';'YLR105C';'YLR106C';'YLR107W';'YLR109W';'YLR110C';'YLR113W';'YLR114C';'YLR115W';'YLR116W';'YLR117C';'YLR118C';'YLR119W';'YLR120C';'YLR121C';'YLR127C';'YLR128W';'YLR129W';'YLR130C';'YLR131C';'YLR132C';'YLR133W';'YLR134W';'YLR135W';'YLR136C';'YLR137W';'YLR138W';'YLR139C';'YLR141W';'YLR142W';'YLR144C';'YLR145W';'YLR146C';'YLR147C';'YLR148W';'YLR150W';'YLR151C';'YLR153C';'YLR154C';'YLR154W-C';'YLR155C';'YLR157C';'YLR158C';'YLR160C';'YLR162W';'YLR163C';'YLR164W';'YLR165C';'YLR166C';'YLR167W';'YLR168C';'YLR170C';'YLR172C';'YLR174W';'YLR175W';'YLR176C';'YLR178C';'YLR179C';'YLR180W';'YLR181C';'YLR182W';'YLR183C';'YLR185W';'YLR186W';'YLR188W';'YLR189C';'YLR190W';'YLR191W';'YLR192C';'YLR193C';'YLR194C';'YLR195C';'YLR196W';'YLR197W';'YLR199C';'YLR200W';'YLR201C';'YLR203C';'YLR204W';'YLR205C';'YLR206W';'YLR207W';'YLR208W';'YLR209C';'YLR210W';'YLR212C';'YLR213C';'YLR214W';'YLR215C';'YLR216C';'YLR218C';'YLR219W';'YLR220W';'YLR221C';'YLR222C';'YLR223C';'YLR226W';'YLR227C';'YLR228C';'YLR229C';'YLR231C';'YLR233C';'YLR234W';'YLR237W';'YLR238W';'YLR239C';'YLR240W';'YLR242C';'YLR243W';'YLR244C';'YLR245C';'YLR246W';'YLR247C';'YLR248W';'YLR249W';'YLR250W';'YLR251W';'YLR254C';'YLR256W';'YLR258W';'YLR259C';'YLR260W';'YLR262C';'YLR262C-A';'YLR263W';'YLR264W';'YLR265C';'YLR266C';'YLR268W';'YLR270W';'YLR272C';'YLR273C';'YLR274W';'YLR275W';'YLR276C';'YLR277C';'YLR284C';'YLR285W';'YLR286C';'YLR287C-A';'YLR288C';'YLR289W';'YLR291C';'YLR292C';'YLR293C';'YLR295C';'YLR298C';'YLR299W';'YLR300W';'YLR301W';'YLR303W';'YLR304C';'YLR305C';'YLR306W';'YLR307W';'YLR308W';'YLR309C';'YLR310C';'YLR312W-A';'YLR313C';'YLR314C';'YLR315W';'YLR316C';'YLR318W';'YLR319C';'YLR320W';'YLR321C';'YLR323C';'YLR324W';'YLR325C';'YLR327C';'YLR328W';'YLR329W';'YLR330W';'YLR332W';'YLR333C';'YLR335W';'YLR336C';'YLR337C';'YLR340W';'YLR341W';'YLR342W';'YLR343W';'YLR344W';'YLR347C';'YLR348C';'YLR350W';'YLR351C';'YLR353W';'YLR354C';'YLR355C';'YLR356W';'YLR357W';'YLR359W';'YLR360W';'YLR361C';'YLR362W';'YLR363C';'YLR364W';'YLR367W';'YLR368W';'YLR369W';'YLR370C';'YLR371W';'YLR372W';'YLR373C';'YLR375W';'YLR376C';'YLR377C';'YLR378C';'YLR380W';'YLR381W';'YLR382C';'YLR383W';'YLR384C';'YLR385C';'YLR386W';'YLR387C';'YLR388W';'YLR389C';'YLR390W';'YLR390W-A';'YLR392C';'YLR393W';'YLR394W';'YLR395C';'YLR396C';'YLR397C';'YLR398C';'YLR399C';'YLR401C';'YLR403W';'YLR404W';'YLR405W';'YLR406C';'YLR409C';'YLR410W';'YLR411W';'YLR412W';'YLR414C';'YLR417W';'YLR418C';'YLR420W';'YLR421C';'YLR423C';'YLR424W';'YLR425W';'YLR427W';'YLR429W';'YLR430W';'YLR431C';'YLR432W';'YLR433C';'YLR435W';'YLR436C';'YLR437C';'YLR438C-A';'YLR438W';'YLR439W';'YLR440C';'YLR441C';'YLR442C';'YLR443W';'YLR445W';'YLR447C';'YLR448W';'YLR449W';'YLR450W';'YLR451W';'YLR452C';'YLR453C';'YLR457C';'YLR459W';'YLR461W';'YLR466W';'YLR467W';'YML001W';'YML004C';'YML005W';'YML006C';'YML007W';'YML008C';'YML009C';'YML010W';'YML011C';'YML012W';'YML013W';'YML014W';'YML015C';'YML016C';'YML017W';'YML019W';'YML021C';'YML022W';'YML023C';'YML024W';'YML025C';'YML026C';'YML027W';'YML028W';'YML029W';'YML030W';'YML031W';'YML032C';'YML034W';'YML035C';'YML036W';'YML038C';'YML041C';'YML042W';'YML043C';'YML046W';'YML047C';'YML048W';'YML049C';'YML050W';'YML051W';'YML052W';'YML054C';'YML055W';'YML056C';'YML057W';'YML058W';'YML058W-A';'YML059C';'YML060W';'YML061C';'YML062C';'YML063W';'YML064C';'YML065W';'YML066C';'YML067C';'YML068W';'YML069W';'YML070W';'YML071C';'YML072C';'YML073C';'YML074C';'YML075C';'YML076C';'YML077W';'YML078W';'YML080W';'YML081C-A';'YML081W';'YML085C';'YML086C';'YML087C';'YML088W';'YML091C';'YML092C';'YML093W';'YML094W';'YML095C';'YML097C';'YML098W';'YML099C';'YML100W';'YML101C';'YML102W';'YML103C';'YML104C';'YML105C';'YML106W';'YML107C';'YML109W';'YML110C';'YML111W';'YML112W';'YML113W';'YML114C';'YML115C';'YML116W';'YML117W';'YML118W';'YML120C';'YML121W';'YML123C';'YML124C';'YML125C';'YML126C';'YML127W';'YML128C';'YML129C';'YML130C';'YML132W';'YMR001C';'YMR002W';'YMR003W';'YMR004W';'YMR005W';'YMR006C';'YMR008C';'YMR009W';'YMR011W';'YMR012W';'YMR013C';'YMR014W';'YMR015C';'YMR016C';'YMR017W';'YMR019W';'YMR020W';'YMR021C';'YMR022W';'YMR023C';'YMR024W';'YMR025W';'YMR026C';'YMR028W';'YMR029C';'YMR030W';'YMR031C';'YMR032W';'YMR033W';'YMR035W';'YMR036C';'YMR037C';'YMR038C';'YMR039C';'YMR040W';'YMR041C';'YMR042W';'YMR043W';'YMR044W';'YMR047C';'YMR048W';'YMR049C';'YMR052W';'YMR053C';'YMR054W';'YMR055C';'YMR056C';'YMR058W';'YMR059W';'YMR060C';'YMR061W';'YMR062C';'YMR063W';'YMR064W';'YMR065W';'YMR066W';'YMR067C';'YMR068W';'YMR069W';'YMR070W';'YMR071C';'YMR072W';'YMR073C';'YMR074C';'YMR075W';'YMR076C';'YMR077C';'YMR078C';'YMR079W';'YMR080C';'YMR081C';'YMR083W';'YMR086W';'YMR087W';'YMR088C';'YMR089C';'YMR091C';'YMR092C';'YMR093W';'YMR094W';'YMR095C';'YMR096W';'YMR097C';'YMR098C';'YMR099C';'YMR100W';'YMR101C';'YMR104C';'YMR105C';'YMR106C';'YMR107W';'YMR108W';'YMR109W';'YMR110C';'YMR112C';'YMR113W';'YMR114C';'YMR115W';'YMR116C';'YMR117C';'YMR119W';'YMR120C';'YMR121C';'YMR123W';'YMR125W';'YMR127C';'YMR128W';'YMR129W';'YMR131C';'YMR133W';'YMR135C';'YMR136W';'YMR137C';'YMR138W';'YMR139W';'YMR140W';'YMR142C';'YMR143W';'YMR145C';'YMR146C';'YMR148W';'YMR149W';'YMR150C';'YMR152W';'YMR153W';'YMR154C';'YMR156C';'YMR157C';'YMR158W';'YMR159C';'YMR161W';'YMR162C';'YMR163C';'YMR164C';'YMR165C';'YMR167W';'YMR168C';'YMR169C';'YMR170C';'YMR171C';'YMR172W';'YMR173W';'YMR174C';'YMR175W';'YMR176W';'YMR177W';'YMR179W';'YMR180C';'YMR182C';'YMR183C';'YMR184W';'YMR186W';'YMR188C';'YMR189W';'YMR190C';'YMR191W';'YMR192W';'YMR193W';'YMR194C-B';'YMR194W';'YMR195W';'YMR197C';'YMR198W';'YMR199W';'YMR200W';'YMR201C';'YMR202W';'YMR203W';'YMR204C';'YMR205C';'YMR207C';'YMR208W';'YMR210W';'YMR211W';'YMR212C';'YMR213W';'YMR214W';'YMR215W';'YMR216C';'YMR217W';'YMR218C';'YMR219W';'YMR220W';'YMR222C';'YMR223W';'YMR224C';'YMR225C';'YMR226C';'YMR227C';'YMR228W';'YMR229C';'YMR230W';'YMR231W';'YMR232W';'YMR233W';'YMR234W';'YMR235C';'YMR236W';'YMR237W';'YMR238W';'YMR239C';'YMR240C';'YMR241W';'YMR242C';'YMR243C';'YMR244C-A';'YMR246W';'YMR247C';'YMR250W';'YMR251W';'YMR251W-A';'YMR255W';'YMR256C';'YMR257C';'YMR258C';'YMR259C';'YMR260C';'YMR261C';'YMR263W';'YMR264W';'YMR266W';'YMR267W';'YMR268C';'YMR269W';'YMR270C';'YMR271C';'YMR272C';'YMR273C';'YMR274C';'YMR275C';'YMR276W';'YMR277W';'YMR278W';'YMR279C';'YMR280C';'YMR281W';'YMR282C';'YMR283C';'YMR284W';'YMR285C';'YMR286W';'YMR287C';'YMR288W';'YMR289W';'YMR290C';'YMR291W';'YMR292W';'YMR293C';'YMR294W';'YMR295C';'YMR296C';'YMR297W';'YMR298W';'YMR299C';'YMR300C';'YMR301C';'YMR302C';'YMR303C';'YMR304W';'YMR305C';'YMR306W';'YMR307W';'YMR308C';'YMR309C';'YMR311C';'YMR312W';'YMR313C';'YMR314W';'YMR315W';'YMR316W';'YMR318C';'YMR319C';'YMR323W';'YMR325W';'YNL001W';'YNL002C';'YNL003C';'YNL004W';'YNL005C';'YNL006W';'YNL007C';'YNL008C';'YNL009W';'YNL012W';'YNL014W';'YNL015W';'YNL016W';'YNL020C';'YNL021W';'YNL023C';'YNL024C-A';'YNL025C';'YNL026W';'YNL027W';'YNL029C';'YNL030W';'YNL031C';'YNL032W';'YNL036W';'YNL037C';'YNL038W';'YNL039W';'YNL041C';'YNL042W';'YNL044W';'YNL045W';'YNL047C';'YNL048W';'YNL049C';'YNL051W';'YNL052W';'YNL053W';'YNL054W';'YNL055C';'YNL056W';'YNL059C';'YNL061W';'YNL062C';'YNL063W';'YNL064C';'YNL065W';'YNL066W';'YNL067W';'YNL068C';'YNL069C';'YNL070W';'YNL071W';'YNL072W';'YNL073W';'YNL074C';'YNL075W';'YNL076W';'YNL077W';'YNL078W';'YNL079C';'YNL080C';'YNL081C';'YNL082W';'YNL083W';'YNL084C';'YNL085W';'YNL087W';'YNL088W';'YNL090W';'YNL091W';'YNL093W';'YNL094W';'YNL096C';'YNL097C';'YNL098C';'YNL099C';'YNL100W';'YNL101W';'YNL102W';'YNL103W';'YNL104C';'YNL106C';'YNL107W';'YNL110C';'YNL111C';'YNL112W';'YNL113W';'YNL116W';'YNL117W';'YNL118C';'YNL119W';'YNL121C';'YNL123W';'YNL124W';'YNL125C';'YNL126W';'YNL127W';'YNL128W';'YNL129W';'YNL130C';'YNL131W';'YNL132W';'YNL133C';'YNL135C';'YNL136W';'YNL137C';'YNL138W';'YNL138W-A';'YNL139C';'YNL141W';'YNL142W';'YNL145W';'YNL147W';'YNL148C';'YNL149C';'YNL151C';'YNL152W';'YNL153C';'YNL154C';'YNL156C';'YNL157W';'YNL158W';'YNL159C';'YNL160W';'YNL161W';'YNL162W';'YNL163C';'YNL164C';'YNL166C';'YNL167C';'YNL169C';'YNL172W';'YNL173C';'YNL175C';'YNL177C';'YNL178W';'YNL180C';'YNL182C';'YNL183C';'YNL185C';'YNL186W';'YNL187W';'YNL188W';'YNL189W';'YNL191W';'YNL192W';'YNL194C';'YNL197C';'YNL199C';'YNL200C';'YNL201C';'YNL202W';'YNL204C';'YNL206C';'YNL207W';'YNL208W';'YNL209W';'YNL210W';'YNL212W';'YNL213C';'YNL214W';'YNL215W';'YNL216W';'YNL218W';'YNL219C';'YNL220W';'YNL221C';'YNL222W';'YNL223W';'YNL224C';'YNL225C';'YNL227C';'YNL229C';'YNL230C';'YNL231C';'YNL232W';'YNL233W';'YNL234W';'YNL236W';'YNL237W';'YNL238W';'YNL239W';'YNL240C';'YNL241C';'YNL242W';'YNL243W';'YNL244C';'YNL245C';'YNL246W';'YNL247W';'YNL248C';'YNL249C';'YNL250W';'YNL251C';'YNL252C';'YNL253W';'YNL254C';'YNL255C';'YNL256W';'YNL257C';'YNL258C';'YNL259C';'YNL260C';'YNL261W';'YNL262W';'YNL263C';'YNL264C';'YNL265C';'YNL267W';'YNL268W';'YNL269W';'YNL270C';'YNL271C';'YNL272C';'YNL273W';'YNL274C';'YNL275W';'YNL277W';'YNL278W';'YNL279W';'YNL280C';'YNL281W';'YNL282W';'YNL283C';'YNL284C';'YNL286W';'YNL287W';'YNL288W';'YNL289W';'YNL290W';'YNL291C';'YNL292W';'YNL293W';'YNL294C';'YNL297C';'YNL298W';'YNL299W';'YNL301C';'YNL302C';'YNL304W';'YNL305C';'YNL306W';'YNL307C';'YNL308C';'YNL309W';'YNL310C';'YNL311C';'YNL312W';'YNL313C';'YNL314W';'YNL315C';'YNL316C';'YNL317W';'YNL318C';'YNL321W';'YNL322C';'YNL323W';'YNL325C';'YNL326C';'YNL327W';'YNL328C';'YNL329C';'YNL330C';'YNL331C';'YNL332W';'YNL333W';'YNL334C';'YNL336W';'YNL339C';'YNR001C';'YNR002C';'YNR003C';'YNR006W';'YNR007C';'YNR008W';'YNR009W';'YNR010W';'YNR011C';'YNR012W';'YNR013C';'YNR015W';'YNR016C';'YNR017W';'YNR018W';'YNR019W';'YNR020C';'YNR022C';'YNR023W';'YNR024W';'YNR026C';'YNR027W';'YNR028W';'YNR030W';'YNR031C';'YNR032C-A';'YNR032W';'YNR033W';'YNR034W';'YNR035C';'YNR036C';'YNR037C';'YNR038W';'YNR039C';'YNR041C';'YNR043W';'YNR044W';'YNR045W';'YNR046W';'YNR047W';'YNR048W';'YNR049C';'YNR050C';'YNR051C';'YNR052C';'YNR053C';'YNR054C';'YNR055C';'YNR056C';'YNR057C';'YNR058W';'YNR059W';'YNR060W';'YNR064C';'YNR067C';'YNR069C';'YNR072W';'YNR074C';'YNR075W';'YNR076W';'YOL001W';'YOL002C';'YOL003C';'YOL004W';'YOL005C';'YOL006C';'YOL007C';'YOL008W';'YOL009C';'YOL010W';'YOL011W';'YOL012C';'YOL013C';'YOL015W';'YOL016C';'YOL017W';'YOL018C';'YOL020W';'YOL021C';'YOL022C';'YOL023W';'YOL025W';'YOL026C';'YOL027C';'YOL028C';'YOL030W';'YOL031C';'YOL032W';'YOL033W';'YOL034W';'YOL038W';'YOL039W';'YOL040C';'YOL041C';'YOL042W';'YOL043C';'YOL044W';'YOL045W';'YOL049W';'YOL051W';'YOL052C';'YOL052C-A';'YOL053W';'YOL054W';'YOL055C';'YOL056W';'YOL057W';'YOL058W';'YOL059W';'YOL060C';'YOL061W';'YOL062C';'YOL063C';'YOL064C';'YOL065C';'YOL066C';'YOL067C';'YOL068C';'YOL069W';'YOL070C';'YOL071W';'YOL072W';'YOL076W';'YOL077C';'YOL077W-A';'YOL078W';'YOL080C';'YOL081W';'YOL082W';'YOL083W';'YOL084W';'YOL086C';'YOL086W-A';'YOL087C';'YOL088C';'YOL089C';'YOL090W';'YOL091W';'YOL092W';'YOL093W';'YOL094C';'YOL095C';'YOL096C';'YOL097C';'YOL100W';'YOL101C';'YOL102C';'YOL103W';'YOL104C';'YOL105C';'YOL108C';'YOL109W';'YOL110W';'YOL111C';'YOL112W';'YOL113W';'YOL115W';'YOL116W';'YOL117W';'YOL119C';'YOL120C';'YOL121C';'YOL122C';'YOL123W';'YOL124C';'YOL125W';'YOL126C';'YOL127W';'YOL128C';'YOL129W';'YOL130W';'YOL132W';'YOL133W';'YOL135C';'YOL136C';'YOL137W';'YOL138C';'YOL139C';'YOL140W';'YOL141W';'YOL142W';'YOL143C';'YOL144W';'YOL145C';'YOL146W';'YOL147C';'YOL148C';'YOL149W';'YOL151W';'YOL152W';'YOL154W';'YOL155C';'YOL156W';'YOL157C';'YOL158C';'YOL159C';'YOL159C-A';'YOL161C';'YOL164W';'YOL165C';'YOR001W';'YOR002W';'YOR003W';'YOR004W';'YOR005C';'YOR006C';'YOR007C';'YOR008C';'YOR009W';'YOR010C';'YOR011W';'YOR014W';'YOR016C';'YOR017W';'YOR018W';'YOR019W';'YOR020C';'YOR021C';'YOR023C';'YOR025W';'YOR026W';'YOR027W';'YOR028C';'YOR030W';'YOR031W';'YOR032C';'YOR033C';'YOR034C';'YOR035C';'YOR036W';'YOR037W';'YOR038C';'YOR039W';'YOR040W';'YOR042W';'YOR043W';'YOR044W';'YOR045W';'YOR046C';'YOR047C';'YOR048C';'YOR049C';'YOR051C';'YOR052C';'YOR054C';'YOR056C';'YOR057W';'YOR058C';'YOR060C';'YOR061W';'YOR063W';'YOR064C';'YOR065W';'YOR066W';'YOR067C';'YOR068C';'YOR069W';'YOR070C';'YOR071C';'YOR073W';'YOR074C';'YOR075W';'YOR076C';'YOR077W';'YOR078W';'YOR079C';'YOR080W';'YOR081C';'YOR083W';'YOR084W';'YOR085W';'YOR086C';'YOR087W';'YOR089C';'YOR090C';'YOR091W';'YOR092W';'YOR094W';'YOR095C';'YOR096W';'YOR098C';'YOR099W';'YOR100C';'YOR101W';'YOR103C';'YOR104W';'YOR106W';'YOR107W';'YOR108W';'YOR109W';'YOR110W';'YOR112W';'YOR113W';'YOR115C';'YOR116C';'YOR117W';'YOR118W';'YOR119C';'YOR120W';'YOR122C';'YOR123C';'YOR124C';'YOR125C';'YOR126C';'YOR127W';'YOR128C';'YOR129C';'YOR130C';'YOR132W';'YOR133W';'YOR134W';'YOR136W';'YOR137C';'YOR138C';'YOR140W';'YOR141C';'YOR142W';'YOR143C';'YOR144C';'YOR145C';'YOR147W';'YOR148C';'YOR149C';'YOR150W';'YOR151C';'YOR153W';'YOR155C';'YOR156C';'YOR157C';'YOR158W';'YOR159C';'YOR160W';'YOR161C';'YOR162C';'YOR163W';'YOR164C';'YOR165W';'YOR166C';'YOR167C';'YOR168W';'YOR171C';'YOR172W';'YOR173W';'YOR174W';'YOR175C';'YOR176W';'YOR177C';'YOR178C';'YOR179C';'YOR180C';'YOR181W';'YOR182C';'YOR184W';'YOR185C';'YOR187W';'YOR188W';'YOR189W';'YOR190W';'YOR191W';'YOR192C';'YOR193W';'YOR194C';'YOR195W';'YOR196C';'YOR197W';'YOR198C';'YOR201C';'YOR202W';'YOR204W';'YOR205C';'YOR206W';'YOR207C';'YOR208W';'YOR209C';'YOR210W';'YOR211C';'YOR212W';'YOR213C';'YOR215C';'YOR216C';'YOR217W';'YOR219C';'YOR221C';'YOR222W';'YOR223W';'YOR224C';'YOR226C';'YOR227W';'YOR228C';'YOR229W';'YOR230W';'YOR231W';'YOR232W';'YOR233W';'YOR234C';'YOR236W';'YOR237W';'YOR239W';'YOR241W';'YOR242C';'YOR243C';'YOR244W';'YOR245C';'YOR246C';'YOR247W';'YOR249C';'YOR250C';'YOR251C';'YOR252W';'YOR253W';'YOR254C';'YOR255W';'YOR256C';'YOR257W';'YOR258W';'YOR259C';'YOR260W';'YOR261C';'YOR262W';'YOR264W';'YOR265W';'YOR266W';'YOR267C';'YOR269W';'YOR270C';'YOR272W';'YOR273C';'YOR274W';'YOR275C';'YOR276W';'YOR278W';'YOR279C';'YOR280C';'YOR281C';'YOR283W';'YOR284W';'YOR285W';'YOR286W';'YOR287C';'YOR288C';'YOR290C';'YOR291W';'YOR293W';'YOR294W';'YOR295W';'YOR297C';'YOR298C-A';'YOR298W';'YOR299W';'YOR301W';'YOR302W';'YOR303W';'YOR304W';'YOR305W';'YOR306C';'YOR307C';'YOR308C';'YOR310C';'YOR311C';'YOR312C';'YOR313C';'YOR315W';'YOR316C';'YOR317W';'YOR319W';'YOR320C';'YOR321W';'YOR322C';'YOR323C';'YOR324C';'YOR326W';'YOR327C';'YOR328W';'YOR329C';'YOR330C';'YOR332W';'YOR334W';'YOR335C';'YOR336W';'YOR337W';'YOR339C';'YOR340C';'YOR341W';'YOR344C';'YOR346W';'YOR347C';'YOR348C';'YOR349W';'YOR350C';'YOR351C';'YOR352W';'YOR353C';'YOR354C';'YOR355W';'YOR356W';'YOR357C';'YOR358W';'YOR359W';'YOR360C';'YOR361C';'YOR362C';'YOR363C';'YOR367W';'YOR368W';'YOR369C';'YOR370C';'YOR371C';'YOR372C';'YOR373W';'YOR374W';'YOR375C';'YOR377W';'YOR380W';'YOR381W';'YOR382W';'YOR383C';'YOR384W';'YOR386W';'YOR388C';'YOR391C';'YOR393W';'YOR394W';'YPL001W';'YPL002C';'YPL003W';'YPL004C';'YPL005W';'YPL006W';'YPL007C';'YPL008W';'YPL009C';'YPL010W';'YPL011C';'YPL012W';'YPL013C';'YPL015C';'YPL016W';'YPL017C';'YPL018W';'YPL019C';'YPL020C';'YPL021W';'YPL022W';'YPL023C';'YPL024W';'YPL026C';'YPL027W';'YPL028W';'YPL029W';'YPL030W';'YPL031C';'YPL032C';'YPL033C';'YPL036W';'YPL037C';'YPL038W';'YPL040C';'YPL042C';'YPL043W';'YPL045W';'YPL046C';'YPL047W';'YPL048W';'YPL049C';'YPL050C';'YPL051W';'YPL052W';'YPL053C';'YPL054W';'YPL055C';'YPL057C';'YPL058C';'YPL059W';'YPL060W';'YPL061W';'YPL063W';'YPL064C';'YPL065W';'YPL066W';'YPL069C';'YPL070W';'YPL072W';'YPL074W';'YPL075W';'YPL076W';'YPL078C';'YPL079W';'YPL081W';'YPL082C';'YPL083C';'YPL084W';'YPL085W';'YPL086C';'YPL087W';'YPL089C';'YPL090C';'YPL091W';'YPL092W';'YPL093W';'YPL094C';'YPL095C';'YPL096C-A';'YPL096W';'YPL097W';'YPL098C';'YPL099C';'YPL100W';'YPL101W';'YPL103C';'YPL104W';'YPL105C';'YPL106C';'YPL110C';'YPL111W';'YPL112C';'YPL115C';'YPL116W';'YPL117C';'YPL118W';'YPL119C';'YPL120W';'YPL121C';'YPL122C';'YPL123C';'YPL124W';'YPL125W';'YPL126W';'YPL127C';'YPL128C';'YPL129W';'YPL130W';'YPL131W';'YPL132W';'YPL133C';'YPL134C';'YPL135W';'YPL137C';'YPL138C';'YPL139C';'YPL140C';'YPL141C';'YPL143W';'YPL144W';'YPL145C';'YPL146C';'YPL147W';'YPL148C';'YPL149W';'YPL151C';'YPL152W';'YPL153C';'YPL154C';'YPL155C';'YPL156C';'YPL157W';'YPL158C';'YPL159C';'YPL160W';'YPL161C';'YPL163C';'YPL164C';'YPL165C';'YPL166W';'YPL167C';'YPL169C';'YPL170W';'YPL171C';'YPL172C';'YPL173W';'YPL174C';'YPL175W';'YPL176C';'YPL177C';'YPL178W';'YPL179W';'YPL180W';'YPL181W';'YPL183C';'YPL183W-A';'YPL184C';'YPL186C';'YPL187W';'YPL188W';'YPL189C-A';'YPL189W';'YPL190C';'YPL192C';'YPL193W';'YPL194W';'YPL195W';'YPL196W';'YPL198W';'YPL200W';'YPL201C';'YPL202C';'YPL203W';'YPL204W';'YPL206C';'YPL207W';'YPL208W';'YPL209C';'YPL210C';'YPL211W';'YPL212C';'YPL213W';'YPL214C';'YPL215W';'YPL217C';'YPL218W';'YPL219W';'YPL220W';'YPL221W';'YPL223C';'YPL224C';'YPL225W';'YPL226W';'YPL227C';'YPL228W';'YPL230W';'YPL231W';'YPL232W';'YPL233W';'YPL234C';'YPL235W';'YPL236C';'YPL237W';'YPL239W';'YPL240C';'YPL241C';'YPL242C';'YPL243W';'YPL244C';'YPL246C';'YPL248C';'YPL249C';'YPL249C-A';'YPL250C';'YPL252C';'YPL253C';'YPL254W';'YPL255W';'YPL256C';'YPL258C';'YPL259C';'YPL262W';'YPL263C';'YPL265W';'YPL266W';'YPL267W';'YPL268W';'YPL269W';'YPL270W';'YPL271W';'YPL273W';'YPL274W';'YPL281C';'YPL282C';'YPL283C';'YPR001W';'YPR002W';'YPR004C';'YPR005C';'YPR006C';'YPR007C';'YPR008W';'YPR009W';'YPR010C';'YPR016C';'YPR017C';'YPR018W';'YPR019W';'YPR020W';'YPR021C';'YPR023C';'YPR024W';'YPR025C';'YPR026W';'YPR028W';'YPR029C';'YPR030W';'YPR031W';'YPR032W';'YPR033C';'YPR034W';'YPR035W';'YPR036W';'YPR036W-A';'YPR037C';'YPR040W';'YPR041W';'YPR042C';'YPR043W';'YPR045C';'YPR046W';'YPR047W';'YPR048W';'YPR049C';'YPR051W';'YPR052C';'YPR054W';'YPR055W';'YPR056W';'YPR057W';'YPR058W';'YPR060C';'YPR061C';'YPR062W';'YPR065W';'YPR066W';'YPR067W';'YPR068C';'YPR069C';'YPR070W';'YPR072W';'YPR073C';'YPR074C';'YPR075C';'YPR079W';'YPR080W';'YPR081C';'YPR082C';'YPR083W';'YPR085C';'YPR086W';'YPR088C';'YPR091C';'YPR093C';'YPR094W';'YPR095C';'YPR096C';'YPR097W';'YPR098C';'YPR100W';'YPR101W';'YPR102C';'YPR103W';'YPR104C';'YPR105C';'YPR106W';'YPR107C';'YPR108W';'YPR110C';'YPR111W';'YPR112C';'YPR113W';'YPR115W';'YPR116W';'YPR118W';'YPR119W';'YPR120C';'YPR121W';'YPR122W';'YPR124W';'YPR125W';'YPR127W';'YPR128C';'YPR129W';'YPR131C';'YPR132W';'YPR133C';'YPR133W-A';'YPR134W';'YPR135W';'YPR137W';'YPR138C';'YPR139C';'YPR140W';'YPR141C';'YPR143W';'YPR144C';'YPR145W';'YPR148C';'YPR149W';'YPR151C';'YPR152C';'YPR153W';'YPR154W';'YPR155C';'YPR156C';'YPR158W';'YPR159W';'YPR160W';'YPR161C';'YPR162C';'YPR163C';'YPR164W';'YPR165W';'YPR166C';'YPR167C';'YPR168W';'YPR169W';'YPR171W';'YPR173C';'YPR175W';'YPR176C';'YPR178W';'YPR179C';'YPR180W';'YPR181C';'YPR182W';'YPR183W';'YPR184W';'YPR185W';'YPR186C';'YPR187W';'YPR188C';'YPR189W';'YPR190C';'YPR191W';'YPR192W';'YPR193C';'YPR194C';'YPR198W';'YPR199C';'YPR200C';'YPR201W';'YPR204W';};
    genes = upper(genes);
end

function genes = kuepfer_all_genes
    % list of 4869 genes screened by Kuepfer et al in Doi:
    % 10.1101/gr.3992505

    genes = {'YAL068C';'YAL067C';'YAL066W';'YAL065C';'YAL062W';'YAL061W';'YAL060W';'YAL059W';'YAL058W';'YAL056W';'YAL055W';'YAL053W';'YAL051W';'YAL049C';'YAL048C';'YAL046C';'YAL045C';'YAL044C';'YAL042W';'YAL043C-a';'YAL040C';'YAL039C';'YAL037W';'YAL036C';'YAL035W';'YAL034C';'YAL031C';'YAL030W';'YAL029C';'YAL028W';'YAL027W';'YAL026C';'YAL023C';'YAL022C';'YAL021C';'YAL020C';'YAL019W';'YAL018C';'YAL017W';'YAL015C';'YAL014C';'YAL013W';'YAL011W';'YAL010C';'YAL009W';'YAL008W';'YAL007C';'YAL004W';'YAL005C';'YAL002W';'YAR002W';'YAR003W';'YAR014C';'YAR015W';'YAR018C';'YAR020C';'YAR023C';'YAR027W';'YAR028W';'YAR029W';'YAR031W';'YAR030C';'YAR035W';'YAR037W';'YAR040C';'YAR042W';'YAR043C';'YAR044W';'YAR047C';'YJL218W';'YJL217W';'YJL216C';'YJL215C';'YJL214W';'YJL212C';'YJL210W';'YJL211C';'YJL209W';'YJL208C';'YJL207C';'YJL206C';'YJL206C-A';'YJL204C';'YJL201W';'YJL200C';'YJL199C';'YJL198W';'YJL197W';'YJL196C';'YJL193W';'YJL192C';'YJL191W';'YJL190C';'YJL189W';'YJL188C';'YJL187C';'YJL186W';'YJL185C';'YJL183W';'YJL181W';'YJL182C';'YJL180C';'YJL179W';'YJL178C';'YJL176C';'YJL172W';'YJL171C';'YJL170C';'YJL169W';'YJL168C';'YJL166W';'YJL165C';'YJL164C';'YJL163C';'YJL162C';'YJL161W';'YJL159W';'YJL158C';'YJL157C';'YJL155C';'YJL154C';'YJL153C';'YJL152W';'YJL151C';'YJL150W';'YJL149W';'YJL148W';'YJL147C';'YJL146W';'YJL145W';'YJL144W';'YJL142C';'YJL140W';'YJL139C';'YJL138C';'YJL135W';'YJL134W';'YJL133W';'YJL131C';'YJL130C';'YJL129C';'YJL127C';'YJL126W';'YJL124C';'YJL123C';'YJL122W';'YJL120W';'YJL121C';'YJL118W';'YJL119C';'YJL117W';'YJL116C';'YJL115W';'YJL112W';'YJL110C';'YJL108C';'YJL107C';'YJL106W';'YJL102W';'YJL100W';'YJL099W';'YJL098W';'YJL096W';'YJL095W';'YJL093C';'YJL092W';'YJL089W';'YJL088W';'YJL084C';'YJL083W';'YJL082W';'YJL080C';'YJL079C';'YJL077C';'YJL075C';'YJL073W';'YJL071W';'YJL068C';'YJL067W';'YJL066C';'YJL064W';'YJL065C';'YJL063C';'YJL062W';'YJL060W';'YJL059W';'YJL058C';'YJL057C';'YJL056C';'YJL055W';'YJL053W';'YJL052W';'YJL051W';'YJL049W';'YJL048C';'YJL047C';'YJL046W';'YJL045W';'YJL044C';'YJL043W';'YJL038C';'YJL037W';'YJL036W';'YJL030W';'YJR073C';'YJR075W';'YJR078W';'YJR079W';'YJR082C';'YJR083C';'YJR088C';'YJR092W';'YJR102C';'YJR103W';'YJR105W';'YJR108W';'YJR110W';'YJR111C';'YJR115W';'YJR127C';'YJR128W';'YJR129C';'YJR130C';'YJR135C';'YJR137C';'YJR146W';'YJR147W';'YJR149W';'YJR152W';'YJR154W';'YLL001W';'YLL002W';'YLL005C';'YLL006W';'YLL009C';'YLL010C';'YLL012W';'YLL013C';'YLL014W';'YLL015W';'YLL016W';'YLL017W';'YLL019C';'YLL020C';'YLL021W';'YLL023C';'YLL024C';'YLL025W';'YLL026W';'YLL027W';'YLL028W';'YLL029W';'YLL032C';'YLL033W';'YLL038C';'YLL039C';'YLL040C';'YLL041C';'YLL042C';'YLL043W';'YLL045C';'YLL046C';'YLL047W';'YLL051C';'YLL052C';'YLL053C';'YLL054C';'YLL055W';'YLL056C';'YLL057C';'YLL058W';'YLL060C';'YLL061W';'YLL062C';'YLL063C';'YLR001C';'YLR003C';'YLR004C';'YLR006C';'YLR011W';'YLR012C';'YLR013W';'YLR014C';'YLR015W';'YLR016C';'YLR017W';'YLR018C';'YLR019W';'YLR020C';'YLR021W';'YLR023C';'YLR024C';'YLR025W';'YLR027C';'YLR028C';'YML037C';'YML035C';'YML034W';'YML035C-A';'YML033W';'YML032C';'YML030W';'YML029W';'YML028W';'YML026C';'YML024W';'YML020W';'YML019W';'YML018C';'YML017W';'YML016C';'YML014W';'YML013W';'YML013C-A';'YML012W';'YML011C';'YML010W-A';'YML009c';'YML008C';'YML007W';'YML006C';'YML005W';'YML004C';'YML003W';'YML002W';'YML001W';'YMR002W';'YMR003W';'YMR006C';'YMR007W';'YMR008C';'YMR009W';'YMR010W';'YMR011W';'YMR012W';'YMR014W';'YMR015C';'YMR016C';'YMR017W';'YMR018W';'YMR019W';'YMR020W';'YMR021C';'YMR022W';'YMR023C';'YMR024W';'YMR025W';'YMR026C';'YMR027W';'YMR029C';'YMR030W';'YMR031W-A';'YMR031C';'YMR032W';'YMR034C';'YMR035W';'YMR036C';'YMR039C';'YMR040W';'YMR041C';'YMR042W';'YMR044W';'YMR140W';'YMR141C';'YMR143W';'YMR144W';'YMR145C';'YMR147W';'YMR148W';'YMR151W';'YMR150C';'YMR152W';'YMR153W';'YMR153C-A';'YMR155W';'YMR156C';'YMR157C';'YMR158W-A';'YMR159C';'YMR161W';'YMR162C';'YMR163C';'YMR164C';'YMR166C';'YMR167W';'YMR169c';'YMR170C';'YMR172C-A';'YMR173W-A';'YMR174c';'YMR175w';'YMR176W';'YMR177W';'YMR178W';'YMR179W';'YMR180C';'YMR182C';'YMR183C';'YMR184W';'YMR185W';'YMR186W';'YMR187C';'YMR188C';'YMR189W';'YMR190C';'YMR191W';'YMR192W';'YMR193W';'YMR194W';'YMR193C-A';'YMR195W';'YMR196W';'YMR198W';'YMR199W';'YMR201C';'YMR202W';'YMR204C';'YMR205C';'YMR206W';'YMR207C';'YMR210W';'YMR214W';'YMR215W';'YMR216C';'YMR219W';'YMR221C';'YMR222C';'YMR223W';'YMR224C';'YMR225C';'YMR226C';'YMR228W';'YMR230W';'YMR231W';'YMR232W';'YMR233W';'YMR234W';'YMR237W';'YMR238W';'YMR241W';'YMR242C';'YMR243C';'YMR244W';'YMR245W';'YMR244C-A';'YMR246W';'YMR247C';'YMR250W';'YMR251W';'YMR251W-A';'YMR252C';'YMR253C';'YMR254C';'YMR255W';'YMR256C';'YMR257C';'YMR258C';'YMR259C';'YMR261C';'YMR262W';'YMR263W';'YMR264W';'YMR265C';'YMR266W';'YMR267W';'YMR269W';'YMR272C';'YMR273C';'YMR274C';'YMR275C';'YMR276W';'YMR278W';'YMR280C';'YMR282C';'YMR283C';'YMR284W';'YMR285C';'YMR286W';'YMR287C';'YMR289W';'YMR291W';'YMR292W';'YMR293C';'YMR294W';'YMR294W-A';'YMR295C';'YMR297W';'YMR299C';'YMR300C';'YMR302C';'YMR303C';'YMR304W';'YMR304C-A';'YMR305C';'YMR306C-A';'YMR307W';'YMR310C';'YNL339C';'YNL338W';'YNL336W';'YNL335W';'YNL334C';'YNL333W';'YNL332W';'YNL330C';'YNL329C';'YNL328C';'YNL327W';'YNL326C';'YNL324W';'YNL325C';'YNL323W';'YNL322C';'YNL321W';'YNL320W';'YNL319W';'YNL318C';'YNL314W';'YNL311C';'YNL309W';'YNL307C';'YNL305C';'YNL304W';'YNL303W';'YNL302C';'YNL301C';'YNL299W';'YNL298W';'YNL296W';'YNL297C';'YNL295W';'YNL294C';'YNL293W';'YNL292W';'YNL291C';'YNL289W';'YNL288W';'YNL286W';'YNL285W';'YNL283C';'YNL281W';'YNL280C';'YNL278W';'YNL277W';'YNL276C';'YNL275W';'YNL273W';'YNL271C';'YNL270C';'YNL269W';'YNL268W';'YNL266W';'YNL265C';'YNL264C';'YNL259C';'YNL257C';'YNL255C';'YNL254C';'YNL253W';'YNL249C';'YNL248C';'YNL246W';'YNL242W';'YNL241C';'YNL239W';'YNL238W';'YNL237W';'YNL236W';'YNL235C';'YNL234W';'YNL233W';'YNL231C';'YNL230C';'YNL229C';'YNL228W';'YNL226W';'YNL227C';'YNL225C';'YNL224C';'YNL223W';'YNL219C';'YNL218W';'YNL217W';'YNL215W';'YNL214W';'YNL213C';'YNL212W';'YNL211C';'YNL208W';'YNL206C';'YNL205C';'YNL204C';'YNL202W';'YNL203C';'YNL201C';'YNL200C';'YNL199C';'YNL198C';'YNL197C';'YNL196C';'YNL195C';'YNL194C';'YNL193W';'YNL192W';'YNL191W';'YNL190W';'YNL187W';'YNL184C';'YNL183C';'YNL179C';'YNL177C';'YNL176C';'YNL175C';'YNL173C';'YNL170W';'YNL171C';'YNL169C';'YNL168C';'YNL167C';'YNL166C';'YNL165W';'YNL164C';'YNL162W';'YNL160W';'YNL159C';'YNL157W';'YNL156C';'YNL155W';'YNL154C';'YNL153C';'YNL148C';'YOR001W';'YOR002W';'YOR003W';'YOR005C';'YOR006C';'YOR007C';'YOR008C';'YOR009W';'YOR010C';'YOR011W';'YOR012W';'YOR013W';'YOR014W';'YOR015W';'YOR016C';'YOR017W';'YOR018W';'YOR019W';'YOR021C';'YOR022C';'YOR023C';'YOR024W';'YOR025W';'YOR026W';'YOR027W';'YOR028C';'YOR029W';'YOR030W';'YOR031W';'YOR032C';'YOR033C';'YOR034C';'YOR035C';'YOR036W';'YOR037W';'YOR038C';'YOR039W';'YOR040W';'YOR041C';'YOR042W';'YOR043W';'YOR044W';'YOR045W';'YOR047C';'YOR049C';'YOR050C';'YOR051C';'YOR052C';'YOR053W';'YOR054c';'YOR055W';'YOR058C';'YOR059C';'YOR061W';'YOR062C';'YOR064C';'YOR065W';'YOR066W';'YOR067C';'YOR068C';'YOR069W';'YOR070C';'YOR071C';'YOR072W';'YOR073W';'YOR076C';'YOR078W';'YOR079C';'YOR080W';'YOR081C';'YOR082C';'YOR083W';'YOR084W';'YOR085W';'YOR086C';'YOR087W';'YOR088W';'YOR089C';'YOR090C';'YOR091W';'YOR092W';'YOR093C';'YOR094W';'YOR097C';'YOR099W';'YOR100C';'YOR101W';'YOR104W';'YOR105W';'YOR106W';'YOR107W';'YOR108W';'YOR109W';'YOR111W';'YOR112W';'YOR113W';'YOR114W';'YOR115C';'YOR118W';'YOR120W';'YOR121C';'YOR123C';'YOR124C';'YOR125C';'YOR126C';'YOR127W';'YOR128C';'YOR129C';'YOR130C';'YOR131C';'YOR132W';'YOR133W';'YOR134W';'YOR135C';'YOR136W';'YOR137C';'YOR138C';'YOR139C';'YOR140W';'YOR141C';'YOR142W';'YOR144C';'YOR150W';'YOR152C';'YOR153W';'YOR154W';'YOR155C';'YOR156C';'YOR158W';'YOR161C';'YOR162C';'YOR163W';'YOR164C';'YOR165W';'YOR166C';'YOR167C';'YOR170W';'YOR171C';'YOR172W';'YOR173W';'YOR175C';'YOR177C';'YOR178C';'YOR182C';'YOR183W';'YOR184W';'YOR185C';'YOR186W';'YOR187W';'YOR188W';'YOR189W';'YOR190W';'YOR191W';'YOR192C';'YOR193W';'YOR195W';'YOR196C';'YOR197W';'YOR198C';'YOR199W';'YOR200W';'YOR201C';'YOR202W';'YOR205C';'YOR208W';'YOR209C';'YOR211C';'YOR212W';'YOR213C';'YOR214C';'YOR215C';'YOR216C';'YOR219C';'YOR220W';'YOR221C';'YOR222W';'YOR223W';'YOR225W';'YOR226C';'YOR227W';'YOR228C';'YOR229W';'YOR230W';'YOR231W';'YOR233W';'YOR234C';'YOR235W';'YOR237W';'YOR238W';'YOR239W';'YOR240W';'YOR241W';'YOR242C';'YOR243C';'YOR245C';'YOR246C';'YOR247W';'YOR248W';'YOR251C';'YOR252W';'YOR253W';'YOR255W';'YOR258W';'YOR263C';'YOR264W';'YOR277C';'YOR279C';'YOR280C';'YOR283W';'YOR284W';'YOR285W';'YOR286W';'YOR288C';'YOR289W';'YOR290C';'YOR291W';'YOR292C';'YOR293W';'YOR295W';'YOR296W';'YOR297C';'YOR298W';'YOR299W';'YOR300W';'YOR301W';'YOR302W';'YOR303W';'YOR304C-A';'YOR304W';'YOR305W';'YOR306C';'YOR307C';'YOR308C';'YOR309C';'YOR311C';'YOR312C';'YOR313C';'YOR314W';'YOR315W';'YOR316C';'YOR318C';'YOR320C';'YOR321W';'YOR322C';'YOR323C';'YOR324C';'YOR325W';'YOR327C';'YOR328W';'YOR330C';'YOR331C';'YOR332W';'YOR333C';'YOR334W';'YOR337W';'YOR338W';'YOR339C';'YOR342C';'YOR343C';'YOR344C';'YOR345C';'YOR346W';'YOR347C';'YOR348C';'YOR349W';'YOR350C';'YOR351C';'YOR352W';'YOR354C';'YOR355W';'YOR356W';'YOR357C';'YOR358W';'YOR359W';'YOR360C';'YOR363C';'YOR365C';'YOR366W';'YOR367W';'YOR368W';'YOR369C';'YOR371C';'YOR374W';'YOR375C';'YOR376W';'YOR377W';'YOR378W';'YOR379C';'YOR380W';'YOR381W';'YOR382W';'YOR383C';'YOR384W';'YOR385W';'YOR386W';'YOL001W';'YOL002C';'YOL003C';'YOL004W';'YOL006C';'YOL007C';'YOL008W';'YOL009C';'YOL011W';'YOL012C';'YOL013C';'YOL014W';'YOL015W';'YOL017W';'YOL018C';'YOL019W';'YOL020W';'YOL023W';'YOL024W';'YOL025W';'YOL027C';'YOL028C';'YOL029C';'YOL030W';'YOL031C';'YOL032W';'YOL033W';'YOL035C';'YOL036W';'YOL037C';'YOL039W';'YOL041C';'YOL042W';'YOL043C';'YOL044W';'YOL045W';'YOL046C';'YOL047C';'YOL048C';'YOL049W';'YOL050C';'YOL051W';'YOL052C';'YOL053C-A';'YOL053W';'YOL054W';'YOL055C';'YOL056W';'YOL057W';'YOL058W';'YOL059W';'YOL060C';'YOL061W';'YOL062C';'YOL063C';'YOL064C';'YOL065C';'YOL067C';'YOL068C';'YOL070C';'YOL071W';'YOL072W';'YOL075C';'YOL076W';'YOL079W';'YOL080C';'YOL081W';'YOL082W';'YOL083W';'YOL084W';'YOL085C';'YPL274W';'YPL273W';'YPL272C';'YPL271W';'YPL270W';'YPL269W';'YPL267W';'YPL265W';'YPL264C';'YPL263C';'YPL262W';'YPL260W';'YPL261C';'YPL259C';'YPL258C';'YPL257W';'YPL256C';'YPL254W';'YPL253C';'YPL250C';'YPL249C';'YPL248C';'YPL247C';'YPL246C';'YPL245W';'YPL244C';'YPL241C';'YPL240C';'YPL239W';'YPL236C';'YPL234C';'YPL232W';'YPL230W';'YPL229W';'YPL227C';'YPL226W';'YPL225W';'YPL223C';'YPL222W';'YPL221W';'YPL220W';'YPL219W';'YPL216W';'YPL215W';'YPL214C';'YPL213W';'YPL212C';'YPL208W';'YPL207W';'YPL206C';'YPL205C';'YPL203W';'YPL202C';'YPL201C';'YPL200W';'YPL199C';'YPL198W';'YPL197C';'YPL196W';'YPL195W';'YPL194W';'YPL193W';'YPL192C';'YPL191C';'YPL189W';'YPL188W';'YPL187W';'YPL185W';'YPL186C';'YPL184C';'YPL181W';'YPL182C';'YPL180W';'YPL179W';'YPL178W';'YPL177C';'YPL176C';'YPL174C';'YPL173W';'YPL172C';'YPL171C';'YPL170W';'YPL168W';'YPL167C';'YPL166W';'YPL165C';'YPL164C';'YPL163C';'YPL162C';'YPL161C';'YPL159C';'YPL157W';'YPL156C';'YPL155C';'YPL154C';'YPL152W';'YPL150W';'YPL149W';'YPL147W';'YPL145C';'YPL144W';'YPL141C';'YPL140C';'YPL139C';'YPL138C';'YPL136W';'YPL135W';'YPL133C';'YPL130W';'YPL129W';'YPL127C';'YPL125W';'YPL123C';'YPL121C';'YPL120W';'YPL119C';'YPL118W';'YPL116W';'YPL115C';'YPL114W';'YPL113C';'YPL112C';'YPL111W';'YPL110C';'YPL109C';'YPL108W';'YPL107W';'YPL106C';'YPL105C';'YPL104W';'YPL103C';'YPL101W';'YPL102C';'YPL100W';'YPL099C';'YPL098C';'YPL097W';'YPL096W';'YPL095C';'YPL092W';'YEL001C';'YEL003W';'YEL004W';'YEL005C';'YEL006W';'YEL007W';'YEL008W';'YEL009C';'YEL010W';'YEL013W';'YEL014C';'YEL015W';'YEL016C';'YEL017C-A';'YEL017W';'YEL018W';'YEL020C';'YEL023C';'YEL024W';'YEL025C';'YEL027W';'YEL028W';'YEL029C';'YEL030W';'YEL031W';'YEL033W';'YEL036C';'YEL037C';'YEL038W';'YEL039C';'YEL040W';'YEL041W';'YEL042W';'YEL043W';'YEL044W';'YEL045C';'YEL046C';'YEL047C';'YEL048C';'YEL049W';'YEL050C';'YEL051W';'YEL052W';'YEL053C';'YEL054C';'YEL056W';'YEL057C';'YEL059W';'YEL060C';'YEL061C';'YEL062W';'YEL063C';'YEL064C';'YEL065W';'YEL066W';'YEL067C';'YEL068C';'YEL071W';'YEL072W';'YER001W';'YER002W';'YER004W';'YER005W';'YER007C-A';'YER007W';'YER010C';'YER011W';'YER014C-A';'YER016W';'YER017C';'YER019W';'YER019C-A';'YER020W';'YER024W';'YER030W';'YER032W';'YER033C';'YER034W';'YER035W';'YER038W-A';'YER039C';'YER040W';'YER041W';'YER042W';'YER044C-A';'YER045C';'YER046W-A';'YER047C';'YER048C';'YER049W';'YER050C';'YER051W';'YER052C';'YER053C';'YER054C';'YER056C';'YER056C-A';'YER057C';'YER058W';'YER059W';'YER060W';'YER060W-A';'YER061C';'YER062C';'YER065C';'YER066C-A';'YER067W';'YER067C-A';'YER068W';'YER068C-A';'YER069W';'YER070W';'YER071C';'YER072W';'YER073W';'YER074W';'YER075C';'YER079W';'YER080W';'YER081W';'YER083C';'YER084W';'YER085C';'YER086W';'YER087C-A';'YHL047C';'YHL046C';'YHL045W';'YHL044W';'YHL043W';'YHL042W';'YHL041W';'YHL040C';'YHL038C';'YHL037C';'YHL036W';'YHL035C';'YHL034C';'YHL033C';'YHL032C';'YHL031C';'YHL030W';'YHL029C';'YHL028W';'YHL027W';'YHL026C';'YHL023C';'YHL022C';'YHL021C';'YHL020C';'YHL019C';'YHL017W';'YHL016C';'YHL014C';'YHL013C';'YHL012W';'YHL010C';'YHL009C';'YHL008C';'YHL007C';'YHL006C';'YHL005C';'YHL003C';'YHR001W-A';'YHR010W';'YHR011W';'YHR012W';'YHR013C';'YHR014W';'YHR015W';'YHR018C';'YHR021C';'YHR022C';'YHR028C';'YHR029C';'YHR030C';'YHR031C';'YHR033W';'YHR034C';'YHR035W';'YHR037W';'YHR038W';'YHR039C';'YHR044C';'YHR046C';'YHR047C';'YHR048W';'YHR049W';'YHR049C-A';'YHR050W';'YHR051W';'YHR057C';'YHR060W';'YHR061C';'YHR066W';'YHR073W';'YHR075C';'YHR076W';'YHR077C';'YHR078W';'YHR079C';'YHR080C';'YHR081W';'YHR082C';'YHR086W';'YHR087W';'YHR091C';'YHR092C';'YHR093W';'YHR094C';'YHR095W';'YHR096C';'YHR097C';'YHR100C';'YHR103W';'YHR104W';'YHR105W';'YHR106W';'YHR108W';'YHR109W';'YHR110W';'YHR111W';'YHR112C';'YHR113W';'YHR114W';'YHR115C';'YHR116W';'YHR117W';'YHR120W';'YHR121W';'YHR123W';'YHR124W';'YHR125W';'YHR126C';'YHR129C';'YHR130C';'YHR132C';'YHR133C';'YHR134W';'YHR135C';'YHR136C';'YHR137W';'YHR138C';'YHR139C';'YHR139C-A';'YIL006W';'YIL007C';'YIL008W';'YIL009W';'YIL010W';'YIL015C-A';'YIL018W';'YIL038C';'YIL042C';'YIL047C';'YIL052C';'YIL054W';'YIL055C';'YIL056W';'YIL059C';'YIL060W';'YIL066C';'YIL067C';'YIL069C';'YIL070C';'YIL071C';'YIL074C';'YIL085C';'YIL089W';'YIL092W';'YIL094C';'YIL001W';'YIL002C';'YIL005W';'YIL011W';'YIL012W';'YIL013C';'YIL014W';'YIL015W';'YIL016W';'YIL017C';'YIL020C';'YIL023C';'YIL024C';'YIL025C';'YIL027C';'YIL028W';'YIL029C';'YIL032C';'YIL034C';'YIL035C';'YIL036W';'YIL037C';'YIL039W';'YIL040W';'YIL041W';'YIL043C';'YIL044C';'YIL045W';'YIL049W';'YIL050W';'YIL053W';'YIL057C';'YIL064W';'YIL065C';'YIL072W';'YIL073C';'YIL076W';'YIL077C';'YIL079C';'YIL084C';'YIL086C';'YIL087C';'YIL088C';'YIL090W';'YIL093C';'YIL095W';'YIL096C';'YIL097W';'YIL098C';'YIL099W';'YIL100W';'YIL101C';'YIL102C';'YIL103W';'YIL105C';'YIL107C';'YIL108W';'YIL110W';'YIL111W';'YIL112W';'YIL113W';'YIL114C';'YIL116W';'YIL117C';'YIL119C';'YIL120W';'YIL121W';'YIL122W';'YIL123W';'YIL124W';'YIL125W';'YIL128W';'YIL130W';'YIL131C';'YIL132C';'YIL133C';'YIL135C';'YIL136W';'YIL137C';'YIL138C';'YIL139C';'YIL140W';'YIL141W';'YIL145C';'YIL146C';'YIL148W';'YIL149C';'YIL152W';'YIL153W';'YIL154C';'YIL155C';'YIL156W';'YIL157C';'YIL158W';'YIL159W';'YIL160C';'YIL161W';'YIL162W';'YIL163C';'YIL164C';'YIL165C';'YIL166C';'YIL167W';'YIL168W';'YIL170W';'YIL173W';'YIR001C';'YIR002C';'YIR003W';'YIR004W';'YIR005W';'YIR007W';'YIR009W';'YIR013C';'YIR014W';'YIR016W';'YAL064C-A';'YBL091C-A';'YBR269C';'YBR271W';'YBR273C';'YBR274W';'YBR277C';'YBR278W';'YBR279W';'YBR281C';'YBR282W';'YBR283C';'YBR284W';'YBR285W';'YBR286W';'YBR290W';'YBR291C';'YBR292C';'YBR293W';'YBR295W';'YBR296C';'YBR297W';'YBR298C';'YBR300C';'YCL001W-A';'YCR020W-B';'YCR024C';'YCR024C-A';'YCR025C';'YCR026C';'YCR027C';'YCR028C';'YCR031C';'YCR034W';'YCR036W';'YCR037C';'YCR043C';'YCR044C';'YCR045C';'YCR049C';'YCR050C';'YCR051W';'YCR059C';'YCR061W';'YCR063W';'YCR065W';'YCR066W';'YCR068W';'YCR071C';'YCR073W-A';'YCR076C';'YCR077C';'YCR079W';'YCR081W';'YCR082W';'YCR085W';'YCR086W';'YCR087C-A';'YCR087W';'YEL012W';'YER031C';'YER046W';'YER063W';'YER066W';'YGR188C';'YGR201C';'YGR204W';'YHL002W';'YHL011C';'YHL039W';'YHR003C';'YHR004C';'YHR005C';'YHR006W';'YHR008C';'YHR009C';'YHR025W';'YHR026W';'YHR041C';'YHR059W';'YHR067W';'YHR127W';'YHR131C';'YHR180W';'YHR185C';'YHR194W';'YLL007C';'YMR154C';'YNL274C';'YOL141W';'YOL143C';'YOL150C';'YOL158C';'YOL159C';'YOL160W';'YOL162W';'YOL163W';'YOR008C-A';'YBR191W';'YCL035C';'YDR071C';'YDR074W';'YER027C';'YER037W';'YGR155W';'YLR192C';'YLR237W';'YLR246W';'YLR334C';'YLR346C';'YLR358C';'YLR361C';'YLR370C';'YLR382C';'YLR394W';'YLR406C';'YML022W';'YML027W';'YML036W';'YML038C';'YML041C';'YML042W';'YML047C';'YML075C';'YML076C';'YML086C';'YMR048W';'YMR135W-A';'YMR137C';'YMR138W';'YMR139W';'YMR160W';'YMR173W';'YMR198W';'YOR298C-A';'YOR364W';'YPL148C';'YPL183C';'YPL183W-A';'YPL189W';'YPL224C';'YAL012W';'YAL016W';'YAL047C';'YAL054C';'YAL058C-A';'YAR050W';'YCL006C';'YCL022C';'YCL023C';'YCL038C';'YCL058C';'YCL074W';'YCL075W';'YCL076W';'YGL199C';'YGL214W';'YGL217C';'YGL218W';'YGL235W';'YGR011W';'YGR018C';'YGR022C';'YGR025W';'YJR069C';'YJR070C';'YJR074W';'YJR077C';'YJR080C';'YJR084W';'YJR087W';'YJR090C';'YJR091C';'YJR094C';'YJR094W-A';'YJR095W';'YJR096W';'YJR097W';'YJR098C';'YJR099W';'YJR100C';'YJR104C';'YJR106W';'YJR107W';'YJR109C';'YJR113C';'YJR116W';'YJR117W';'YJR118C';'YJR119C';'YJR120W';'YJR121W';'YJR122W';'YJR124C';'YJR125C';'YJR126C';'YJR131W';'YJR133W';'YJR134C';'YJR139C';'YJR140C';'YJR142W';'YJR144W';'YJR145C';'YJR148W';'YJR150C';'YJR153W';'YKL005C';'YKL030W';'YLR146C';'YLR286C';'YLR343W';'YML067C';'YML068W';'YML072C';'YMR136W';'YMR172W';'YOR300W';'YOR309C';'YJL003W';'YJL004C';'YJL006C';'YJL007C';'YJL012C';'YJL013C';'YJL016W';'YJL017W';'YJL020C';'YJL021C';'YJL022W';'YJL023C';'YJL024C';'YJL027C';'YJL028W';'YJL029C';'YJR001W';'YJR004C';'YJR005W';'YJR008W';'YJR009C';'YJR010C-A';'YJR010W';'YJR011C';'YJR014W';'YJR015W';'YJR018W';'YJR019C';'YJR020W';'YJR021C';'YJR024C';'YJR025C';'YJR026W';'YJR030C';'YJR031C';'YJR032W';'YJR033C';'YJR034W';'YJR035W';'YJR036C';'YJR037W';'YJR038C';'YJR039W';'YJR040W';'YJR043C';'YJR044C';'YJR047C';'YJR048W';'YJR049C';'YJR050W';'YJR051W';'YJR052W';'YJR053W';'YJR054W';'YJR055W';'YJR056C';'YJR058C';'YJR059W';'YJR060W';'YJR061W';'YJR062C';'YJR063W';'YJR066W';'YBR189W';'YCR095C';'YCR102W-A';'YDL133C-A';'YDR058C';'YDR174W';'YDR202C';'YDR205W';'YDR445C';'YDR537C';'YFR039C';'YGL219C';'YGR028W';'YGR032W';'YGR038W';'YGR040W';'YGR050C';'YGR053C';'YGR063C';'YGR086C';'YGR089W';'YGR092W';'YGR093W';'YGR106C';'YGR110W';'YGR117C';'YGR238C';'YGR239C';'YGR248W';'YGR250C';'YJL129C';'YJL132W';'YJL136C';'YJL137C';'YJL139C';'YJL140W';'YJL141C';'YJL151C';'YJL160C';'YJL161W';'YJL163C';'YJL165C';'YJL172W';'YJL175W';'YJL177W';'YJL184W';'YJL189W';'YJL191W';'YJL196C';'YJL200C';'YJL206C';'YJL213W';'YKL096W-A';'YKL115C';'YKL139W';'YKL194C';'YKL201C';'YKL202W';'YKL204W';'YKL215C';'YKL220C';'YKR010C';'YKR019C';'YKR023W';'YKR027W';'YKR028W';'YKR029C';'YKR034W';'YKR036C';'YKR039W';'YKR040C';'YKR041W';'YKR046C';'YKR053C';'YML035C';'YFL013W-A';'YFL014W';'YFL019C';'YFL042C';'YFR019W';'YFR024C';'YFR025C';'YFR030W';'YHR146W';'YHR171W';'YJL042W';'YJL070C';'YJL078C';'YJL094C';'YJL101C';'YJL105W';'YJL128C';'YKR085C';'YKR094C';'YKR095W';'YKR096W';'YKR102W';'YLR110C';'YLR390W-A';'YLR439W';'YLR442C';'YLR455W';'YML066C';'YML115C';'YMR037C';'YMR074C';'YMR104C';'YMR118C';'YMR119W';'YNR052C';'YNR055C';'YNR069C';'YOL125W';'YOL147C';'YOL153C';'YPL268W';'YPR007C';'YPR008W';'YPR013C';'YPR022C';'YPR023C';'YPR024W';'YPR026W';'YPR031W';'YPR037C';'YPR043W';'YPR050C';'YPR064W';'YPR067W';'YPR078C';'YAR002C-A';'YBR083W';'YBR084C-A';'YBR090C';'YBR100W';'YBR112C';'YBR125C';'YBR131W';'YBR150C';'YBR168W';'YBR169C';'YBR270C';'YBR272C';'YBR275C';'YBR276C';'YBR280C';'YBR287W';'YBR288C';'YBR289W';'YBR294W';'YBR301W';'YCL026C-A';'YCR028C-A';'YCR030C';'YCR032W';'YCR033W';'YCR046C';'YCR047C';'YCR048W';'YCR053W';'YCR060W';'YCR062W';'YCR067C';'YCR069W';'YCR073C';'YCR075C';'YCR083W';'YCR084C';'YCR088W';'YCR089W';'YDL194W';'YDR007W';'YDR048C';'YFR011C';'YFR013W';'YNL051W';'YNL052W';'YNL056W';'YNL065W';'YNL066W';'YNL067W';'YNL068C';'YNL070W';'YNL071W';'YNL072W';'YNL073W';'YNL074C';'YNL076W';'YNL078W';'YNL079C';'YNL080C';'YNL082W';'YNL083W';'YNL085W';'YNL087W';'YNL090W';'YNL091W';'YNL093W';'YNL095C';'YNL097C';'YNL099C';'YNL100W';'YNL104C';'YNL105W';'YNL106C';'YNL107W';'YNL115C';'YNL119W';'YNL120C';'YNL121C';'YNL125C';'YNL130C';'YAL024C';'YBR299W';'YCR107W';'YDR242W';'YDR326C';'YDR417C';'YDR444W';'YDR461W';'YDR493W';'YDR500C';'YDR502C';'YDR506C';'YDR512C';'YDR515W';'YER089C';'YFL001W';'YFL003C';'YFL004W';'YFL007W';'YFL010C';'YFL010W-A';'YFL012W';'YFL013C';'YFL016C';'YFL033C';'YGR122C-A';'YGR162W';'YGR252W';'YGR254W';'YGR255C';'YGR257C';'YGR258C';'YGR271W';'YGR272C';'YGR273C';'YGR276C';'YGR289C';'YGR291C';'YGR292W';'YGR295C';'YHR132W-A';'YIL030C';'YIL058W';'YIL092W';'YIR023W';'YIR030C';'YIR032C';'YIR043C';'YIR044C';'YJR003C';'YJR055W';'YKL053C-A';'YKR106W';'YMR191W';'YMR322C';'YNL138W';'YNL140C';'YNL142W';'YNL315C';'YOL151W';'YOL152W';'YOL155C';'YOR265W';'YOR266W';'YOR267C';'YOR268C';'YOR269W';'YOR270C';'YOR271C';'YOR273C';'YOR274W';'YOR275C';'YOR276W';'YOR298C-A';'YOR302W';'YOR303W';'YPL004C';'YPL017C';'YPL027W';'YPL034W';'YPL036W';'YPL078C';'YPL137C';'YBR020W';'YBR075W';'YDR417C';'YFL033C';'YFL063W';'YJL103C';'YML073C';'YNL011C';'YNL014W';'YNL047C';'YNL053W';'YNL055C';'YNL059C';'YNL069C';'YNL086W';'YNL089C';'YNL096C';'YNL109W';'YNL111C';'YNL147W';'YNL220W';'YNL268W';'YNL284C';'YNL315C';'YNR033W';'YOL148C';'YOL151W';'YOL152W';'YPL158C';'YPL183W-A';'YPL194W';'YPR011C';'YPR021C';'YPR083W';'YPR091C';'YPR118W';'YPR133W-A';'YPR151C';'YCR090C';'YCR091W';'YCR092C';'YCR094W';'YCR098C';'YCR099C';'YCR100C';'YCR101C';'YCR102C';'YCR105W';'YCR106W';'YDL130W-A';'YDR363W-A';'YDR525W-A';'YDR535C';'YDR536W';'YDR538W';'YDR539W';'YDR540C';'YDR541C';'YER039C-A';'YER091C-A';'YER144C';'YER188W';'YFL034C-A';'YFR032C';'YFR032C-A';'YFR033C';'YFR034C';'YFR035C';'YFR036W';'YFR038W';'YFR040W';'YFR041C';'YFR043C';'YFR044C';'YFR045W';'YFR046C';'YFR047C';'YFR048W';'YFR049W';'YFR053C';'YFR054C';'YFR055W';'YFR056C';'YFR057W';'YGR219W';'YGR220C';'YGR221C';'YGR222W';'YGR223C';'YGR224W';'YGR225W';'YGR226C';'YGR227W';'YGR228W';'YGR229C';'YGR230W';'YGR231C';'YGR232W';'YGR233C';'YGR234W';'YGR235C';'YGR236C';'YGR237C';'YGR240C';'YGR241C';'YGR242W';'YGR243W';'YGR244C';'YGR247W';'YGR249W';'YGR256W';'YGR259C';'YGR260W';'YGR261C';'YGR262C';'YGR263C';'YGR266W';'YGR268C';'YGR269W';'YGR270W';'YGR275W';'YGR279C';'YGR281W';'YGR282C';'YGR283C';'YGR284C';'YGR285C';'YGR286C';'YGR287C';'YGR288W';'YGR290W';'YHR021W-A';'YHR039C-B';'YHR079C-B';'YIL009C-A';'YIR017C';'YIR018W';'YIR019C';'YIR020C';'YIR020W-B';'YIR021W';'YIR024C';'YIR025W';'YIR026C';'YIR027C';'YIR028W';'YIR029W';'YIR031C';'YIR033W';'YIR034C';'YIR035C';'YIR036C';'YIR037W';'YIR038C';'YIR039C';'YIR042C';'YKL033W-A';'YKL162C-A';'YKR035W-A';'YKR066C';'YKR067W';'YKR069W';'YKR070W';'YKR072C';'YKR073C';'YKR074W';'YKR075C';'YKR076W';'YKR077W';'YKR078W';'YKR080W';'YKR082W';'YKR084C';'YKR087C';'YKR088C';'YKR089C';'YKR090W';'YKR091W';'YKR092C';'YKR093W';'YKR097W';'YKR098C';'YKR099W';'YKR100C';'YKR101W';'YKR103W';'YKR104W';'YKR105C';'YLL018C-A';'YLR262C-A';'YLR422W';'YLR423C';'YLR425W';'YLR426W';'YLR427W';'YLR428C';'YLR429W';'YLR431C';'YLR432W';'YLR433C';'YLR434C';'YLR435W';'YLR436C';'YLR437C';'YLR438W';'YLR441C';'YLR443W';'YLR444C';'YLR445W';'YLR446W';'YLR447C';'YLR448W';'YLR449W';'YLR450W';'YLR452C';'YLR453C';'YLR454W';'YLR456W';'YLR460C';'YLR461W';'YML009C';'YML010C-B';'YML021C';'YML081C-A';'YMR060C';'YMR158C-B';'YMR169C';'YMR174C';'YMR175W';'YMR194C-A';'YMR326C';'YNR032C-A';'YNR050C';'YNR051C';'YNR056C';'YNR057C';'YNR058W';'YNR059W';'YNR060W';'YNR061C';'YNR062C';'YNR063W';'YNR064C';'YNR065C';'YNR066C';'YNR067C';'YNR068C';'YER101C';'YER103W';'YER106W';'YER108C';'YER109C';'YER110C';'YER111C';'YER113C';'YER114C';'YER115C';'YER116C';'YER117W';'YER118C';'YER119C';'YER119C-A';'YER120W';'YER121W';'YER122C';'YER123W';'YER124C';'YER128W';'YER129W';'YER130C';'YER131W';'YER132C';'YER134C';'YER135C';'YER137C';'YER139C';'YER140W';'YER141W';'YER142C';'YER143W';'YER145C';'YER149C';'YER150W';'YER151C';'YER152C';'YER153C';'YER154W';'YER155C';'YER156C';'YER158C';'YER161C';'YER162C';'YER163C';'YER164W';'YER166W';'YER167W';'YER169W';'YER170W';'YER173W';'YER174C';'YER175C';'YER176W';'YER177W';'YER178W';'YER179W';'YER180C';'YER181C';'YER182W';'YER183C';'YER184C';'YER185W';'YER186C';'YER187W';'YMR052C-A';'YMR052W';'YMR053C';'YMR054W';'YMR055C';'YMR056C';'YMR057C';'YMR058W';'YMR062C';'YMR063W';'YMR064W';'YMR065W';'YMR066W';'YMR067C';'YMR068W';'YMR069W';'YMR070W';'YMR071C';'YMR072W';'YMR073C';'YMR075C-A';'YMR075W';'YMR077C';'YMR078C';'YMR080C';'YMR081C';'YMR082C';'YMR083W';'YMR084W';'YMR085W';'YMR086C-A';'YMR086W';'YMR087W';'YMR088C';'YMR089C';'YMR090W';'YMR091C';'YMR092C';'YNR070W';'YNR071C';'YNR072W';'YNR073C';'YNR074C';'YNR075W';'YOL013W-A';'YOL086C';'YOL087C';'YOL088C';'YOL089C';'YOL090W';'YOL091W';'YOL092W';'YOL093W';'YOL095C';'YOL096C';'YOL098C';'YOL099C';'YOL100W';'YOL101C';'YOL103W';'YOL104C';'YOL105C';'YOL106W';'YOL107W';'YOL108C';'YOL109W';'YOL110W';'YOL111C';'YOL112W';'YOL113W';'YOL114C';'YOL115W';'YOL116W';'YOL117W';'YOL118C';'YOL119C';'YOL121C';'YOL122C';'YOL124C';'YOL126C';'YOL128C';'YOL129W';'YOL131W';'YOL132W';'YOL136C';'YOL137W';'YOL138C';'YBR232C';'YDR424C';'YEL011W';'YER064C';'YER077C';'YER078C';'YER088C';'YER090W';'YER091C';'YER092W';'YER093C-A';'YER095W';'YER096W';'YER097W';'YER098W';'YGR134W';'YGR210C';'YHL024W';'YHL025W';'YHR016C';'YHR017W';'YHR032W';'YHR045W';'YHR064C';'YHR140W';'YHR162W';'YHR168W';'YHR181W';'YHR191C';'YHR193C';'YLL030C';'YLL044W';'YLL048C';'YLL049W';'YLL059C';'YLR030W';'YLR031W';'YLR032W';'YLR034C';'YLR035C';'YLR036C';'YLR037C';'YLR038C';'YLR039C';'YLR040C';'YLR041W';'YLR050C';'YLR052W';'YMR142C';'YMR158W';'YMR171C';'YMR181C';'YMR209C';'YMR271C';'YMR279C';'YMR306W';'YMR311C';'YMR312W';'YMR313C';'YMR315W';'YMR316C-A';'YMR316C-B';'YMR316W';'YMR317W';'YMR318C';'YMR319C';'YMR320W';'YNL250W';'YNL252C';'YNL279W';'YNL300W';'YNL316C';'YOL016C';'YOR096W';'YOR306C';'YOR317W';'YPL132W';'YPL134C';'YML070W';'YML071C';'YML074C';'YML090W';'YML094W';'YML095C';'YML095C-A';'YML096W';'YML097C';'YML099C';'YML100W';'YML100W-A';'YML101C';'YML102C-A';'YML102W';'YML103C';'YML104C';'YML106W';'YML107C';'YML108W';'YML109W';'YML110C';'YML111W';'YML112W';'YML113W';'YML116W';'YML117W';'YML117W-A';'YML118W';'YML119W';'YML120C';'YML121W';'YML122C';'YML123C';'YML124C';'YML128C';'YML129C';'YML131W';'YMR004W';'YMR095C';'YMR096W';'YMR097C';'YMR098C';'YMR099C';'YMR100W';'YMR101C';'YMR102C';'YMR103C';'YMR105C';'YMR106C';'YMR107W';'YMR109W';'YMR110C';'YMR111C';'YMR114C';'YMR115W';'YMR116C';'YMR119W-A';'YMR120C';'YMR121C';'YMR122C';'YMR123W';'YMR124W';'YMR125W';'YMR126C';'YMR127C';'YMR129W';'YMR130W';'YMR132C';'YMR133W';'YMR135C';'YFL006W';'YFL011W';'YFL015C';'YFL018C';'YFL020C';'YFL021W';'YFL023W';'YFL025C';'YFL026W';'YFL027C';'YFL028C';'YFL030W';'YFL031W';'YFL032W';'YFL034W';'YFL035C-B';'YFL036W';'YFL040W';'YFL041W';'YFL043C';'YFL044C';'YFL046W';'YFL047W';'YFL048C';'YFL049W';'YFL050C';'YFL051C';'YFL052W';'YFL053W';'YFL054C';'YFL055W';'YFL056C';'YFR001W';'YFR006W';'YFR007W';'YFR008W';'YFR009W';'YFR010W';'YFR012W';'YFR014C';'YFR015C';'YFR016C';'YFR017C';'YFR018C';'YFR020W';'YFR021W';'YFR022W';'YFR023W';'YFR024C-A';'YFR026C';'YFR031C-A';'YKL096W';'YKL097C';'YKL098W';'YKL100C';'YKL101W';'YKL102C';'YKL103C';'YKL105C';'YKL106W';'YKL107W';'YKL109W';'YKL110C';'YKL113C';'YKL114C';'YKL116C';'YKL117W';'YKL118W';'YKL119C';'YKL120W';'YKL121W';'YKL123W';'YKL124W';'YKL126W';'YKL127W';'YKL128C';'YKL129C';'YKL130C';'YKL131W';'YKL132C';'YKL133C';'YKL134C';'YKL135C';'YKL136W';'YKL137W';'YKL138C';'YKL140W';'YKL142W';'YKL143W';'YKL146W';'YKL147C';'YKL148C';'YKL149C';'YKL150W';'YKL151C';'YKL155C';'YKL156W';'YKL157W';'YKL158W';'YKL159C';'YKL160W';'YKL161C';'YKL162C';'YKL163W';'YKL164C';'YKL166C';'YKL167C';'YKL168C';'YKL169C';'YKL170W';'YKL171W';'YKL174C';'YKL175W';'YKL176C';'YKL177W';'YKL178C';'YKL179C';'YKL183W';'YKL184W';'YKL185W';'YKL187C';'YKL188C';'YKL190W';'YKL191W';'YKL197C';'YKL198C';'YKL199C';'YKL200C';'YKL205W';'YKL206C';'YKL207W';'YKL208W';'YKL211C';'YKL212W';'YKL213C';'YKL214C';'YKL216W';'YKL217W';'YKL218C';'YKL221W';'YKL222C';'YKR001C';'YKR003W';'YKR005C';'YKR006C';'YKR007W';'YKR009C';'YKR011C';'YKR012C';'YKR013W';'YKR014C';'YKR015C';'YKR016W';'YKR017C';'YKR018C';'YKR020W';'YKR021W';'YKR024C';'YKR026C';'YKR030W';'YKR031C';'YKR032W';'YKR033C';'YKR035C';'YKR042W';'YKR043C';'YKR044W';'YKR045C';'YKR047W';'YKR048C';'YKR049C';'YKR050W';'YKR051W';'YKR052C';'YKR054C';'YKR055W';'YKR056W';'YKR057W';'YKR058W';'YKR059W';'YKR060W';'YKR061W';'YKR064W';'YKR065C';'YLR228C';'YLR231C';'YLR232W';'YLR233C';'YLR234W';'YLR235C';'YLR236C';'YLR238W';'YLR239C';'YLR240W';'YLR241W';'YLR242C';'YLR244C';'YLR247C';'YLR248W';'YLR250W';'YLR251W';'YLR252W';'YLR253W';'YLR254C';'YLR255C';'YLR257W';'YLR258W';'YLR260W';'YLR261C';'YLR262C';'YLR263W';'YLR264W';'YLR265C';'YLR266C';'YLR267W';'YLR268W';'YLR269C';'YLR270W';'YLR271W';'YLR273C';'YLR278C';'YLR279W';'YLR280C';'YLR281C';'YLR282C';'YLR283W';'YLR284C';'YLR285W';'YLR287C';'YLR287C-A';'YLR288C';'YLR289W';'YLR290C';'YLR292C';'YLR294C';'YLR295C';'YLR296W';'YLR297W';'YLR299W';'YLR300W';'YLR303W';'YLR304C';'YLR306W';'YLR307W';'YLR308W';'YLR309C';'YLR311C';'YLR312C';'YLR312W-A';'YLR313C';'YLR315W';'YLR318W';'YLR319C';'YLR320W';'YLR322W';'YLR324W';'YLR325C';'YLR326W';'YLR327C';'YLR328W';'YLR329W';'YLR330W';'YLR331C';'YLR332W';'YLR333C';'YLR335W';'YLR337C';'YLR338W';'YLR341W';'YLR342W';'YLR344W';'YLR345W';'YLR348C';'YLR349W';'YLR350W';'YLR351C';'YLR352W';'YLR353W';'YLR354C';'YLR356W';'YLR357W';'YLR360W';'YLR362W';'YLR363C';'YLR364W';'YLR365W';'YLR366W';'YLR367W';'YLR368W';'YLR369W';'YLR371W';'YLR372W';'YLR373C';'YLR374C';'YLR375W';'YLR376C';'YLR377C';'YLR380W';'YLR381W';'YLR384C';'YLR385C';'YLR386W';'YLR387C';'YLR388W';'YLR389C';'YLR390W';'YLR391W';'YLR392C';'YLR393W';'YLR395C';'YLR396C';'YLR398C';'YLR399C';'YLR400W';'YLR401C';'YLR402W';'YLR403W';'YLR404W';'YLR405W';'YLR407W';'YLR408C';'YLR410W';'YLR412W';'YLR413W';'YLR414C';'YLR415C';'YLR416C';'YLR417W';'YLR418C';'YLR421C';'YML089C';'YML088W';'YML087C';'YML086C';'YML084W';'YML083C';'YML082W';'YML081W';'YML080W';'YML079W';'YML078W';'YML063W';'YML062C';'YML061C';'YML060W';'YML059C';'YML058W';'YML057W';'YML058C-A';'YML056C';'YML055W';'YML054C';'YML053C';'YML052W';'YML051W';'YML050W';'YML048W';'YML048W-A';'YNL001W';'YNL003C';'YNL004W';'YNL005C';'YNL008C';'YNL009W';'YNL010W';'YNL012W';'YNL013C';'YNL015W';'YNL016W';'YNL020C';'YNL021W';'YNL022C';'YNL023C';'YNL024C';'YNL025C';'YNL027W';'YNL028W';'YNL029C';'YNL030W';'YNL031C';'YNL032W';'YNL034W';'YNL035C';'YNL037C';'YNL040W';'YNL041C';'YNL043C';'YNL044W';'YNL045W';'YNL046W';'YNL049C';'YNL050C';'YNR001C';'YNR002C';'YNR004W';'YNR005C';'YNR006W';'YNR007C';'YNR008W';'YNR009W';'YNR010W';'YNR012W';'YNR013C';'YNR014W';'YNR015W';'YNR018W';'YNR019W';'YNR020C';'YNR021W';'YNR022C';'YNR024W';'YNR025C';'YNR027W';'YNR028W';'YNR029C';'YNR030W';'YNR031C';'YNR032W';'YNR034W';'YNR036C';'YNR037C';'YNR039C';'YNR040W';'YNR041C';'YNR042W';'YNR045W';'YNR047W';'YNR048W';'YNR049C';'YPR006C';'YPR009W';'YPR012W';'YPR014C';'YPR015C';'YPR017C';'YPR018W';'YPR020W';'YPR027C';'YPR028W';'YPR029C';'YPR030W';'YPR032W';'YPR036W';'YPR038W';'YPR039W';'YPR040W';'YPR042C';'YPR044C';'YPR045C';'YPR046W';'YPR047W';'YPR049C';'YPR051W';'YPR052C';'YPR053C';'YPR054W';'YPR057W';'YPR058W';'YPR059C';'YPR060C';'YPR061C';'YPR062W';'YPR063C';'YPR065W';'YPR066W';'YPR068C';'YPR069C';'YPR070W';'YPR071W';'YPR072W';'YPR073C';'YPR074C';'YPR075C';'YPR076W';'YPR077C';'YPR079W';'YPR084W';'YPR087W';'YPR089W';'YPR090W';'YPR092W';'YPR093C';'YPR095C';'YPR096C';'YPR098C';'YPR097W';'YPR099C';'YPR100W';'YPR101W';'YPR106W';'YPR109W';'YPR111W';'YPR114W';'YPR115W';'YPR116W';'YPR117W';'YPR119W';'YPR120C';'YPR121W';'YPR122W';'YPR123C';'YPR124W';'YPR125W';'YPR126C';'YPR127W';'YPR128C';'YPR129W';'YPR130C';'YPR131C';'YPR132W';'YPR134W';'YPR135W';'YPR138C';'YPR139C';'YPR140W';'YPR141C';'YPR145W';'YPR146C';'YPR147C';'YPR148C';'YPR149W';'YPR150W';'YPR152C';'YPR153W';'YPR154W';'YPR155C';'YPR156C';'YPR157W';'YPR158W';'YPR159W';'YPR160W';'YPR163C';'YPR164W';'YPR166C';'YPR167C';'YPR170C';'YPR171W';'YPR172W';'YPR173C';'YPR174C';'YPR179C';'YPR184W';'YPR185W';'YPR188C';'YPR189W';'YPR191W';'YPR192W';'YPR193C';'YPR194C';'YPR195C';'YPR196W';'YPR197C';'YPR198W';'YPR199C';'YPR200C';'YPR201W';'YKL001C';'YKL002W';'YKL003C';'YKL006W';'YKL007W';'YKL008C';'YKL009W';'YKL010C';'YKL011C';'YKL015W';'YKL016C';'YKL017C';'YKL020C';'YKL023W';'YKL025C';'YKL026C';'YKL027W';'YKL029C';'YKL031W';'YKL032C';'YKL034W';'YKL037W';'YKL038W';'YKL039W';'YKL040C';'YKL041W';'YKL043W';'YKL044W';'YKL046C';'YKL047W';'YKL048C';'YKL050C';'YKL051W';'YKL053W';'YKL054C';'YKL055C';'YKL056C';'YKL057C';'YKL061W';'YKL062W';'YKL063C';'YKL064W';'YKL065C';'YKL066W';'YKL067W';'YKL068W';'YKL069W';'YKL070W';'YKL071W';'YKL072W';'YKL073W';'YKL074C';'YKL075C';'YKL076C';'YKL077W';'YKL079W';'YKL080W';'YKL081W';'YKL084W';'YKL085W';'YKL086W';'YKL087C';'YKL090W';'YKL091C';'YKL092C';'YKL093W';'YKL094W';'YLR042C';'YLR043C';'YLR044C';'YLR046C';'YLR047C';'YLR048W';'YLR049C';'YLR053C';'YLR054C';'YLR055C';'YLR056W';'YLR057W';'YLR058C';'YLR059C';'YLR061W';'YLR062C';'YLR063W';'YLR064W';'YLR065C';'YLR067C';'YLR068W';'YLR069C';'YLR070C';'YLR072W';'YLR073C';'YLR074C';'YLR077W';'YLR079W';'YLR080W';'YLR081W';'YLR082C';'YLR083C';'YLR084C';'YLR085C';'YLR087C';'YLR089C';'YLR090W';'YLR091W';'YLR092W';'YLR093C';'YLR094C';'YLR095C';'YLR096W';'YLR097C';'YLR098C';'YLR099C';'YLR102C';'YLR104W';'YLR107W';'YLR108C';'YLR109W';'YLR111W';'YLR112W';'YLR113W';'YLR114C';'YLR118C';'YLR119W';'YLR120C';'YLR121C';'YLR122C';'YLR123C';'YLR124W';'YLR125W';'YLR420W';'YLR451W';'YLR126C';'YLR128W';'YLR130C';'YLR131C';'YLR133W';'YLR134W';'YLR135W';'YLR136C';'YLR137W';'YLR138W';'YLR139C';'YLR142W';'YLR143W';'YLR144C';'YLR148W';'YLR149C';'YLR150W';'YLR151C';'YLR152C';'YLR154C';'YLR164W';'YLR165C';'YLR168C';'YLR169W';'YLR170C';'YLR171W';'YLR172C';'YLR173W';'YLR174W';'YLR176C';'YLR177W';'YLR178C';'YLR179C';'YLR180W';'YLR181C';'YLR182W';'YLR183C';'YLR184W';'YLR185W';'YLR187W';'YLR188W';'YLR189C';'YLR190W';'YLR191W';'YLR193C';'YLR194C';'YLR199C';'YLR200W';'YLR201C';'YLR202C';'YLR203C';'YLR204W';'YLR205C';'YLR206W';'YLR207W';'YLR209C';'YLR210W';'YLR211C';'YLR213C';'YLR214W';'YLR216C';'YLR217W';'YLR218C';'YLR219W';'YLR220W';'YLR221C';'YLR224W';'YLR225C';'YLR226W';'YLR227C';'YNL146W';'YNL145W';'YNL144C';'YNL143C';'YNL141W';'YNL139C';'YNL136W';'YNL135C';'YNL134C';'YNL133C';'YNL129W';'YNL128W';'YNL127W';'YNL123W';'YNL122C';'YNL117W';'YNL116W';'YNL108C';'YNL101W';'YNL098C';'YNL094W';'YNL092W';'YNL089C';'YNL084C';'YNL081C';'YNL077W';'YNL069C';'YNL064C';'YNL063W';'YNL057W';'YNL058C';'YNL054W';'YPL091W';'YPL090C';'YPL089C';'YPL088W';'YPL087W';'YPL086C';'YPL084W';'YPL081W';'YPL080C';'YPL079W';'YPL078C';'YPL077C';'YPL074W';'YPL072W';'YPL073C';'YPL071C';'YPL070W';'YPL069C';'YPL068C';'YPL067C';'YPL066W';'YPL065W';'YPL064C';'YPL062W';'YPL061W';'YPL060W';'YPL059W';'YPL058C';'YPL057C';'YPL056C';'YPL055C';'YPL054W';'YPL053C';'YPL052W';'YPL051W';'YPL050C';'YPL049C';'YPL048W';'YPL047W';'YPL046C';'YPL045W';'YPL042C';'YPL041C';'YPL040C';'YPL039W';'YPL038W';'YPL037C';'YPL035C';'YPL033C';'YPL032C';'YPL031C';'YPL030W';'YPL029W';'YPL026C';'YPL025C';'YPL024W';'YPL023C';'YPL022W';'YPL021W';'YPL019C';'YPL018W';'YPL015C';'YPL014W';'YPL013C';'YPL009C';'YPL008W';'YPL006W';'YPL005W';'YPL003W';'YPL002C';'YPL001W';'YPR001W';'YPR002W';'YPR003C';'YPR004C';'YPR005C';'YBL001C';'YBL002W';'YBL003C';'YBL005W';'YBL006C';'YBL007C';'YBL008W';'YBL009W';'YBL010C';'YBL011W';'YBL012C';'YBL013W';'YBL015W';'YBL016W';'YBL017C';'YBL019W';'YBL021C';'YBL022C';'YBL024W';'YBL025W';'YBL027W';'YBL028C';'YBL029W';'YBL031W';'YBL032W';'YBL033C';'YBL036C';'YBL037W';'YBL038W';'YBL039C';'YBL042C';'YBL043W';'YBL044W';'YBL045C';'YBL046W';'YBL047C';'YBL048W';'YBL049W';'YBL051C';'YBL052C';'YBL053W';'YBL054W';'YBL055C';'YBL056W';'YBL057C';'YBL058W';'YBL059W';'YBL060W';'YBL061C';'YBL062W';'YBL063W';'YBL064C';'YBL065W';'YBL066C';'YBL067C';'YBL068W';'YBL069W';'YBL070C';'YBL071C';'YBL072C';'YBL075C';'YBL078C';'YBL079W';'YBL080C';'YBL081W';'YBL082C';'YBL083C';'YBL085W';'YBL086C';'YBL087C';'YBL088C';'YBL089W';'YBL090W';'YBL091C';'YBL093C';'YBL094C';'YBL095W';'YBL096C';'YBL098W';'YBL099W';'YBL100C';'YBL101C';'YBL102W';'YBL103C';'YBL104C';'YBL106C';'YBL107C';'YBR001C';'YBR003W';'YBR005W';'YBR006W';'YBR007C';'YBR008C';'YBR009C';'YBR010W';'YBR012C';'YBR013C';'YBR014C';'YBR015C';'YBR016W';'YBR018C';'YBR019C';'YBR020W';'YBR021W';'YBR022W';'YBR023C';'YBR024W';'YBR025C';'YBR026C';'YBR027C';'YBR028C';'YBR030W';'YBR031W';'YBR032W';'YBR033W';'YBR034C';'YBR035C';'YBR036C';'YBR037C';'YBR040W';'YBR041W';'YBR042C';'YBR043C';'YBR044C';'YBR045C';'YBR046C';'YBR047W';'YBR048W';'YBR050C';'YBR051W';'YBR052C';'YBR053C';'YBR054W';'YBR056W';'YBR057C';'YBR058C';'YBR059C';'YBR061C';'YBR062C';'YBR063C';'YBR064W';'YBR065C';'YBR066C';'YBR067C';'YBR068C';'YBR069C';'YBR071W';'YBR072W';'YBR073W';'YBR074W';'YBR075W';'YBR076W';'YBR077C';'YBR078W';'YBR081C';'YBR082C';'YBR084W';'YBR085W';'YBR090C-A';'YBR092C';'YBR093C';'YBR094W';'YBR095C';'YBR097W';'YBR098W';'YBR099C';'YBR100W';'YBR101C';'YBR103W';'YBR104W';'YBR105C';'YBR106W';'YBR107C';'YBR108W';'YBR111C';'YBR113W';'YBR114W';'YBR115C';'YBR116C';'YBR119W';'YBR120C';'YBR121C';'YBR122C';'YBR126C';'YBR127C';'YBR128C';'YBR129C';'YBR130C';'YBR133C';'YBR134W';'YBR137W';'YBR138C';'YBR139W';'YBR141C';'YBR144C';'YBR145W';'YBR146W';'YBR147W';'YBR148W';'YBR149W';'YBR151W';'YBR156C';'YBR157C';'YBR158W';'YBR159W';'YBR161W';'YBR162C';'YBR162W-A';'YBR163W';'YBR164C';'YBR165W';'YBR166C';'YBR170C';'YBR171W';'YBR172C';'YBR173C';'YBR174C';'YBR175W';'YBR176W';'YBR177C';'YBR178W';'YBR179C';'YBR180W';'YBR181C';'YBR182C';'YBR183W';'YBR184W';'YBR185C';'YBR186W';'YBR187W';'YBR188C';'YBR194W';'YBR195C';'YBR197C';'YBR199W';'YBR200W';'YBR201W';'YBR203W';'YBR204C';'YBR205W';'YBR206W';'YBR207W';'YBR208C';'YBR209W';'YBR210W';'YBR212W';'YBR213W';'YBR214W';'YBR215W';'YBR216C';'YBR217W';'YBR218C';'YBR219C';'YBR220C';'YBR221C';'YBR222C';'YBR223C';'YBR224W';'YBR225W';'YBR226C';'YBR227C';'YBR228W';'YBR229C';'YBR230C';'YBR231C';'YBR233W';'YBR235W';'YBR238C';'YBR239C';'YBR240C';'YBR241C';'YBR242W';'YBR244W';'YBR245C';'YBR246W';'YBR248C';'YBR249C';'YBR250W';'YBR251W';'YBR255W';'YBR258C';'YBR259W';'YBR260C';'YBR261C';'YBR262C';'YBR263W';'YBR264C';'YBR266C';'YBR267W';'YBR268W';'YCL001W';'YCL002C';'YCL005W';'YCL006C';'YCL007C';'YCL008C';'YCL009C';'YCL010C';'YCL011C';'YCL012W';'YCL013W';'YCL014W';'YCL016C';'YCL023C';'YCL024W';'YCL025C';'YCL026C';'YCL027W';'YCL028W';'YCL029C';'YCL030C';'YCL032W';'YCL033C';'YCL034W';'YCL036W';'YCL037C';'YCL039W';'YCL040W';'YCL042W';'YCL044C';'YCL045C';'YCL046W';'YCL047C';'YCL048W';'YCL049C';'YCL050C';'YCL051W';'YCL055W';'YCL056C';'YCL057W';'YCL060C';'YCL061C';'YCL062W';'YCL063W';'YCL064C';'YCL069W';'YCR001W';'YCR002C';'YCR003W';'YCR004C';'YCR005C';'YCR006C';'YCR007C';'YCR008W';'YCR009C';'YCR010C';'YCR011C';'YCR014C';'YCR015C';'YCR016W';'YCR017C';'YCR019W';'YCR020C';'YCR020C-A';'YCR021C';'YCR022C';'YCR023C';'YDL001W';'YDL002C';'YDL005C';'YDL006W';'YDL009C';'YDL010W';'YDL011C';'YDL012C';'YDL013W';'YDL018C';'YDL019C';'YDL020C';'YDL021W';'YDL022W';'YDL023C';'YDL024C';'YDL025C';'YDL026W';'YDL027C';'YDL032W';'YDL033C';'YDL034W';'YDL035C';'YDL036C';'YDL037C';'YDL038C';'YDL039C';'YDL040C';'YDL041W';'YDL042C';'YDL044C';'YDL045W-A';'YDL046W';'YDL047W';'YDL048C';'YDL049C';'YDL050C';'YDL051W';'YDL052C';'YDL053C';'YDL054C';'YDL056W';'YDL057W';'YDL059C';'YDL061C';'YDL062W';'YDL063C';'YDL065C';'YDL066W';'YDL067C';'YDL068W';'YDL069C';'YDL070W';'YDL071C';'YDL072C';'YDL073W';'YDL074C';'YDL075W';'YDL076C';'YDL077C';'YDL078C';'YDL079C';'YDL080C';'YDL081C';'YDL082W';'YDL083C';'YDL085W';'YDL086W';'YDL088C';'YDL089W';'YDL090C';'YDL091C';'YDL093W';'YDL094C';'YDL095W';'YDL096C';'YDL099W';'YDL100C';'YDL101C';'YDL104C';'YDL106C';'YDL107W';'YDL109C';'YDL110C';'YDL112W';'YDL113C';'YDL114W';'YDL115C';'YDL116W';'YDL117W';'YDL118W';'YDL119C';'YDL121C';'YDL122W';'YDL123W';'YDL124W';'YDL125C';'YDL127W';'YDL128W';'YDL129W';'YDL130W';'YDL131W';'YDL133W';'YDL134C';'YDL134C-A';'YDL135C';'YDL136W';'YDL137W';'YDL138W';'YDL142C';'YDL144C';'YDL146W';'YDL149W';'YDL151C';'YDL154W';'YDL155W';'YDL156W';'YDL157C';'YDL159W';'YDL160C';'YDL161W';'YDL162C';'YDL167C';'YDL168W';'YDL169C';'YDL170W';'YDL171C';'YDL172C';'YDL173W';'YDL174C';'YDL175C';'YDL176W';'YDL177C';'YDL178W';'YDL179W';'YDL180W';'YDL181W';'YDL182W';'YDL183C';'YDL184C';'YDL185W';'YDL186W';'YDL187C';'YDL188C';'YDL189W';'YDL190C';'YDL191W';'YDL192W';'YDL197C';'YDL198C';'YDL199C';'YDL200C';'YDL201W';'YDL202W';'YDL203C';'YDL204W';'YDL206W';'YDL210W';'YDL211C';'YDL213C';'YDL214C';'YDL215C';'YDL216C';'YDL218W';'YDL219W';'YDL222C';'YDL223C';'YDL224C';'YDL225W';'YDL226C';'YDL227C';'YDL229W';'YDL230W';'YDL231C';'YDL232W';'YDL233W';'YDL234C';'YDL236W';'YDL237W';'YDL238C';'YDL239C';'YDL240W';'YDL241W';'YDL242W';'YDL243C';'YDR001C';'YDR003W';'YDR004W';'YDR005C';'YDR006C';'YDR008C';'YDR009W';'YDR010C';'YDR011W';'YDR014W';'YDR015C';'YDR017C';'YDR018C';'YDR019C';'YDR020C';'YDR022C';'YDR024W';'YDR025W';'YDR026C';'YDR027C';'YDR028C';'YDR029W';'YDR030C';'YDR031W';'YDR032C';'YDR033W';'YDR034C';'YDR035W';'YDR036C';'YDR042C';'YDR043C';'YDR046C';'YDR049W';'YDR050C';'YDR051C';'YDR055W';'YDR056C';'YDR057W';'YDR059C';'YDR061W';'YDR063W';'YDR065W';'YDR066C';'YDR067C';'YDR068W';'YDR069C';'YDR070C';'YDR072C';'YDR073W';'YDR075W';'YDR076W';'YDR077W';'YDR078C';'YDR079W';'YDR080W';'YDR083W';'YDR084C';'YDR085C';'YDR089W';'YDR090C';'YDR092W';'YDR093W';'YDR094W';'YDR095C';'YDR096W';'YDR097C';'YDR098C';'YDR099W';'YDR100W';'YDR101C';'YDR102C';'YDR103W';'YDR104C';'YDR105C';'YDR107C';'YDR108W';'YDR109C';'YDR110W';'YDR111C';'YDR112W';'YDR114C';'YDR115W';'YDR116C';'YDR117C';'YDR119W';'YDR120C';'YDR121W';'YDR122W';'YDR123C';'YDR124W';'YDR125C';'YDR126W';'YDR127W';'YDR128W';'YDR129C';'YDR130C';'YDR131C';'YDR132C';'YDR133C';'YDR134C';'YDR135C';'YDR136C';'YDR137W';'YDR138W';'YDR139C';'YDR140W';'YDR142C';'YDR143C';'YDR144C';'YDR146C';'YDR147W';'YDR148C';'YDR149C';'YDR150W';'YDR151C';'YDR152W';'YDR153C';'YDR154C';'YDR155C';'YDR156W';'YDR157W';'YDR158W';'YDR159W';'YDR161W';'YDR162C';'YDR163W';'YDR165W';'YDR169C';'YDR171W';'YDR173C';'YDR175C';'YDR176W';'YDR178W';'YDR179C';'YDR179W-A';'YDR181C';'YDR183W';'YDR184C';'YDR185C';'YDR186C';'YDR191W';'YDR192C';'YDR193W';'YDR194C';'YDR195W';'YDR197W';'YDR198C';'YDR199W';'YDR200C';'YDR203W';'YDR204W';'YDR206W';'YDR207C';'YDR209C';'YDR210W';'YDR213W';'YDR214W';'YDR215C';'YDR216W';'YDR217C';'YDR218C';'YDR219C';'YDR220C';'YDR221W';'YDR222W';'YDR223W';'YDR225W';'YDR226W';'YDR227W';'YDR229W';'YDR230W';'YDR231C';'YDR233C';'YDR234W';'YDR237W';'YDR239C';'YDR241W';'YDR244W';'YDR245W';'YDR247W';'YDR248C';'YDR249C';'YDR250C';'YDR251W';'YDR252W';'YDR253C';'YDR254W';'YDR255C';'YDR256C';'YDR257C';'YDR258C';'YDR259C';'YDR260C';'YDR261C';'YDR262W';'YDR263C';'YDR264C';'YDR265W';'YDR266C';'YDR268W';'YDR270W';'YDR272W';'YDR273W';'YDR274C';'YDR275W';'YDR276C';'YDR277C';'YDR278C';'YDR279W';'YDR281C';'YDR282C';'YDR283C';'YDR284C';'YDR285W';'YDR286C';'YDR287W';'YDR289C';'YDR291W';'YDR293C';'YDR294C';'YDR295C';'YDR296W';'YDR297W';'YDR298C';'YDR300C';'YDR304C';'YDR305C';'YDR306C';'YDR307W';'YDR309C';'YDR310C';'YDR312W';'YDR313C';'YDR314C';'YDR315C';'YDR316W';'YDR317W';'YDR318W';'YDR319C';'YDR320C';'YDR321W';'YDR322W';'YDR323C';'YDR329C';'YDR330W';'YDR332W';'YDR333C';'YDR334W';'YDR335W';'YDR336W';'YDR337W';'YDR338C';'YDR340W';'YDR344C';'YDR345C';'YDR346C';'YDR347W';'YDR348C';'YDR349C';'YDR350C';'YDR351W';'YDR352W';'YDR354W';'YDR357C';'YDR358W';'YDR359C';'YDR360W';'YDR363W';'YDR364C';'YDR368W';'YDR369C';'YDR370C';'YDR371W';'YDR372C';'YDR374C';'YDR375C';'YDR377W';'YDR378C';'YDR379W';'YDR380W';'YDR382W';'YDR383C';'YDR384C';'YDR385W';'YDR386W';'YDR387C';'YDR388W';'YDR389W';'YDR391C';'YDR392W';'YDR393W';'YDR395W';'YDR399W';'YDR400W';'YDR401W';'YDR402C';'YDR403W';'YDR405W';'YDR406W';'YDR408C';'YDR409W';'YDR410C';'YDR411C';'YDR414C';'YDR415C';'YDR418W';'YDR419W';'YDR420W';'YDR421W';'YDR422C';'YDR423C';'YDR425W';'YDR426C';'YDR428C';'YDR430C';'YDR431W';'YDR432W';'YDR433W';'YDR435C';'YDR436W';'YDR438W';'YDR439W';'YDR440W';'YDR441C';'YDR442W';'YDR443C';'YDR446W';'YDR447C';'YDR448W';'YDR450W';'YDR451C';'YDR452W';'YDR453C';'YDR455C';'YDR456W';'YDR457W';'YDR458C';'YDR459C';'YDR462W';'YDR463W';'YDR465C';'YDR466W';'YDR467C';'YDR469W';'YDR470C';'YDR471W';'YDR474C';'YDR475C';'YDR476C';'YDR477W';'YDR479C';'YDR480W';'YDR481C';'YDR482C';'YDR483W';'YDR484W';'YDR485C';'YDR486C';'YDR488C';'YDR490C';'YDR491C';'YDR492W';'YDR494W';'YDR495C';'YDR496C';'YDR497C';'YDR500C';'YDR501W';'YDR503C';'YDR504C';'YDR505C';'YDR506C';'YDR507C';'YDR508C';'YDR509W';'YDR511W';'YDR512C';'YDR513W';'YDR514C';'YDR516C';'YDR517W';'YDR518W';'YDR519W';'YDR521W';'YDR522C';'YDR523C';'YDR524C';'YDR525W';'YDR528W';'YDR529C';'YDR530C';'YDR532C';'YDR533C';'YDR534C';'YGL002W';'YGL003C';'YGL004C';'YGL005C';'YGL006W';'YGL007W';'YGL009C';'YGL010W';'YGL012W';'YGL013C';'YGL014W';'YGL015C';'YGL016W';'YGL017W';'YGL019W';'YGL020C';'YGL021W';'YGL023C';'YGL024W';'YGL025C';'YGL026C';'YGL027C';'YGL028C';'YGL029W';'YGL031C';'YGL032C';'YGL033W';'YGL034C';'YGL035C';'YGL036W';'YGL037C';'YGL038C';'YGL039W';'YGL041C';'YGL042C';'YGL043W';'YGL045W';'YGL046W';'YGL049C';'YGL050W';'YGL051W';'YGL053W';'YGL054C';'YGL056C';'YGL057C';'YGL058W';'YGL059W';'YGL060W';'YGL062W';'YGL063W';'YGL064C';'YGL066W';'YGL067W';'YGL070C';'YGL071W';'YGL072C';'YGL076C';'YGL077C';'YGL078C';'YGL079W';'YGL080W';'YGL081W';'YGL082W';'YGL083W';'YGL084C';'YGL085W';'YGL086W';'YGL087C';'YGL088W';'YGL089C';'YGL090W';'YGL094C';'YGL095C';'YGL096W';'YGL101W';'YGL104C';'YGL105W';'YGL107C';'YGL108C';'YGL109W';'YGL110C';'YGL114W';'YGL115W';'YGL117W';'YGL118C';'YGL121C';'YGL124C';'YGL125W';'YGL126W';'YGL127C';'YGL129C';'YGL131C';'YGL132W';'YGL133W';'YGL135W';'YGL136C';'YGL138C';'YGL139W';'YGL140C';'YGL141W';'YGL143C';'YGL144C';'YGL146C';'YGL147C';'YGL148W';'YGL149W';'YGL151W';'YGL152C';'YGL153W';'YGL154C';'YGL156W';'YGL157W';'YGL158W';'YGL159W';'YGL160W';'YGL161C';'YGL162W';'YGL163C';'YGL164C';'YGL165C';'YGL166W';'YGL167C';'YGL168W';'YGL170C';'YGL173C';'YGL174W';'YGL175C';'YGL176C';'YGL177W';'YGL179C';'YGL180W';'YGL181W';'YGL194C';'YGL195W';'YGL196W';'YGL197W';'YGL198W';'YGL200C';'YGL202W';'YGL203C';'YGL205W';'YGL206C';'YGL208W';'YGL209W';'YGL210W';'YGL211W';'YGL212W';'YGL213C';'YGL215W';'YGL216W';'YGL220W';'YGL221C';'YGL222C';'YGL223C';'YGL224C';'YGL226C-A';'YGL226W';'YGL227W';'YGL228W';'YGL229C';'YGL230C';'YGL231C';'YGL232W';'YGL234W';'YGL236C';'YGL237C';'YGL240W';'YGL241W';'YGL242C';'YGL243W';'YGL244W';'YGL246C';'YGL248W';'YGL249W';'YGL250W';'YGL251C';'YGL252C';'YGL253W';'YGL254W';'YGL255W';'YGL256W';'YGL257C';'YGL258W';'YGL259W';'YGL260W';'YGL261C';'YGL262W';'YGL263W';'YGR001C';'YGR003W';'YGR004W';'YGR006W';'YGR007W';'YGR008C';'YGR010W';'YGR012W';'YGR014W';'YGR015C';'YGR016W';'YGR017W';'YGR019W';'YGR020C';'YGR021W';'YGR023W';'YGR026W';'YGR027C';'YGR031W';'YGR033C';'YGR034W';'YGR035C';'YGR036C';'YGR037C';'YGR039W';'YGR041W';'YGR042W';'YGR043C';'YGR044C';'YGR045C';'YGR049W';'YGR051C';'YGR052W';'YGR054W';'YGR055W';'YGR056W';'YGR057C';'YGR058W';'YGR059W';'YGR061C';'YGR062C';'YGR064W';'YGR066C';'YGR067C';'YGR068C';'YGR069W';'YGR070W';'YGR071C';'YGR072W';'YGR076C';'YGR077C';'YGR078C';'YGR079W';'YGR080W';'YGR081C';'YGR084C';'YGR085C';'YGR087C';'YGR088W';'YGR096W';'YGR097W';'YGR100W';'YGR101W';'YGR102C';'YGR104C';'YGR105W';'YGR107W';'YGR108W';'YGR109C';'YGR111W';'YGR112W';'YGR118W';'YGR121C';'YGR122W';'YGR123C';'YGR124W';'YGR125W';'YGR126W';'YGR127W';'YGR129W';'YGR130C';'YGR131W';'YGR132C';'YGR133W';'YGR135W';'YGR136W';'YGR137W';'YGR138C';'YGR139W';'YGR141W';'YGR142W';'YGR143W';'YGR144W';'YGR146C';'YGR148C';'YGR149W';'YGR150C';'YGR151C';'YGR152C';'YGR153W';'YGR154C';'YGR157W';'YGR159C';'YGR160W';'YGR161C';'YGR163W';'YGR164W';'YGR165W';'YGR166W';'YGR167W';'YGR168C';'YGR169C';'YGR170W';'YGR171C';'YGR173W';'YGR174C';'YGR176W';'YGR177C';'YGR178C';'YGR180C';'YGR181W';'YGR182C';'YGR183C';'YGR184C';'YGR187C';'YGR189C';'YGR192C';'YGR193C';'YGR194C';'YGR196C';'YGR197C';'YGR199W';'YGR200C';'YGR202C';'YGR203W';'YGR205W';'YGR206W';'YGR207C';'YGR208W';'YGR209C';'YGR212W';'YGR213C';'YGR214W';'YGR215W';'YGR217W';'YHR142W';'YHR143W';'YHR147C';'YHR150W';'YHR151C';'YHR152W';'YHR153C';'YHR154W';'YHR155W';'YHR156C';'YHR157W';'YHR158C';'YHR159W';'YHR160C';'YHR161C';'YHR163W';'YHR167W';'YHR176W';'YHR177W';'YHR178W';'YHR179W';'YHR182W';'YHR183W';'YHR184W';'YHR189W';'YHR195W';'YHR198C';'YHR199C';'YHR200W';'YHR202W';'YHR203C';'YHR204W';'YHR206W';'YHR207C';'YHR209W';'YHR210C';'YMR038C';'YOR147W';'YOR179C';'YOR180C';'YER014W';'YER028C';'YER044C';'YER055C';'YER087W';'YIL134W';'YKL135C';'YEL024W';'YLR308W';'YPL006W';'YOR364W';'YOR054c';'YBR132C';'YDR269C';'YDR271C';'YDR290W';'YGL199C';'YGL214W';'YGL217C';'YGL218W';'YGL235W';'YGR011W';'YGR018C';'YGR022C';'YGR025W';'YOL086C';'YHR210C';'YHR209W';'YHR198C';'YOR158W';'YHR177W';};
    genes = upper(genes);
end

function genes = kuepfer_glucose_aerobic
    % list of 59 genes found to be essential in defined media + glucose +
    % aerobic - those scored as 0 in the supplemental data to Doi:
    % 10.1101/gr.3992505

    genes = {'YAL012W';'YAR015W';'YBL033C';'YBR115C';'YBR209W';'YCR053W';'YCR066W';'YDL072C';'YDL075W';'YDR007W';'YDR008C';'YDR127W';'YDR158W';'YDR234W';'YDR300C';'YDR354W';'YDR477W';'YEL024W';'YER052C';'YER057C';'YER068C-A';'YER069W';'YER086W';'YER090W';'YFL032W';'YGL026C';'YGL148W';'YGL154C';'YGL234W';'YGR061C';'YGR155W';'YGR168C';'YGR204W';'YHL003C';'YHL005C';'YHR018C';'YHR025W';'YIL069C';'YIL094C';'YIR034C';'YJL088W';'YJR122W';'YJR139C';'YKL135C';'YKL211C';'YKR058W';'YLR304C';'YML022W';'YMR300C';'YNL218W';'YNL220W';'YNL316C';'YNR050C';'YOL058W';'YOL143C';'YOR128C';'YPL006W';'YPR060C';'YPR067W';};
    genes = upper(genes);
end

function genes = kuepfer_galactose_aerobic
    % list of 307 genes found to be essential in defined media + galactose +
    % aerobic - those scored as 0 in the supplemental data to Doi:
    % 10.1101/gr.3992505

    genes = {'YAL009W';'YAL012W';'YAL039C';'YAL048C';'YAR015W';'YBL002W';'YBL012C';'YBL022C';'YBL033C';'YBL038W';'YBL045C';'YBL090W';'YBL099W';'YBL100C';'YBR003W';'YBR018C';'YBR019C';'YBR020W';'YBR037C';'YBR112C';'YBR115C';'YBR120C';'YBR122C';'YBR132C';'YBR179C';'YBR209W';'YBR251W';'YBR268W';'YBR282W';'YCR003W';'YCR004C';'YCR024C';'YCR028C-A';'YCR046C';'YCR053W';'YCR066W';'YCR071C';'YCR084C';'YDL044C';'YDL045W-A';'YDL049C';'YDL057W';'YDL062W';'YDL063C';'YDL067C';'YDL068W';'YDL069C';'YDL072C';'YDL075W';'YDL104C';'YDL107W';'YDL113C';'YDL135C';'YDL146W';'YDL167C';'YDL181W';'YDL198C';'YDL202W';'YDR007W';'YDR008C';'YDR009W';'YDR042C';'YDR065W';'YDR078C';'YDR079W';'YDR114C';'YDR115W';'YDR127W';'YDR158W';'YDR175C';'YDR194C';'YDR197W';'YDR204W';'YDR230W';'YDR231C';'YDR234W';'YDR237W';'YDR268W';'YDR283C';'YDR295C';'YDR296W';'YDR298C';'YDR300C';'YDR322W';'YDR337W';'YDR354W';'YDR477W';'YDR507C';'YDR518W';'YDR523C';'YDR529C';'YEL024W';'YEL029C';'YEL036C';'YEL050C';'YER017C';'YER050C';'YER052C';'YER057C';'YER058W';'YER068C-A';'YER070W';'YER086W';'YER087W';'YER090W';'YER103W';'YER110C';'YER122C';'YER141W';'YER153C';'YER154W';'YER169W';'YFL016C';'YFL032W';'YFL036W';'YGL012W';'YGL026C';'YGL038C';'YGL064C';'YGL095C';'YGL107C';'YGL129C';'YGL135W';'YGL148W';'YGL154C';'YGL206C';'YGL218W';'YGL220W';'YGL223C';'YGL234W';'YGL240W';'YGR061C';'YGR062C';'YGR076C';'YGR102C';'YGR112W';'YGR150C';'YGR155W';'YGR168C';'YGR171C';'YGR180C';'YGR204W';'YGR215W';'YGR219W';'YGR220C';'YGR222W';'YGR255C';'YGR257C';'YGR262C';'YHL005C';'YHL038C';'YHR011W';'YHR025W';'YHR038W';'YHR051W';'YHR091C';'YHR116W';'YHR120W';'YHR147C';'YHR168W';'YHR177W';'YIL069C';'YIL094C';'YIR021W';'YIR034C';'YJL003W';'YJL063C';'YJL088W';'YJL096W';'YJL102W';'YJL180C';'YJL209W';'YJR004C';'YJR122W';'YJR139C';'YJR144W';'YKL003C';'YKL016C';'YKL134C';'YKL135C';'YKL138C';'YKL155C';'YKL169C';'YKL170W';'YKL194C';'YKL211C';'YKR057W';'YKR058W';'YKR061W';'YKR085C';'YLL018C-A';'YLR067C';'YLR069C';'YLR091W';'YLR114C';'YLR139C';'YLR202C';'YLR203C';'YLR204W';'YLR238W';'YLR255C';'YLR257W';'YLR295C';'YLR304C';'YLR312W-A';'YLR369W';'YLR382C';'YML022W';'YML061C';'YML088W';'YML110C';'YML129C';'YMR021C';'YMR024W';'YMR064W';'YMR066W';'YMR071C';'YMR072W';'YMR084W';'YMR089C';'YMR097C';'YMR098C';'YMR135W-A';'YMR150C';'YMR151W';'YMR158W';'YMR184W';'YMR193W';'YMR228W';'YMR231W';'YMR256C';'YMR257C';'YMR267W';'YMR282C';'YMR286W';'YMR287C';'YMR300C';'YNL003C';'YNL005C';'YNL073W';'YNL081C';'YNL084C';'YNL160W';'YNL170W';'YNL177C';'YNL184C';'YNL213C';'YNL220W';'YNL225C';'YNL252C';'YNL284C';'YNL316C';'YNR036C';'YNR037C';'YNR041C';'YNR050C';'YOL023W';'YOL058W';'YOL095C';'YOL096C';'YOL100W';'YOL143C';'YOR033C';'YOR065W';'YOR125C';'YOR128C';'YOR130C';'YOR150W';'YOR158W';'YOR184W';'YOR187W';'YOR199W';'YOR200W';'YOR201C';'YOR205C';'YOR211C';'YOR241W';'YOR304C-A';'YOR305W';'YOR318C';'YPL005W';'YPL006W';'YPL013C';'YPL029W';'YPL031C';'YPL040C';'YPL059W';'YPL078C';'YPL097W';'YPL104W';'YPL118W';'YPL132W';'YPL148C';'YPL172C';'YPL173W';'YPL215W';'YPL248C';'YPL271W';'YPR047W';'YPR060C';'YPR067W';'YPR072W';'YPR099C';'YPR100W';'YPR116W';'YPR124W';'YPR166C';};
    genes = upper(genes);
end

function genes = kuepfer_glycerol_aerobic
    % list of 291 genes found to be essential in defined media + glycerol +
    % aerobic - those scored as 0 in the supplemental data to Doi:
    % 10.1101/gr.3992505
    
    genes = {'YAL009W';'YAL012W';'YAL039C';'YAL048C';'YAR015W';'YBL002W';'YBL012C';'YBL022C';'YBL033C';'YBL038W';'YBL044W';'YBL045C';'YBL090W';'YBL099W';'YBL100C';'YBR003W';'YBR037C';'YBR112C';'YBR115C';'YBR120C';'YBR122C';'YBR132C';'YBR179C';'YBR209W';'YBR251W';'YBR268W';'YBR282W';'YCR003W';'YCR004C';'YCR024C';'YCR028C-A';'YCR046C';'YCR053W';'YCR066W';'YCR071C';'YDL044C';'YDL045W-A';'YDL049C';'YDL057W';'YDL062W';'YDL063C';'YDL067C';'YDL068W';'YDL069C';'YDL072C';'YDL075W';'YDL104C';'YDL107W';'YDL113C';'YDL135C';'YDL146W';'YDL167C';'YDL181W';'YDL198C';'YDL202W';'YDR007W';'YDR008C';'YDR042C';'YDR065W';'YDR078C';'YDR079W';'YDR114C';'YDR127W';'YDR158W';'YDR194C';'YDR197W';'YDR204W';'YDR230W';'YDR231C';'YDR234W';'YDR237W';'YDR268W';'YDR283C';'YDR295C';'YDR296W';'YDR298C';'YDR300C';'YDR322W';'YDR337W';'YDR354W';'YDR470C';'YDR477W';'YDR507C';'YDR518W';'YDR523C';'YDR529C';'YEL024W';'YEL029C';'YEL036C';'YEL050C';'YER017C';'YER050C';'YER052C';'YER058W';'YER070W';'YER086W';'YER087W';'YER090W';'YER103W';'YER110C';'YER122C';'YER141W';'YER153C';'YER154W';'YER169W';'YFL016C';'YFL032W';'YFL036W';'YGL012W';'YGL026C';'YGL038C';'YGL064C';'YGL095C';'YGL107C';'YGL129C';'YGL135W';'YGL148W';'YGL154C';'YGL218W';'YGL220W';'YGL234W';'YGL240W';'YGR061C';'YGR062C';'YGR076C';'YGR102C';'YGR112W';'YGR150C';'YGR155W';'YGR168C';'YGR171C';'YGR180C';'YGR204W';'YGR215W';'YGR219W';'YGR220C';'YGR222W';'YGR255C';'YGR257C';'YHL005C';'YHL038C';'YHR011W';'YHR025W';'YHR038W';'YHR051W';'YHR091C';'YHR116W';'YHR120W';'YHR147C';'YHR168W';'YHR177W';'YIL094C';'YIL125W';'YIR021W';'YIR034C';'YJL003W';'YJL063C';'YJL088W';'YJL096W';'YJL102W';'YJL180C';'YJL209W';'YJR004C';'YJR122W';'YJR139C';'YJR144W';'YKL003C';'YKL134C';'YKL138C';'YKL155C';'YKL169C';'YKL170W';'YKL194C';'YKL211C';'YKR057W';'YKR058W';'YKR061W';'YKR085C';'YLL018C-A';'YLR067C';'YLR069C';'YLR091W';'YLR114C';'YLR139C';'YLR202C';'YLR203C';'YLR204W';'YLR238W';'YLR255C';'YLR257W';'YLR295C';'YLR304C';'YLR312W-A';'YLR369W';'YLR377C';'YML022W';'YML061C';'YML088W';'YML110C';'YML129C';'YMR024W';'YMR064W';'YMR066W';'YMR071C';'YMR072W';'YMR084W';'YMR089C';'YMR097C';'YMR098C';'YMR135W-A';'YMR150C';'YMR151W';'YMR158W';'YMR184W';'YMR193W';'YMR228W';'YMR256C';'YMR257C';'YMR267W';'YMR282C';'YMR286W';'YMR287C';'YMR300C';'YNL003C';'YNL005C';'YNL073W';'YNL081C';'YNL084C';'YNL160W';'YNL170W';'YNL177C';'YNL184C';'YNL213C';'YNL220W';'YNL225C';'YNL252C';'YNL284C';'YNL316C';'YNR036C';'YNR037C';'YNR041C';'YOL002C';'YOL023W';'YOL058W';'YOL095C';'YOL096C';'YOL100W';'YOL143C';'YOR033C';'YOR065W';'YOR125C';'YOR128C';'YOR130C';'YOR150W';'YOR158W';'YOR184W';'YOR187W';'YOR199W';'YOR200W';'YOR201C';'YOR205C';'YOR211C';'YOR241W';'YOR304C-A';'YOR305W';'YOR318C';'YPL005W';'YPL006W';'YPL013C';'YPL029W';'YPL031C';'YPL040C';'YPL059W';'YPL078C';'YPL097W';'YPL104W';'YPL118W';'YPL132W';'YPL148C';'YPL172C';'YPL173W';'YPL215W';'YPL271W';'YPR047W';'YPR060C';'YPR067W';'YPR072W';'YPR099C';'YPR100W';'YPR124W';'YPR166C';};
    genes = upper(genes);
end

function genes = kuepfer_ethanol_aerobic
    % list of 322 genes found to be essential in defined media + galactose +
    % aerobic - those scored as 0 in the supplemental data to Doi:
    % 10.1101/gr.3992505

    genes = {'YAL009W';'YAL039C';'YAL048C';'YAR015W';'YBL002W';'YBL012C';'YBL021C';'YBL022C';'YBL033C';'YBL038W';'YBL044W';'YBL045C';'YBL090W';'YBL099W';'YBL100C';'YBR003W';'YBR037C';'YBR112C';'YBR115C';'YBR120C';'YBR122C';'YBR132C';'YBR179C';'YBR209W';'YBR251W';'YBR268W';'YBR282W';'YCR003W';'YCR004C';'YCR024C';'YCR028C-A';'YCR046C';'YCR066W';'YCR071C';'YCR084C';'YDL044C';'YDL045W-A';'YDL049C';'YDL057W';'YDL062W';'YDL063C';'YDL067C';'YDL068W';'YDL069C';'YDL072C';'YDL075W';'YDL104C';'YDL107W';'YDL113C';'YDL135C';'YDL146W';'YDL167C';'YDL181W';'YDL198C';'YDL202W';'YDR007W';'YDR008C';'YDR042C';'YDR065W';'YDR078C';'YDR079W';'YDR114C';'YDR115W';'YDR127W';'YDR158W';'YDR175C';'YDR194C';'YDR197W';'YDR204W';'YDR226W';'YDR230W';'YDR231C';'YDR234W';'YDR237W';'YDR268W';'YDR283C';'YDR295C';'YDR296W';'YDR298C';'YDR300C';'YDR322W';'YDR337W';'YDR354W';'YDR377W';'YDR470C';'YDR477W';'YDR507C';'YDR518W';'YDR523C';'YDR529C';'YEL024W';'YEL027W';'YEL029C';'YEL036C';'YEL050C';'YER017C';'YER050C';'YER052C';'YER058W';'YER065C';'YER068C-A';'YER069W';'YER070W';'YER086W';'YER087W';'YER090W';'YER103W';'YER110C';'YER122C';'YER141W';'YER153C';'YER154W';'YER169W';'YFL016C';'YFL032W';'YFL036W';'YGL012W';'YGL026C';'YGL038C';'YGL064C';'YGL095C';'YGL107C';'YGL129C';'YGL135W';'YGL148W';'YGL154C';'YGL206C';'YGL218W';'YGL220W';'YGL223C';'YGL234W';'YGL237C';'YGL240W';'YGR061C';'YGR062C';'YGR076C';'YGR102C';'YGR112W';'YGR150C';'YGR168C';'YGR171C';'YGR180C';'YGR183C';'YGR204W';'YGR208W';'YGR215W';'YGR219W';'YGR220C';'YGR222W';'YGR255C';'YGR257C';'YGR262C';'YHL003C';'YHL005C';'YHL038C';'YHR011W';'YHR018C';'YHR025W';'YHR038W';'YHR039C-B';'YHR051W';'YHR091C';'YHR116W';'YHR120W';'YHR147C';'YHR168W';'YHR177W';'YIL069C';'YIL094C';'YIL125W';'YIR021W';'YIR034C';'YJL003W';'YJL063C';'YJL088W';'YJL096W';'YJL102W';'YJL180C';'YJL209W';'YJR004C';'YJR122W';'YJR139C';'YJR144W';'YKL003C';'YKL016C';'YKL134C';'YKL135C';'YKL138C';'YKL148C';'YKL155C';'YKL169C';'YKL170W';'YKL194C';'YKL211C';'YKR057W';'YKR058W';'YKR061W';'YKR085C';'YKR097W';'YLL018C-A';'YLR067C';'YLR069C';'YLR091W';'YLR114C';'YLR139C';'YLR202C';'YLR203C';'YLR204W';'YLR238W';'YLR255C';'YLR257W';'YLR295C';'YLR304C';'YLR308W';'YLR312W-A';'YLR369W';'YLR377C';'YLR382C';'YLR396C';'YML022W';'YML061C';'YML088W';'YML110C';'YML129C';'YMR021C';'YMR024W';'YMR064W';'YMR066W';'YMR071C';'YMR072W';'YMR084W';'YMR089C';'YMR097C';'YMR098C';'YMR135W-A';'YMR150C';'YMR151W';'YMR158W';'YMR184W';'YMR193W';'YMR228W';'YMR231W';'YMR256C';'YMR257C';'YMR267W';'YMR282C';'YMR286W';'YMR287C';'YMR300C';'YNL001W';'YNL003C';'YNL005C';'YNL073W';'YNL081C';'YNL084C';'YNL160W';'YNL170W';'YNL177C';'YNL184C';'YNL213C';'YNL220W';'YNL225C';'YNL252C';'YNL284C';'YNL315C';'YNL316C';'YNR036C';'YNR037C';'YNR041C';'YNR050C';'YOL023W';'YOL058W';'YOL095C';'YOL096C';'YOL100W';'YOL126C';'YOL143C';'YOL148C';'YOR033C';'YOR065W';'YOR125C';'YOR128C';'YOR150W';'YOR158W';'YOR199W';'YOR200W';'YOR201C';'YOR205C';'YOR211C';'YOR241W';'YOR304C-A';'YOR305W';'YOR318C';'YOR381W';'YPL005W';'YPL006W';'YPL013C';'YPL029W';'YPL031C';'YPL040C';'YPL059W';'YPL078C';'YPL097W';'YPL104W';'YPL118W';'YPL132W';'YPL148C';'YPL172C';'YPL173W';'YPL215W';'YPL262W';'YPL271W';'YPR047W';'YPR060C';'YPR067W';'YPR072W';'YPR099C';'YPR100W';'YPR116W';'YPR124W';'YPR166C';'YPR191W';};
    genes = upper(genes);
end

function genes = alwaysEssential
    % list of 47 genes found to be essential in all tests by Kuepfer et al
    % - most are auxotrophs, and five or dubious ORFs
    genes = {'YAR015W';'YBL033C';'YBR115C';'YBR209W';'YCR066W';'YDL072C';'YDL075W';'YDR007W';'YDR008C';'YDR127W';'YDR158W';'YDR234W';'YDR300C';'YDR354W';'YDR477W';'YEL024W';'YER052C';'YER086W';'YER090W';'YFL032W';'YGL026C';'YGL148W';'YGL154C';'YGL234W';'YGR061C';'YGR168C';'YGR204W';'YHL005C';'YHR025W';'YIL094C';'YIR034C';'YJL088W';'YJR122W';'YJR139C';'YKL211C';'YKR058W';'YLR304C';'YML022W';'YMR300C';'YNL220W';'YNL316C';'YOL058W';'YOL143C';'YOR128C';'YPL006W';'YPR060C';'YPR067W';};
end