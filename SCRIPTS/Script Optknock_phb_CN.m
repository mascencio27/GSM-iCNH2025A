%% This code is an adaptation of the Cobra Toolbox OptKnock Tutorial developed by Sebastian N. Martinez
% source: https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialOptKnock.html
% Adapted by: Martha L. Ascencio-Galvan, University of Valle - Cali - Colombia

%%
changeCobraSolver('gurobi','all');

%%
% loading iCNH2025A model
load('iCNH2025A.mat');

biomass = 'Growth';
%% 
% Define the maximum number of solutions to find (i.e., maximum number of 
% remvable reactions that lead to the overproduction of the metabolite of interest)

threshold = 5;
%% 
% Define the set of reactions that will be used to search knockouts. Note 
% : only reactions in this set will be deleted

selectedRxnList = {'R_GLCabcpp'; 'R_GLCptspp';'R_HEX1';'R_FBP';'R_FBA';'R_GAPD';'R_TPI';'R_PGK';'R_PGM';'R_ENO';
                   'R_PYK';'R_PDH';'R_PGI';'R_G6PDH2r';'R_PGL';'R_EDD';'R_EDA';'R_ACACT1r';'R_AACOAR_syn';'R_CELLS';'R_CAT';
                   'R_ACTDa';'R_ACTD_1';'R_ACONTb';'R_ACONTa';'R_ACOADH2';'R_ACLSb_1';'R_ACLSa_1';'R_ACLS';'R_ACCOAL';'R_2DHGLCK';
                   'R_CYO1_KT';'R_CYO1b';'R_CYTBD';'R_CYTBDpp';'R_CYTCAA3pp';'R_CYTCBB3pp';'R_DHAD1';'R_F6PP';'R_FDH';'R_GALUi';
                   'R_GLUDy';'R_H2ASE_syn';'R_HADPCOADH';'R_HCO3E';'R_LEUt3pp';'R_MI1PP';'R_PIt2r';'R_PItex';'R_PItpp';'R_PPS';'R_PY5CCR';
                   'R_PYDXNtr';'R_THD2pp';'R_THRA2';'R_VALTA';'R_VPAMTr';'R_PPC';'R_PDH';'R_CS';'R_ACONTa';'R_ICDHyr';'R_ICL';'R_AKGDH';
                   'R_SUCDi';'R_FUM';'R_MDH';'R_MALS';'R_ASPTA6'};

%%
exchangeRxns = model.rxns(findExcRxns(model));

model = changeRxnBounds(model,exchangeRxns,0,'l'); 

% (millimoles per gram dry cell weight per hour, the default flux 
%(mmol gDW-1 hr-1)
% (millimoles per gram dry cell weight per hour, the default flux
% units used in the COBRA Toolbox), enter:

% Sugars

model = changeRxnBounds(model,'R_EX_glc__D_e',-0.12,'b');

% Cofactors
model  = changeRxnBounds(model,'R_EX_pi_e',-1,'l');
model  = changeRxnBounds(model,'R_EX_nh4_e',-0.089,'l');
model  = changeRxnBounds(model,'R_EX_so4_e',-1,'l');
model  = changeRxnBounds(model,'R_EX_o2_e',-1,'l');
model  = changeRxnBounds(model,'R_EX_co2_e',5.21,'u');


% Trace elements
model = changeRxnBounds(model,'R_EX_mn2_e',-1,'l');
model = changeRxnBounds(model,'R_EX_cu2_e',-1,'l');
model = changeRxnBounds(model,'R_EX_zn2_e',-1,'l');
model = changeRxnBounds(model,'R_EX_fe3_e',-1,'l');
model = changeRxnBounds(model,'R_EX_k_e',-1,'l');
model = changeRxnBounds(model,'R_EX_mg2_e',-1,'l');
model = changeRxnBounds(model,'R_EX_ca2_e',-1,'l');
model = changeRxnBounds(model,'R_EX_cl_e',-1,'l');
model = changeRxnBounds(model,'R_EX_mobd_e',-1,'l');
model = changeRxnBounds(model,'R_EX_na1_e',-1,'l');
model = changeRxnBounds(model,'R_EX_fe2_e',-1,'l');
model = changeRxnBounds(model,'R_EX_cobalt2_e',-1,'l');

model = changeRxnBounds(model,'R_SK_phb_c',0.0,'l'); % set lower bound to 0.0
model = changeRxnBounds(model,'R_SK_phb_c',1000,'u'); % set upper bound to 1000

% By setting the lower bound of the oxygen uptake reaction to such a 
% large number, it is practically unbounded. 

% Set optimization objective to Growth

model = changeObjective(model, 'Growth', 1);
model.c(findRxnIDs(model,'R_SK_phb_c')) = 0.001;

%Perform FBA with Growth as the objective, 

FBAsolution = optimizeCbModel(model,'max',0,1); % FBA

Growth = FBAsolution.x(findRxnIDs(model,'Growth'))
PHA = FBAsolution.x(findRxnIDs(model,'R_SK_phb_c'))
%% 
% Then, calculates the production of metabolites before running optKnock.

% determine PHB production and growth rate
fbaWT = optimizeCbModel(model,'max',0,1);
phbFluxWT = fbaWT.x(findRxnIDs(model,'R_SK_phb_c'));
growthRateWT = fbaWT.x(findRxnIDs(model,'Growth'));
fprintf('The production of phb before optimization is %.4f \n', phbFluxWT);
fprintf('The growth rate before optimization is %.4f \n', growthRateWT);

%%
fprintf('\n...EXAMPLE 1: Finding optKnock sets of size 3 or less...\n\n')
% Set optKnock options
% The exchange of PHB will be the objective of the outer problem
options = struct('targetRxn', 'R_SK_phb_c', 'numDel', 3);
% We will impose that biomass be at least 30% of the biomass of wild-type
constrOpt = struct('rxnList', {{biomass}},'values', 0.3*fbaWT.x(findRxnIDs(model,'Growth')), 'sense', 'G');
% We will try to find 10 optKnock sets of a maximun length of 2
previousSolutions = cell(10, 1);
contPreviousSolutions = 1;
nIter = 1;
while nIter < threshold
    fprintf('...Performing optKnock analysis...')
    if isempty(previousSolutions{1})
        optKnockSol = OptKnock(model, selectedRxnList, options, constrOpt);
    else
        optKnockSol = OptKnock(model, selectedRxnList, options, constrOpt, previousSolutions);
    end
    
    % determine PHB production and growth rate after optimization
    phbFluxM1 = optKnockSol.fluxes(strcmp(model.rxns, 'R_SK_phb_c'));
    growthRateM1 = optKnockSol.fluxes(strcmp(model.rxns,biomass));
 
    setM1 = optKnockSol.rxnList;
    
    if ~isempty(setM1)
        previousSolutions{contPreviousSolutions} = setM1;
        contPreviousSolutions = contPreviousSolutions + 1;
        %printing results
        fprintf('optKnock found a optKnock set of large %d composed by ', length(setM1));
        for j = 1:length(setM1)
            if j == 1
                fprintf('%s', setM1{j});
            elseif j == length(setM1)
                fprintf(' and %s', setM1{j});
            else
                fprintf(', %s', setM1{j});
            end
        end
        fprintf('\n');
        fprintf('The production of phb after optimization is %.4f \n', phbFluxM1);
        fprintf('The growth rate after optimization is %.4f \n', growthRateM1);
                fprintf('...Performing coupling analysis...\n');
        [type, maxGrowth, maxProd, minProd] = analyzeOptKnock(model, setM1, 'R_SK_phb_c');
        fprintf('The solution is of type: %s\n', type);
        fprintf('The maximun growth rate given the optKnock set is %.4f\n', maxGrowth);
        fprintf(['The maximun and minimun production of phb given the optKnock set is ' ...
                 '%.4f and %.4f, respectively \n\n'], minProd, maxProd);
        singleProductionEnvelope(model, setM1, 'R_SK_phb_c', biomass, 'savePlot', 1, 'showPlot', 1, ...
                                 'fileName', ['PHB_SK2_' num2str(nIter)], 'outputFolder', 'OptKnockResults');
    else
        if nIter == 1
            fprintf('optKnock was not able to found an optKnock set\n');
        else
            fprintf('optKnock was not able to found additional optKnock sets\n');
        end
        break;
    end
    nIter = nIter + 1;
end
%cd(currectDirectory);