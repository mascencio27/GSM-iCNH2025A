% Script to determine aerobic growth rate of C. necator on glucose
% Taken from "What is flux balance  analysis? - Supplementary tutorial" 
% by J. D. Orth, I. Thiele, & B. 0. Palsson, Nature Biotechnology,
% Volume 28, Number 3, March 201+0

clear;

%solverOK = changeCobraSolver('glpk','LP');
solverOK = changeCobraSolver('glpk','all');

% loading iCNH2025A model

load('iCNH2025A.mat');

%saving in excel format

exchangeRxns = model.rxns(findExcRxns(model));

model = changeRxnBounds(model,exchangeRxns,0,'l'); 

% (millimoles per gram dry cell weight per hour, the default flux 
%(mmol gDW-1 hr-1)
% (millimoles per gram dry cell weight per hour, the default flux
% units used in the COBRA Toolbox), enter:

% Sugars

model_1 = changeRxnBounds(model,'R_EX_glc__D_e',-0.20,'b');
%model_1 = changeRxnBounds(model_1,'R_EX_fru_e',0,'b');

% Cofactors
model_1 = changeRxnBounds(model_1,'R_EX_pi_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_nh4_e',-0.14,'l');
model_1 = changeRxnBounds(model_1,'R_EX_so4_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_o2_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_co2_e',5.21,'u');


%Trace elements
model_1 = changeRxnBounds(model_1,'R_EX_mn2_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_cu2_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_zn2_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_fe3_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_k_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_mg2_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_ca2_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_cl_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_mobd_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_na1_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_fe2_e',-1,'l');
model_1 = changeRxnBounds(model_1,'R_EX_cobalt2_e',-1,'l');

% By setting the lower bound of the oxygen uptake reaction to such a 
% large number, it is practically unbounded. 

% Set optimization objective to Growth

model_1 = changeObjective(model_1, 'Growth', 1);
model_1.c(findRxnIDs(model_1,'R_SK_phb_c')) = 0.001;


%Perform FBA with Growth as the objective, 


FBAsolution = optimizeCbModel(model_1,'max',0,1); % FBA 

%FBAsolution = optimizeCbModel(model_1,'max',0,0); %without loops.

%Inspection of the flux distribution vector FBAsolution.x shows that 
%there is high flux

%printFluxVector(model_1, FBAsolution.x, true, true)
printFluxVector(model_1, FBAsolution.x)

Growth = FBAsolution.x(findRxnIDs(model_1,'Growth'))
PHA = FBAsolution.x(findRxnIDs(model_1,'R_SK_phb_c'))


%codigo para minización de fluxes

x0=randn(2737,1);
f=@sumsqr
Aeq=model_1.S;
beq=model_1.b;
lb=model_1.lb;
ub=model_1.ub;
options=optimset('Algorithm','sqp');
[x,fval]=fmincon(@objfun,x0,[],[],Aeq,beq,lb,ub,{},options);

%%FVA
[minFlux, maxFlux] = fluxVariability(model_1, 100,'max');