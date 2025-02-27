%clear;
 
%solverOK = changeCobraSolver('glpk','LP');
%addpath('MC3')
%[SCM, DEM, ZFR, UR, CR, RCR] = mc_checkmodel ('SBML', 2, 'raeu_2909v6.xml');

clear;
 
solverOK = changeCobraSolver('glpk','LP');
addpath('MC3')
[SCM, DEM, ZFR, UR, CR, RCR] = mc_checkmodel ('SBML', 2, 'raeu_2909v6.xml');
