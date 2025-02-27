function model = load_model(metabolic_model)
% Load a genome scale model for the selected organism.
% Options:
% metabolic_model: 'iRH2022A',
% Author: Victor A. Lopez-Agudelo, Junio 12, 2019.


%% MODEL PATH
  
    raeu2023_MODEL= 'model_tic162_2';

%% CHOICES

switch metabolic_model


    case 'raeu2023'

                   load(raeu2023_MODEL);
                   model = model;
                   model.lb(model.lb < 0) = -1000;
                   model.ub(model.ub > 1) = 1000;
                   exchangeRxns  = model.rxns(findExcRxns(model));
                   model = changeRxnBounds(model,exchangeRxns,-1,'l'); % substrate input
                   %obj_functions = 'Growth';
                   %model = changeObjective(model,obj_functions,1); % Set Biomass as Objective Function
                   model = changeObjective(model,'Growth',1);

                                
                   
end


end
