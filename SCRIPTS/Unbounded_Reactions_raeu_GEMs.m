%% Identifying Unbounded Reactions in Cupriavidus necator Genome Scale Models
% 15 / 02 / 2023
% This methodology to identify Thermodynamically infeasible cycles
% was adapted from the mathematical description of Schellenberger paper :
% Schellenberger, J., Lewis, N. E., & Palsson, B. �. (2011).
% Elimination of thermodynamically infeasible loops in steady-state metabolic models.
% Biophysical journal, 100(3), 544-553.
% Victor A. Lopez-Agudelo - University of Antioquia - Colombia

clc; clear

addpath ('models')

METABOLIC_MODEL = {'raeu2023'}; % model's names

NURs = zeros(length(METABOLIC_MODEL), 2);

solverOK = changeCobraSolver('gurobi7','LP');

%% Counting and Identifying Unbounded Reactions

for k = 1: length(METABOLIC_MODEL)

model = load_model(METABOLIC_MODEL{k});

[minFVA,maxFVA] = fluxVariability(model,100); % FVA

count_minFVA = length(minFVA(minFVA < -900));
count_maxFVA = length(maxFVA(maxFVA > 900));

count_minmax = 0;

%% Unbounded Rxns
URs_max = model.rxns(maxFVA > 900);
URs_min = model.rxns(minFVA < -900);
URs = union(URs_max,URs_min);
%% positions of Unbounded Rxns
positions_vector = (1:1:length(model.rxns))';
id_URs_max = positions_vector(maxFVA > 900);
id_URs_min = positions_vector(minFVA < -900);
id_URs = union(id_URs_max, id_URs_min); % positions of URs
%% Stoichiometric Matrix of Unbounded Reactions with metabolites and reactions
MatrixS = full(model.S);
MatrixS = MatrixS(:,id_URs); % Stoichiometric Matrix of Unbounded Reactions
%% positions of Metabolites
mets_positions_vector = zeros(length(model.mets),1);
for q = 1:length(model.mets)
    find_met_vector = find(model.S(q,:)); %
    empty = isempty(find_met_vector); % is empty "find_met_vector"?.
    if (empty == 0)
        mets_positions_vector(q) = q;
    end
end
mets_positions_vector(mets_positions_vector == 0) = [];
MatrixS = MatrixS(mets_positions_vector,:);
NullSpace = null(MatrixS,'r'); % vectors that form the null space of the S matrix.
UBrxnList = model.rxns(id_URs); % Unbounded Reaction list sorted
[mrow,ncol] = size(NullSpace);
% Finding the LB and UB for the Unbounded Reactions
model.lb = model.lb(id_URs);
model.ub = model.ub(id_URs);

%% writing text in file with possible TICs
folder = cd;  %
filename = fullfile(folder, sprintf('TICs_%s.txt', METABOLIC_MODEL{k})); % create file with Thermodynamical infeasible cycles
fileID = fopen(filename, 'a');
if fileID == -1, error('Cannot open file: %s', filename); end

TIC_count = 0; % TIC counter
 for z = 1:ncol

     fprintf(fileID, '###########\n'); % No. of TIC
     rxncount = 0; % counter of added reaction in List
     for v = 1:length(UBrxnList)
         if (model.lb(v) < 0) & (model.ub(v) > 0) & (NullSpace(v,z) ~= 0)
           TIC = UBrxnList(v); % Extreme pathways formed by unbounded reactions (TICs)
           fprintf(fileID, '%s \n', TIC{1:end});
           rxncount = rxncount + 1;
         end
         if (model.lb(v) < 0) && (model.ub(v) == 0) & (NullSpace(v,z) < 0)
           TIC = UBrxnList(v); % Extreme pathways formed by unbounded reactions (TICs)
           fprintf(fileID, '%s \n', TIC{1:end});
           rxncount = rxncount + 1;
         end
         if (model.lb == 0) & (model.ub(v) > 0) & (NullSpace(v,z) > 0)
           TIC = UBrxnList(v); % Extreme pathways formed by unbounded reactions (TICs)
           fprintf(fileID, '%s \n', TIC{1:end});
           rxncount = rxncount + 1;
         end
     end
     if (rxncount > 1)
          TIC_count = TIC_count + 1;
          fprintf(fileID, 'TIC_#%d.\n', TIC_count); % No. of TIC
     elseif (rxncount == 0)
          fprintf(fileID, '#####THIS IS NOT A TIC#####\n');
     else
          fprintf(fileID, '#####THIS IS NOT A TIC#####\n'); % That's not a TIC

     end
 end
fprintf(fileID, '###########\n');
fclose(fileID); % closing file

for i = 1:length(model.rxns)

    if (minFVA(i) < -900) && (maxFVA(i) > 900)

      count_minmax = count_minmax + 1;

    end

end

count_UR = (count_minFVA - count_minmax) + (count_maxFVA - count_minmax) + count_minmax

NURs(k,1) = count_UR;
NURs(k,2) = (count_UR/length(model.rxns)*100);

end

%% Summary Results

Number_of_Unbounded_Reactions = table (NURs(:,1),NURs(:,2),'RowNames',METABOLIC_MODEL,'VariableNames',{'Unbounded_Rxns','Percentage_URs'})
file_excel = 'Results_UBs.xlsx';
writetable(Number_of_Unbounded_Reactions,file_excel,'Sheet',1,'WriteRowNames',true)
