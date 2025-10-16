
%% ANOVA Testing for total phenol

clear
clc 

% Defining data
Data = 'Minerals';
% define dimension of HSD test#
D = [1];
switch Data
    case 'Compounds'
        T = readtable("Seaweed_ALL.xlsx", 'Sheet', 'Matlab_Phenol');
        % Perform two-way ANOVA and label variables
        [~,~,stats] = anovan(T.Concentration, {T.Method, T.Material, T.Solvent}, 'Model', 'interaction', 'varnames', {'Method', 'Material', 'Solvent'});
    
    case'Minerals'
        % Create table from minerals spreadsheet
        T = readtable("Seaweed_Minerals_Matlab.xlsx", 'Sheet', 'Matlab_Minerals');
        % Perform two-way ANOVA and label variables
        [~,~,stats] = anovan(T.Concentration, {T.Method, T.Material, T.Solvent, T.Element}, 'Model', 'interaction', 'varnames', {'Method', 'Material', 'Solvent', 'Element'});

        % Display the results
        significant_pairs = pairs(p_values < 0.05, :);
        significant_p_values = p_values(p_values < 0.05);
        
        % Display the significant results
        disp('Significant Pairwise Comparisons (p < 0.05):');
        disp(table(significant_pairs(:,1), significant_pairs(:,2), significant_p_values, 'VariableNames', {'SType1', 'SType2', 'pValue'}));
end


% Create a table
% T = readtable("Seaweed_ALL.xlsx", 'Sheet', 'Matlab_DPPH');
% Or for Minerals:
% T = readtable("Seaweed_Minerals_Matlab.xlsx", 'Sheet', 'Matlab_Minerals');

% Create a new variable that represents the combination of 'Material' and 'Solvent'
% This group will become the 6th column in T
% T.Material_Solvent = strcat(T.Material, "_", T.Solvent);


% % Perform Tukey's HSD test
% % Dimension: Which "groups" of "factors" are significantly different
% [results,~,~,gnames] = multcompare(stats, 'Dimension', D, 'CType', 'hsd');
% 
% tbl = array2table(results,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% tbl.("Group A")=gnames(tbl.("Group A"));
% tbl.("Group B")=gnames(tbl.("Group B"))



