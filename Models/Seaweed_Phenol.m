
%% ANOVA Testing for total phenol
% Assuming your data is in a table T
% with columns 'method', 'subgroup', 'concentration', and 'optical_depth'
clear
clc 

% T = readtable("Seaweed_ALL.xlsx", 'Sheet', 'Matlab_Phenol')
% head(T)
% 
% % Convert categorical variables to categorical data type
% T.method = categorical(T.method);
% T.subgroups = categorical(T.subgroups);
% 
% % Fit the model
% lm = fitlm(T, 'concentration ~ method + subgroups + method:subgroups');
% 
% % Perform ANOVA and display table
% anova(lm, 'summary')


% Define your data
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
    case 'Ahmad MC'
        T = readtable("Seaweed_Ahmed_MatlabCopy.xlsx", 'Sheet', 'Sheet3');
        % Perform one-way ANOVA
        [~,~,stats] = anovan(T.FinalMC, {T.Method, T.Temp, T.AirflowRate}, 'Model', 'interaction', 'varnames', {'Method', '[Temp','AirflowRate]'});        
        % Display the results
        disp('ANOVA Table:');
        % disp(tbl);        
        % % If you want to perform multiple comparisons 
        [results,~,~,gnames] = multcompare(stats, 'Dimension', D, 'CType', '');
         tbl = array2table(results,"VariableNames", ...
            ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
        tbl.("Group A")=gnames(tbl.("Group A"));
        tbl.("Group B")=gnames(tbl.("Group B"))
        significantResults = tbl(0 < tbl{:, 6} < 0.05, :);
    case 'Ahmad TPC' 
        T = readtable("Seaweed_Ahmed_MatlabCopy.xlsx", Sheet="Sheet3");
        T=T(1:33,:);
        % Perform one-way ANOVA
        [p, tbl, stats] = anova1(T.UVabs, T.Temp);
        % Display the results
        disp('ANOVA Table:');
        disp(tbl);
        [results,~,~,gnames] = multcompare(stats, 'Dimension', [1], 'CType', '');

        tbl = array2table(results,"VariableNames", ...
            ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
        tbl.("Group A")=gnames(tbl.("Group A"));
        tbl.("Group B")=gnames(tbl.("Group B"))

        
        % If you want to perform multiple comparisons
        multcompare(stats);
    case 'Ahmad TPC (tTest)'
        % Perform pairwise t-tests  
        T = readtable("Seaweed_Ahmed_MatlabCopy.xlsx", Sheet="Sheet4");
        [pairs, p_values] = deal([]);
        UVabs = zeros(12,3);
        UVabs(:,:) = [T.UVabs1, T.UVabs2, T.UVabs3];
        SType =  repmat(1:11, 3, 1);    
        SType = SType(:);

for i = 1:11
    for j = i+1:11
        [~, p] = ttest2(UVabs(SType == i), UVabs(SType == j));
        pairs = [pairs; i, j];
        p_values = [p_values; p];
    end
end

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


