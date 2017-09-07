function [design, allBalancedVar] = bkitSummary(allData)
% This function reads in the output of bkitRead (the output from all participants) and outputs a summary
% for the dataset. The summary includes the design of the experiment. The design is 
% what variables are in the dataset and what levels they have. A dataset passed to it
% can have multiple designs in it. It will also find what variables are balanced together
% in the dataset
%
% terms I used here:
% variable: things like camera azimuth, elevation, gender etc. They are variables in the
% experiments, but matlab table column titles are also called variables. And here
% the table variables are indeed the experiment variables.
% level: e.g. in camera elevation in Zhong's, there are 2 levels, -20 and 20.
% design: if two participants have two different designs for the same variable, that means the 
%   two participants experienced different levels in the same variable.
%
% input: a cell array which is the dataset output from bkitRead.
% output: design - a table. Each row of the table contains the information of one design,
% the table variables are the variables in the experiment. And each entry has the corresponding levels
% in each variable in the corresponding design.
%
%  **It is assumed here that all designs in the same dataset passed to it at one time
%  has the same set of variables that are balanced together.
%  **It also assumes that all designs in the dataset that is passed to it at one time
%  all contains the same variables in all sessions (although the level combinations in each
%  variable might be different.

% things to add: design should include all variables, not just the ones that have more uniques

% for all participant, find their designs (levels of variables) and store
% them in an array.
singleSubj = allData{1};
NofColumn = width(singleSubj);
% summary is a (length(allData)+1)*NofColumn cell array, the first row stores the names of variables,
% the second row stores the corresponding combination of variables in the first row. 
summary = cell(size(allData,2)+1,NofColumn);
for ii = 1:length(allData)
    singleSubj = allData{ii};
    summary(1,:) = singleSubj.Properties.VariableNames;
    for jj = 1:NofColumn
        summary{ii+1,jj} = unique(singleSubj{:,jj});
    end
end

% if the levels in a variable is different for every session, then it must not be an important variable
% (e.g. session ID), so I will exclude that variable from the counting.
ColumnNumtoExclude = [];
for jj = 1:size(summary,2)
    UniqueLevel = uniquecell(summary(2:end,jj));
    % record the # of unique levels of each variable, if it is almost equal to the number of sessions
    % do not include this variable.
    if max(size(UniqueLevel)) >= (size(allData,2)-3)
        ColumnNumtoExclude = [ColumnNumtoExclude jj]; % keep track of the column number to exclude it later
    end
end

allColumn = 1:size(summary,2);
selectedColumn = allColumn(~ismember(allColumn,ColumnNumtoExclude));
keepVariable = summary(:,selectedColumn);

design = [];
for jj = 1:size(keepVariable,2)
    diffDesign = uniquecell(keepVariable(2:end,jj));
    % pad the columns missing a few cells with Nans
    % if the row number of the whole design array is smaller than the number of cells in the diffDesign 
    if (size(design,1) < size(diffDesign,1))
        % add rows of NaN to the end of the design array
        numofRowtoAdd = size(diffDesign,1) - size(design,1);
        design = [design; NaN([numofRowtoAdd,size(design,2)])];
    % pad diffDesign with zero, if it is the shorter one
    elseif (size(design,1) > size(diffDesign,1))
        numofRowtoAdd = size(design,1) - size(diffDesign,1);
        diffDesign = [diffDesign; NaN([numofRowtoAdd,1])];
    end
    design = [design diffDesign];
end

% add titles to design matrix
design = array2table(design,'VariableNames',keepVariable(1,:));

%% count the occurence of each level in variables to determine
% balanced/random design (assume all designs in the same dataset always have the same 
% set of variables as balanced)
singleSubjData = allData{1};
% delete all columns with no variance in them
cleanedData1 = singleSubjData(:,selectedColumn);
cleanedData = [];
for ii = 1:length(cleanedData1.Properties.VariableNames)
    if var(cleanedData1{:,ii})~=0
        cleanedData = [cleanedData cleanedData1(:,ii)];
    end
end

allVariable = cleanedData.Properties.VariableNames;
nofAllVariable = length(allVariable);
allBalancedVar = {};

 % this should be the data of one participant with a particular design
jj = 1;
while isempty(allBalancedVar) && (jj<=nofAllVariable)
    firstVarData = cleanedData.(char(allVariable(jj)));
    firstVarLevel = unique(firstVarData);
    % for each level in the current selected variable
    for ee = (jj+1):nofAllVariable
        varToCompareTo = allVariable(ee);
        varToCompareToLevel = unique(cleanedData.(char(varToCompareTo)));
        % every level in firstVar should have all combinations of the varToCompareTo's levels (or multiples of it)
        % for the two vars to be considered balanced together.
        trueOrFalse = true;
        kk = 1;
        while trueOrFalse && (kk <= length(firstVarLevel))
            % select all rows that has this level in the participant's dataset
            rowsWithSlctLv = cleanedData(firstVarData==firstVarLevel(kk),:);
            % for every variable selected, count if it has a full set of levels or multiples of it,
            % true means it is
            trueOrFalse = allCombForOneLevel(varToCompareToLevel,rowsWithSlctLv.(char(varToCompareTo)));
            kk = kk + 1;
        end
        % add it to a list of balanced variables if it is true
        if trueOrFalse
            allBalancedVar = [allBalancedVar allVariable(ee)];
        end
    end
    if ~isempty(allBalancedVar)
        allBalancedVar = [allVariable(jj) allBalancedVar];
    end
    jj = jj + 1;
end
end
%% helper functions

function result = allCombForOneLevel(level,data)
    allCount = zeros(1,length(level));
    for ii = 1:length(level)
        allCount(ii) = sum(data==level(ii));
    end
    if var(allCount)==0
        result = true;
    else
        result = false;
    end
end