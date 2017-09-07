% This script is created to aid in creating the table for a summary. 
% It generates the formatted result of 1- to N-back serial dependency analysis of
% the dataset, using the format of the summary table, the result T can be copied into the
% summary directly.

% please enter the correct name of the script you are using on the 9th and 12th line.
% (bkitSequenceCW,bkitSequenceFTV,bkitSequenceVFA)

clear T
N = 5;
T = table(bkitSequenceFunc(1));
for i=2:N
    % enter function name here
    x = table(bkitSequenceFunc(i));
    x.Properties.VariableNames{'Var1'} = ['a', num2str(i)];
    T = [T ,x];
end
