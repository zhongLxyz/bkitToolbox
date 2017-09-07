function record = bkitGenerate(parameter, design, allBalancedVar, ntrial)
% This function generates simulation of a experiment session using parameter values
% fitted to the GLM from bkitModel, along with the designs of that experiment from 
% bkitSummary. It will calculate the simulated response of a 'participant' from 
% the given values.
%
% input: parameter - a table of model parameters from a single session(participant), output
% from bkitModel (must only select the data from one session)
% design - a table with information about one particular design of an experiment, output from bkitSummary
% allBalancedVar - a cell array, its entries are names of balanced variables (in strings)
%     in that design of the experiment.
% ntrials - the number of trials you want to generate, if there are balanced variables in the design, 
%     ntrials would better be multiples of the number of trials in a full block
%
% **If the variables in the parameter table has not been specified in the current function before,
%     an error will occur during the execution of the function. The user need to edit the function
%     and specify what kind of processing has to be done for that new variable.
% **If a variable is included in the parameter, but not the design, then the level of the variable will
%     show up as 99999 in the final output. Please look at the input parameter and adjust accordingly.
% **assumption: this function assumes that, if there is any balanced variable, AND if angular velocity
%     does not just contain one level, then angular velocity must be among one of the balanced vars.

selectionPool = []; % a new variable for putting in the given levels in order, the 1 is
slctPoolVarNames = {};
% for the intercept
% the first column of parameter is always intercept, so ignore this column for now
for jj = 2:width(parameter)
    kk = 1;
    while (kk <= width(design)) && (~strcmpi(parameter.Properties.VariableNames{jj}, design.Properties.VariableNames{kk}))
        kk = kk + 1;
    end
    if kk <= width(design) 
        selectionPool = [selectionPool design{1,kk}];
        slctPoolVarNames = [slctPoolVarNames {parameter.Properties.VariableNames{jj}}];
    else
        selectionPool = [selectionPool {99999}]; % assign a indicator value to remind processing later
        slctPoolVarNames = [slctPoolVarNames {parameter.Properties.VariableNames{jj}}];
    end
end
selectionPool = array2table(selectionPool,'VariableNames',slctPoolVarNames);

balancedVarName = {};
if ~isempty(allBalancedVar)
    % for all variable in the balanced names, collect them and generate all comb of their levels
    balancedVarLevel = {};
    designVarName = design.Properties.VariableNames;
    for ii = 1:length(allBalancedVar)
        for jj = 1:width(design)
            if strcmpi(allBalancedVar{ii},designVarName{jj})
                balancedVarLevel = [balancedVarLevel design{:,jj}];
                balancedVarName = [balancedVarName {allBalancedVar{ii}}];
            end
        end
    end
    % put all variables with the same name in allcomb and generate a list of trials
    allComb = allcomb(balancedVarLevel{:});
    allComb = array2table(allComb,'VariableNames',balancedVarName);
    angVelComb = [];
    % treat the balanced variables to prepare them for later calculation according to what they are
    % separate the angularVel variable column from the rest, cuz it will not enter X
    if contains('AngularVelocity',balancedVarName)
        idx = find(strcmpi('AngularVelocity',balancedVarName));
        angVelComb = sign(allComb{:,idx});
        allComb(:,idx) = [];
        balancedVarName(:,idx) = [];
    elseif length([design.AngularVelocity{:}])==1
        angVelComb = design.AngularVelocity{:};
    end
    rawBalancedVar = allComb{:,:};
    if ~isempty(angVelComb)
        if contains('CameraAzimuth',balancedVarName)
            h = (allComb.CameraAzimuth + 90) * pi/180;
            fColumn = 1.232*sin(h) + 0.33*sin(3*h) + 0.14*sin(5*h) + 0.052*sin(7*h) + 0.01*sin(9*h);
            allComb.CameraAzimuth = angVelComb .* fColumn;
        end
        if sum(strcmpi('CameraElevation',balancedVarName))==1
            allComb.CameraElevation = angVelComb .* allComb.CameraElevation;
        end
    end
end
% find the non-balanced vars
nonBalancedVarName = parameter.Properties.VariableNames(~ismember(parameter.Properties.VariableNames,allBalancedVar));
nonBalancedVarName = nonBalancedVarName(2:end);

angularVel = 0;
for ii = 1:ntrial
    X = [1]; % X is the matrix for keeping all setting in a trial, 1 is for dealing with the intercept
    
    angularVel1Back = angularVel;
    if ~isempty(allBalancedVar)
        % randomly select a row from allcomb, tap it onto X
        % if there is balanced variable
        index = randsample(size(allComb,1),1);
        randomRow = allComb{index,:};
        X = [X randomRow];
    end
    if contains('AngularVelocity',allBalancedVar)
        angularVel = angVelComb(index,:);
    else
        angularVel = sign(design.AngularVelocity{:}(randi(length(design.AngularVelocity{:}))));
    end
    % dealing with history factors
    if ii>1
        f1 = f;
        elev1Back = elev;
        resp1Back = resp*2-1;
        azimu1Back = azimuth;
    else
        f1 = 0;
        elev1Back = 0;
        resp1Back = 0;
        azimu1Back = 0;
    end
    
    % dealing with the non-balanced variables, add them to X
    if ~isempty(allBalancedVar)
        rawVars = rawBalancedVar(index,:);
    else
        rawVars = [];
    end
    for jj = 1:length(nonBalancedVarName)
        slctVarName = nonBalancedVarName{jj};
        if strcmpi(slctVarName,'CameraAzimuth')
            % create variable for camera azimuth
            Azlevel = selectionPool.CameraAzimuth{1};
            azimuth = Azlevel(randi(length(Azlevel)));
            h = (azimuth + 90) * pi/180;
            f = 1.232*sin(h) + 0.33*sin(3*h) + 0.14*sin(5*h) + 0.052*sin(7*h) + 0.01*sin(9*h);
            CameraAzimuth = angularVel*f;
            X = [X CameraAzimuth];
            rawVars = [rawVars azimuth];
        elseif strcmpi(slctVarName,'CameraElevation')
            Elevlevel = selectionPool.CameraElevation{1};
            elev = Elevlevel(randi(length(Elevlevel)));
            X = [X angularVel*elev];
            rawVars = [rawVars elev];
        elseif strcmpi(slctVarName,'CCW1Back')
            X = [X resp1Back];
            rawVars = [rawVars resp1Back];
        elseif strcmpi(slctVarName,'FTV1Back')
            X = [X angularVel1Back*resp1Back*f1];
            rawVars = [rawVars resp1Back*f1];
        elseif strcmpi(slctVarName,'VFA1Back')
            if ii>1
                vfa1Back = findVFA(angularVel1Back,azimu1Back,round(resp1Back),elev1Back);
                if isempty(vfa1Back)
                    x = findVFA(angularVel1Back,azimu1Back,round(resp1Back),elev1Back);
                end
            else
                vfa1Back = 0;
            end
            X = [X vfa1Back];
            rawVars = [rawVars vfa1Back];
        else
            % when the variable hasn't been specified in the script yet
            fprintf('ERROR: The calculation for Variable "%s" in model parameters has not been specified, please code in its calculation in the script',slctVarName);
            return
        end
    end
    if ~isempty(allBalancedVar)
        if contains('CameraAzimuth',balancedVarName)
            azimuth = randomRow(strcmpi('CameraAzimuth',balancedVarName));
            f = fColumn(index);
        end
        if contains('CameraElevation',balancedVarName)
            elev = randomRow(strcmpi('CameraElevation',balancedVarName));
        end
    end
    if width(parameter)~=length(X)
        disp('x');
    end
    g = parameter{:,:} * X';
    resp = (1./(1+exp(-g))); % resp is the predicted cw/ccw response
    record(ii,:) = [angularVel*45 rawVars resp];
end

record = array2table(record,'VariableNames', ['AngularVelocity' balancedVarName nonBalancedVarName 'Response']);
