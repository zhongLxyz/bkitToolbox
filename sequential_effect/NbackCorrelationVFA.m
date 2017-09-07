%function meanResults = NbackCorrelationVFA(allData,N,varargin)

% this function reads in the dataset from a viewing perspective condition and finds the
% correlation between each trial and its previous 1~N trials (see below for
% more information)
% This is one measure of the serial dependency effect. 

% outputs: results in one table for each session. Each column is the
% correlation r value between each trial
% and its N-back trials. e.g. column 1 is the r value of each trial and its
% 1-back (previous) trial, column 2 is the r of each trial and its 2-back
% trial etc. The first row is the r-values of the correlation analysis, and
% the 2nd row is the p-values.
% in the end it outputs a table which has the average values across all
% sessions.

% modify this section to specify the dataset:
allData = bkitRead('VFA');
variable = 'ccw';  % what variable you want to look at ('ccw','ftv','vfa')
N = 6;  % N-back response;

% results is a 3D array, in each layer it stores the 1 to N-back
% corelation results of one participant.
results = zeros(2, N, length(allData));

for i=1:length(allData)
    previousTrial = zeros(2,N);
    singleSubj = allData{i};
    fprintf('%s:',singleSubj.Session{1});

    response = singleSubj.Response;           % 0: clockwise,  1: ccw
    angularVel = singleSubj.AngularVelocity;     % spinning direction
    azimuth = singleSubj.CameraAzimuth;
    elevation = singleSubj.CameraElevation;             % camera elevation

    nWithNans = [];
    
    if strcmp('ftv',variable)
        d = findFTV(angularVel,azimuth,response,elevation,'vfa');  % d is the facing direction perceived. 1: FTV, 0: FA
        for j = 1:N
            [RHO,PVAL] = corr(d(1:end-j),d(j+1:end));
            previousTrial(:,j) = [RHO; PVAL];
        end
    elseif strcmp('ccw',variable)
        d = findCCW(angularVel,azimuth,response,elevation,'vfa');
        for j = 1:N
            [RHO,PVAL] = corr(d(1:end-j),d(j+1:end));
            previousTrial(:,j) = [RHO; PVAL];
        end
    elseif strcmp('vfa',variable)
        d = findVFA(angularVel,azimuth,response,elevation,'vfa');
        for j = 1:N
            [RHO,PVAL] = corr(d(1:end-j),d(j+1:end));
            previousTrial(:,j) = [RHO; PVAL];
            if isequal(ones(size(d(j+1:end))),d(j+1:end))
                nWithNans = [nWithNans j];
            end
        end
    end
    results(:,:,i) = previousTrial;
    indivTable = array2table(previousTrial,'RowNames',{'r-value','p-value'});
    disp(indivTable);
    if not(isempty(nWithNans))
        disp('NOTE: for the following N, all response are "above" in the second group, therefore the calculation is not possible.');
        disp(nWithNans);
    end
    disp('-------------------------------------')
end

meanResults = nanmean(results,3);
meanTable = array2table(meanResults,'RowNames',{'r-value','p-value'});
disp('average value across sessions:')
disp(meanTable)