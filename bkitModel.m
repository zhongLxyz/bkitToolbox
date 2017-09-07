function M = bkitModel(allData, type, varargin)
% This function finds the model parameters from glm for the specified datasets,
% it basically finds what variable setting in a walker has an effect on the 
% perception of the participants, and how big the effect is.
% for each individual dataset, it outputs a set of model parameters.
%
% If the datasets are from facing directin or viewing perspective condition, it converts the responses to coding
% for cw/ccw and then do the usual analysis.
%
% input arguments: 
% allData: a cell array of dataset, output from bkitRead
% type: a char array, it specifies what condition the datasets are from (all dataset in each run of the script
% have to come from the same condition) - ('ftv', 'ccw' or 'vfa').
% varargin: char arrays separated by comma, they are the names of variables you would like to include in the modeling，
%    (can be among: 'CameraElevation','CameraAzimuth','CCW1Back','FTV1Back','VFA1Back' etc.)
%
% output: 
% M - a table containing all the parameter values fitted from the GLM, listed under corresponding
%    variale names.

for i = 1:length(allData)
    % singleSubj is a table that contains the relevant data for a single participant
    singleSubj = allData{i};
    %gender = singleSubj.Gender;
    %dynamic = (singleSubj.Speed) * 2 - 1; % convert to -1 and 1 (instead of 0 and 1)
    h = (singleSubj.CameraAzimuth + 90) * pi/180;
    % f ranges from -1 to 1. a measure related to the real facing
    % direction (its value can be close to -1, 1 or somewhere in between)
    f = 1.232*sin(h) + 0.33*sin(3*h) + 0.14*sin(5*h) + 0.052*sin(7*h) + 0.01*sin(9*h);
    elevation = singleSubj.CameraElevation;
    angularVel = singleSubj.AngularVelocity;
	
	% turn the response into coding for cw/ccw.
    y = findCCW(angularVel,singleSubj.CameraAzimuth,singleSubj.Response,elevation,type);
    
    correct = double(y==sign(angularVel));
    correct(correct==0) = -1; % -1 and 1
    correct1B = [0; correct(1:end-1)]; % 1-back correct
    
    % ***alter the elevation values to -1 and 1 for Zhong's 2nd experiment for scaling correctly.
    elevation = sign(elevation);
    
    % history factors
    response = y*2-1; % scale the response by turning it into -1 and 1 (y range from 0 to 1)
    resp1Back = [0; response(1:end-1)];   % previous response (cw/ccw) , y1
    % dependent variable in the model should range from 0 to 1, but everything else should be -1 and 1
    f1 = [0; f(1:end-1)];
    FTV1Back = correct1B.*f1; % one-back perceived facing direction history (which of the two facing direction was perceived in the last trial)
    realElev1Back = [0; elevation(1:end-1)]; % one-back actual elevation
    VFA1Back = correct1B.*realElev1Back; % one-back perveived elevation history
    
    % build the factors matrix
    X = [];
    varInX = [];
    if sum(strcmpi('CameraElevation',varargin))==1 % compare case-insensitive %
        X = [X sign(angularVel).*elevation];
        varInX = [varInX {'CameraElevation'}];
    end
    if sum(strcmpi('CameraAzimuth',varargin))==1
        X = [X sign(angularVel).*f];
        varInX = [varInX {'CameraAzimuth'}];
    end
    if sum(strcmpi('CCW1Back',varargin))==1
        X = [X resp1Back];
        varInX = [varInX {'CCW1Back'}];
    end
    if sum(strcmpi('FTV1Back',varargin))==1
        X = [X sign(angularVel).*FTV1Back];
        varInX = [varInX {'FTV1Back'}];
    end
    if sum(strcmpi('VFA1Back',varargin))==1
        X = [X sign(angularVel).*VFA1Back];
        varInX = [varInX {'VFA1Back'}];
    end
    
    % modify X and y to stay away from extreme values (change the number of
    % zeros according to how many factors there are in total)
    v = unique(elevation);
    for j = 1:length(v)
        if var(y(X(:,1)==v(j))) == 0 % if all y values corresponding to a certain elevation has the same answer.
            X = [X; [v(j), 0,0,0,0]; [v(j), 0,0,0,0]];
            y = [y; 0; 1];
        end
    end
    
    % modify X and y to stay away from extreme values (for 'vfa' data where there are a lot of perfect separations)
%     v = unique(X(:,3));
%     for j = 1:length(v)
%         if var(y(X(:,3)==v(j))) == 0
%             X = [X; [0, 0, v(j),0, 0]; [0, 0, v(j),0,0]];
%             y = [y; 0; 1];
%         end
%     end
    
    [b, dev, stats] = glmfit(X, y, 'binomial', 'link', 'logit');
    model = [b',dev];
    M(i,:) = model; 
    indivpVals(i,:) = [rot90(stats.p), 1]; % the 1 is added for later convenience of marking the significant coefficient
end

% label the coefficients that are significant from the calculation
parameters = cell(size(M));
for ii=1:size(parameters,1)
    for jj=1:size(parameters,2)
        if indivpVals(ii,jj) < 0.05
            % convert the values into string, assign the value to tTestDecisions{i,j}.
            parameters{ii,jj} = strcat(num2str(M(ii,jj)),' (*)');
        else
            parameters{ii,jj} = num2str(M(ii,jj));
        end
    end
end

T = array2table(parameters,'VariableNames',{'Constant_term','Elevation','Azimuth','Previous_perceived_rotation_direction','Previous_perceived_facing_direction','Previous_perceived_viewing_perspective','Deviance_of_the_fit'});   
disp('=====Model Coefficient=====')
disp('(each line is the model parameter for an individual dataset)')
disp('((*) indicates the predictor is significant in that individual dataset)')
disp(T);

avg = mean(M,1);
avgT = array2table(avg,'VariableNames',{'Constant_term','Elevation','Azimuth','Previous_perceived_rotation_direction','Previous_perceived_facing_direction','Previous_perceived_viewing_perspective','Deviance_of_the_fit'});   
disp('=====Mean model coefficient=====')
disp(avgT)

% H tells you which predictor is significant across the whole set
[H,p,ci,stats] = ttest(M(:,1:end-1),0);
disp('===Across all datasets, which predictor is significant? (1 is significnat, 0 is not)')
HT = array2table(H,'VariableNames',{'Constant_term','Elevation','Azimuth','Previous_perceived_rotation_direction','Previous_perceived_facing_direction','Previous_perceived_viewing_perspective'});   
disp(HT)

varInX = [{'Intercept'} varInX];
M = array2table(M(:,1:end-1),'VariableNames',varInX);
