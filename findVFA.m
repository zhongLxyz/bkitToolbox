function vfa = findVFA(angularVel,azimuth,response,elevation,varargin)
% This function finds the trials in which the responses are VFA.
% input: angularVelocity, elevation and response columns in the dataset. The
% response column can be coded in 0 and 1, or -1 and 1 . Angular Velocity
% can be coded in 45 and -45, or 1 and -1. 
% optional argument: what the response column is coding for.
% 'cw' = cw/ccw, 'ftv'=perceived facing direction, 'vfa' = perceived
% elevation.
% output: a column vector in which VFA=1 and Viewed-from-below=0 for corresponding trials.
% It can also process outputs from PsychBench 

if nargin < 5
    type = 'ccw'; % cw/ccw and left/right can be converted to vfa in the same way
else
    type = varargin;
end

if ismember(0,response) % if response is coded in 0 and 1
   response = (response-0.5)*2; % turn into -1 and 1
end

if strcmp(type,'ccw')
    vfa = cw(angularVel, elevation, response);
elseif strcmp(type,'vfa')
    response = -1 * response;
    vfa = response/2+0.5;
elseif strcmp(type,'ftv')
    vfa = ftv(elevation, azimuth, response);
end
end

function vfa = cw(angularVel, elevation, response)

correct = response==sign(angularVel);
vfa = zeros(length(elevation));
vfa(correct & elevation>0) = 1;
vfa(~correct & elevation<0) = 1;
vfa = vfa';

end

function vfa = ftv(elevation, azimuth, response)

facingDir = azimuth; % the actual facing direction of the walker
facingDir(facingDir>-90 & facingDir<90)=1;
facingDir(facingDir<-90 | facingDir>90)=-1;
correct = facingDir~=response;
vfa = zeros(length(azimuth));
vfa(correct & elevation>0)=1;
vfa(~correct & elevation<0)=1;
vfa = vfa';

end