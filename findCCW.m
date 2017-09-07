% This function finds the trials in which the responses are CCW.
% input: angularVelocity, azimuth, elevation and response columns in the dataset. The
% response column can be coded in 0 and 1, or -1 and 1, Angular Velocity
% can be coded in 45 and -45, or 1 and -1.
% output: a column vector in which FTV=1 and FA=0 for corresponding trials.

% optional argument: what the response column is coding for.
% 'cw' = cw/ccw, 'ftv'=perceived facing direction, 'vfa' = perceived
% elevation.

% It can also process outputs from PsychBench 
% (last modified: June 14, 2017)

function ccw = findCCW(angularVel,azimuth,response,elevation,varargin)

if nargin < 5
    type = 'ccw'; % cw/ccw and left/right can be converted to ftv in the same way
else
    type = varargin{1};
end

if ismember(0,response) % if response is coded in 0 and 1
   response = (response-0.5)*2; % turn into -1 and 1
end

if strcmp(type,'ccw')
    ccw = response/2+0.5;
elseif strcmp(type,'ftv')
    ccw = ftv(angularVel, azimuth, response);
elseif strcmp(type,'vfa')
    ccw = vfa(angularVel,response, elevation);
end

end

function ccw = ftv(angularVel, azimuth, response)

ccw = zeros(length(response),1);

facingDir = azimuth; % the actual facing direction of the walker
facingDir(facingDir>-90 & facingDir<90)=1;
facingDir(facingDir<-90 | facingDir>90)=-1;
correct = facingDir~=response;

ccw(correct & sign(angularVel)>0)=1;
ccw(~correct & sign(angularVel)<0)=1;
end

function ccw = vfa(angularVel,response,elevation)

ccw = zeros(length(response),1);
correct = response~=sign(elevation);

ccw(correct & sign(angularVel)>0)=1;
ccw(~correct & sign(angularVel)<0)=1;

end
