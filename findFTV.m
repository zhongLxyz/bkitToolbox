% This function finds the trials in which the responses are FTV.
% input: angularVelocity, azimuth, elevation and response columns in the dataset. The
% response column can be coded in 0 and 1, or -1 and 1. Angular Velocity
% can be coded in 45 and -45, or 1 and -1. 
% output: a column vector in which FTV=1 and FA=0 for corresponding trials.

% The last varargin codes for what the original
% column codes for, it can be 'vfa', 'cw', or 'ftv'.

% It can also process outputs from PsychBench 
% (last modified: June 14, 2017)

function ftv = findFTV(angularVel,azimuth,response,elevation,varargin)

if nargin < 5
    type = 'ccw'; % cw/ccw and left/right can be converted to ftv in the same way
else
    type = varargin;
end

if ismember(0,response) % if response is coded in 0 and 1
   response = (response-0.5)*2; % turn into -1 and 1
end

if strcmp(type,'ccw')
    ftv = cw(angularVel, azimuth, response);
elseif strcmp(type,'ftv')
    response = -1 * response;
    ftv = response/2+0.5;
elseif strcmp(type,'vfa')
    ftv = vfa(azimuth,response, elevation);
end
end

function ftv = cw(angularVel, azimuth, response)

ftv = zeros(length(response),1);
correct = response==sign(angularVel);

% for azimuth angles in the front circle, if the response is the same as
% the actual spinning direction, then it is FTV.
ftv(azimuth>-90 & azimuth<90 & correct) = 1;
% for -azimuth angles in the back circle, if the response is the opposite as
% the actual spinning direction, then it is FTV.
ftv((azimuth<-90 & ~correct) | (azimuth>90 & ~correct)) = 1;

% for 90 and -90, the perceived facing direction is the same as the
% response.
ftv(azimuth==-90 & response==1) = 1;
ftv(azimuth==90 & response==-1) = 1;
end

function ftv = vfa(azimuth,response,elevation)

ftv = zeros(length(response),1);
correct = response~=sign(elevation);

% for azimuth angles in the front circle, if the response is the same as
% the actual spinning direction, then it is FTV.
ftv(azimuth>-90 & azimuth<90 & correct) = 1;
% for -azimuth angles in the back circle, if the response is the opposite as
% the actual spinning direction, then it is FTV.
ftv((azimuth<-90 & ~correct) | (azimuth>90 & ~correct)) = 1;

% for 90 and -90, the perceived facing direction is the same as the
% response.

% is this required???
ftv(azimuth==-90 & response==1) = 1;
ftv(azimuth==90 & response==-1) = 1;
%

end
