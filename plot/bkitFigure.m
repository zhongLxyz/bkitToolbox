% This function plots the relevant graphs for bmlkit data.
% arguments: 
% how to use it: the first argument has to be the dataset read from bkitRead. The other arguments are optional,
% including 
% (1) the name of the plot you want. 'FTVelevation', 'VFAelevation', 'VFAazimuth' for the 3 kinds of plots.
%   FTVelevation: for each session, plot the percentage of FTV responses at each elevation angles.
%   VFAelevation: for each session, plot the percentage of VFA responses at each elevation angles.
%   VFAazimuth: for each session, plot the percentage of VFA responses at each azimuth angles.
% (2) specification of what kind of responses is coded in the response column of the dataset.
% 'ftv', 'cw', 'vfa' - it specifies what the kind of responses is coded in the response column.

function bkitFigure(allData, varargin)

% 'type' specify what responses is coded in the response column.
if ismember('ccw',varargin)
    type = 'ccw';
elseif ismember('ftv',varargin)
    type = 'ftv';
elseif ismember('vfa',varargin)
    type = 'vfa';
else 
    type = 'ccw';
end

for j = 1:length(allData)
    data = allData{j}; % data stores a table of values for one participant
    elevation = data.CameraElevation;
    azimuth = data.CameraAzimuth;
    angularVel = data.AngularVelocity;
    response = (data.Response - 0.5) * 2; % turn into -1 and 1
    if iscell(data.SessionId(1))
        id = data.SessionId{1};
    else
        id = data.SessionId(1);
    end
    disp(id)
    
    %% plot the rate of FTV responses vs. the elevation
    if ismember('FTVelevation', varargin)
        v_elev = unique(elevation);
        for i = 1:length(v_elev) % for every elevation value
            % collect the number of FTV responses for it
            trialsNum = elevation == v_elev(i); % trial #s having this elevation
            percFTV(i) = mean(findFTV(angularVel(trialsNum), azimuth(trialsNum), response(trialsNum), elevation(trialsNum),type));
        end
        disp(j)
        disp(percFTV)
        plot(v_elev,percFTV);
        hold on
        scatter(v_elev,percFTV);
        
        if ischar(id)
            title(sprintf('FTV bias at each elevation angle for participant %s',id))
        else
            title(sprintf('FTV bias at each elevation angle for participant %d',id))
        end
        xlabel('elevation')
        ylim([0 1])
        ylabel('percentage of FTV responses')
        hold off
%         if ischar(id)
%             saveas(figure,sprintf('FTVelevation%s.jpg', id));
%         else
%             saveas(figure,sprintf('FTVelevation%d.jpg', id));
%         end
        input('press ENTER to continue')
    end
    
    %% Plot the VFA bias (% of VFA response at each elevation angles)
    if ismember('VFAelevation', varargin)
        v_elev = unique(elevation);
        for i = 1:length(v_elev)
            trialsNum = elevation == v_elev(i); % trial #s having this elevation
            percVFA(i) = mean(findVFA(angularVel(trialsNum), azimuth(trialsNum), response(trialsNum), elevation(trialsNum),type));
        end
        plot(v_elev, percVFA)
        hold on
        scatter(v_elev, percVFA)
        if ischar(id)
            title(sprintf('VFA bias at each elevation angle for participant %s',id))
        else
            title(sprintf('VFA bias at each elevation angle for participant %d',id))
        end
        xlabel('elevation')
        ylim([0 1])
        ylabel('percentage of VFA responses')
        hold off
%         if ischar(id)
%             saveas(figure,sprintf('VFAelevation%s.jpg', id));
%         else
%             saveas(figure,sprintf('VFAelevation%d.jpg', id));
%         end
        input('press ENTER to continue');
    end
    
    %% Plot the % of VFA responses at each azimuth angle
    if ismember('VFAazimuth', varargin)
        v_azimuth = unique(azimuth);
        for i = 1:length(v_azimuth)
            trialsNum = azimuth == v_azimuth(i); % trial #s having this elevation
            percVFA(i) = mean(findVFA(angularVel(trialsNum), azimuth(trialsNum), response(trialsNum), elevation(trialsNum),type));
        end
        plot(v_azimuth, percVFA)
        hold on
        scatter(v_azimuth, percVFA)
        ylim([0 1])
        if ischar(id)
            title(sprintf('VFA bias at each azimuth angle for participant %s',id))
        else
            title(sprintf('VFA bias at each azimuth angle for participant %d',id))
        end
        hold off
%         if ischar(id)
%             saveas(figure,sprintf('VFAazimuth%s.jpg', id));
%         else
%             saveas(figure,sprintf('VFAazimuth%d.jpg', id));
%         end
        input('press ENTER to continue');
    end
end
