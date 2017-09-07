function [data] = bkitRead(name)
%% this function reads in data files from all participants (sessions), the output is in
% the forms of a cell array, with each cell containing
% a table of relevant data for one participant.
% input the csv file name or the folder name to process the data in that folder.
%
% if the input is a folder, then it assumes each file in the folder is
% the data of an individual participant; if the input is a file (can 
% contain the path, or just a name in the current folder), it
% assumes the csv file contains data of all participants (with different
% session IDs in the first column).if the input is a file then it must 
% contain its extension.
%
% If the design uses 20, 40... and their corresponding rotation direction
% to represent 30, 0... degree azimuth, this function
% automatically converts the azimuths to 30, 0, -30 degrees for the output.
% (for BMLkit: For responses coded in '0i' and 'i0', 'i0' is recoded as 1, which means
% the response corresponding to the left arrow key is 1)

[pathstr,filename,ext] = fileparts(name);

if isempty(ext) % if the extension is none, so 'name' is a folder name.
    
    % locate the folder and find all data files in it
    
    lookFor = [name, '\', '*.csv'];
	if isempty(lookFor)
		lookFor = [name, '\', '*.xlsx']; % some files are in .xlsx format
	end
    files = dir(lookFor);
	
	% decide if the data is output from the PsychBench or the BMLkit by seeing how many columns there are
	firstDataDir = [name, '\', files(1).name];
	firstData = readtable(firstDataDir);
    % if the data is from Psychbench (the first variable name has 'trial' in
    % it)
	if ~isempty(strfind(firstData.Properties.VariableNames{1},'trial'))
		% if the files are the data output from the Psychbench
		for i = 1:length(files)
			% read the csv file and turn it into a table
			filename = files(i).name;
			fullName = [name, '\', filename];
			singleSubj = readtable(fullName,'ReadVariableNames',false);
			
			% change the header of the columns into correct ones
			singleSubj.Properties.VariableNames{'Var2'} = 'Walkers';
			singleSubj.Properties.VariableNames{'Var3'} = 'CameraAzimuth';
			singleSubj.Properties.VariableNames{'Var4'} = 'AngularVelocity';
			singleSubj.Properties.VariableNames{'Var5'} = 'CameraElevation';
			singleSubj.Properties.VariableNames{'Var7'} = 'Answer';
			singleSubj.Properties.VariableNames{'Var10'} = 'ResponseLatency';
			
			% delete the unnecessary rows
			singleSubj(1:2,:) = [];
            % delete rows that have blanks in the azimuth slot
            singleSubj(strcmp(singleSubj.CameraAzimuth,'[]')==1,:) = []; 
			
			% if azimuth is other numbers instead of 30, 0..., convert it
			% back to the usual way of representing azimuth.
			%singleSubj.CameraAzimuth = psyBenConvertAzimuth(singleSubj.CameraAzimuth);
			
			% recode the responses as 0 and 1 (or -1 and 1)
			Response = recodeResponse(singleSubj.Answer);
            Response = table(Response);
            
            % add a column of session number (or just the file name)
            session = cell(size(Response));
            for ii=1:length(session)
                session(:) = {filename};
            end
            session = table(session,'VariableName',{'Session'});
            
            % convert data type
            singleSubj.Walkers = cellfun(@str2double,singleSubj.Walkers);
            if iscell(singleSubj.ResponseLatency(1))
                singleSubj.ResponseLatency = cellfun(@str2double,singleSubj.ResponseLatency);
            end
            singleSubj.CameraAzimuth = cellfun(@str2double,singleSubj.CameraAzimuth);
            singleSubj.AngularVelocity = cellfun(@str2double,singleSubj.AngularVelocity);
            singleSubj.CameraElevation = cellfun(@str2double,singleSubj.CameraElevation);
			
			% keep the relevant columns for later analysis
			cleanData = [session, singleSubj(:,{'Walkers','ResponseLatency','CameraAzimuth','AngularVelocity','CameraElevation'}), Response];
			cleanData = {cleanData};
			% output data is a cell array with each cell containing a table of data
			% for each csv file.
			data(i) = cleanData;
			trialCounts(i) = height(Response);
		end
	else
		% if the files are the data output from the BMLkit 
		for i = 1:length(files)
			% read the csv file and turn it into a table
			filename = files(i).name;
			fullName = [name, '\', filename];
			singleSubj = readtable(fullName);
			
			catchTrials = singleSubj.Object_0_walker_1_box_ ~= 0;
			singleSubj(catchTrials,:) = [];
			
			% change the header of the answer column if needed
			if ismember('Answer_br__',singleSubj.Properties.VariableNames)
				singleSubj.Properties.VariableNames{'Answer_br__'} = 'Answer';
			end
			
			% if azimuth uses 20, 40... instead of 30, 0..., convert it
			% back to the usual way of representing azimuth.
			singleSubj.CameraAzimuth = bkitConvertAzimuth(singleSubj.CameraAzimuth);
			
			Response = recodeResponse(singleSubj.Answer);
			Response = table(Response);
			
			% keep the relevant columns for later analysis
			cleanData = [singleSubj(:,{'SessionId','Gender','Speed','CameraAzimuth','AngularVelocity','CameraElevation'}), Response];
			cleanData = {cleanData};
			% output data is a cell array with each cell containing a table of data
			% for each csv file.
			data(i) = cleanData;
			trialCounts(i) = height(Response);
        end
    end
    %% if the provided 'name' is a csv file containing all participants' data.
else 
    allData = readtable(name);
    ids = unique(allData.SessionId);
    i = 1; 
    
    for k = 1:length(ids)
        subjId = ids(k);
        singleSubj = allData(allData.SessionId == subjId,:);
        catchTrials = singleSubj.Object_0_walker_1_box_ ~= 0;
        singleSubj(catchTrials,:) = [];
        
        % if there are no trials left after deleting the catchTrials, do
        % not continue into the reading process.
        if height(singleSubj) > 40
            % change the header of the answer column if needed
            if ismember('Answer_br__',singleSubj.Properties.VariableNames)
                singleSubj.Properties.VariableNames{'Answer_br__'} = 'Answer';
            end
            
            % if azimuth uses 20, 40... instead of 30, 0..., convert it
            % back to the usual way of representing azimuth.
            singleSubj.CameraAzimuth = bkitConvertAzimuth(singleSubj.CameraAzimuth);
            
            Response = recodeResponse(singleSubj.Answer);
            Response = table(Response);
            
            % keep the relevant columns for later analysis
            cleanData = [singleSubj(:,{'SessionId','Gender','Speed','CameraAzimuth','AngularVelocity','CameraElevation'}), Response];
            cleanData = {cleanData};
            
            data(i) = cleanData;
            trialCounts(i) = height(Response);
            i = i + 1;
        end
    end
end

% delete data that has too few trials in them.
limitSize = mode(trialCounts)-5;
for i = length(data):-1:1
    if height(data{i}) < limitSize
        data(i) = [];
    end
end

function converted = bkitConvertAzimuth(origAzimuth)
% this function convert the azimuth back to -30, 0 and 30 in Zhong's 3rd
% experiment.
% input: origAzimuth - original column of azimuth for one participant
% output: converted - the converted column of azimuth

origAzimuth(origAzimuth==-40 | origAzimuth==-20) = -30;
origAzimuth(origAzimuth==-10| origAzimuth==10) = 0;
origAzimuth(origAzimuth==20| origAzimuth==40) = 30;
origAzimuth(origAzimuth==-140| origAzimuth==-160) = -150;
origAzimuth(origAzimuth==-170| origAzimuth==170) = 180;
origAzimuth(origAzimuth==160| origAzimuth==140) = 150;
converted = origAzimuth;

