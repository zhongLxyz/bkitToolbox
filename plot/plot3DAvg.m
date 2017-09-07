% This script does the following things:
% 1. creating 3D bar graph of bias size for each participant in each condition (3*3 bar graph).
%    including an average plot for all participants.
% 2. printing a table containing the average bias size across all participants for each conditions,
%    along with their SD
% 3. print a table with t-test results comparing whether there is a significant difference between
%    the size of (the same) biases in different conditions.

% reading in the datasets
% ** please change the name of the 3 .csv files to correspond to each
% experiment
%cw = bkitRead('cw.csv');
%facing = bkitRead('facing.csv');
%above = bkitRead('above.csv');

cw = bkitRead('CCW');
facing = bkitRead('FTV');
above = bkitRead('VFA');

% initialize a 3D matrix to store the data for one participant on each layer. 
all = zeros(3,3,length(cw));

figure1 = figure;
for i=1:length(cw)
    
    fprintf('%s: \n',cw{i}.Session{1});
    
    % for no history factors: order the columns as ccw bias, FTV bias and VFA bias
    [angularVelCW, azimuthCW, responseCW, elevationCW] = returnvalues(cw{i});
    [angularVelFacing, azimuthFacing, responseFacing, elevationFacing] = returnvalues(facing{i});
    [angularVelAbove, azimuthAbove, responseAbove, elevationAbove] = returnvalues(above{i});
    
    % finding the size of each bias (# of responses/all)
    bar = zeros(3,3);
    % for all the cw bias
    bar(1,1) = mean(findCCW(angularVelCW, azimuthCW, responseCW, elevationCW, 'ccw'));
    bar(2,1) = mean(findCCW(angularVelFacing, azimuthFacing, responseFacing, elevationFacing, 'ftv'));
    bar(3,1) = mean(findCCW(angularVelAbove, azimuthAbove, responseAbove, elevationAbove, 'vfa'));
    
    % for all the FTV bias
    bar(1,2) = mean(findFTV(angularVelCW, azimuthCW, responseCW, elevationCW, 'ccw'));
    bar(2,2) = mean(findFTV(angularVelFacing, azimuthFacing, responseFacing, elevationFacing,'ftv'));
    bar(3,2) = mean(findFTV(angularVelAbove, azimuthAbove, responseAbove, elevationAbove,'vfa'));
    
    % for all the VFA bias
    bar(1,3) = mean(findVFA(angularVelCW, azimuthCW, responseCW, elevationCW, 'ccw'));
    bar(2,3) = mean(findVFA(angularVelFacing, azimuthFacing, responseFacing, elevationFacing,'ftv'));
    bar(3,3) = mean(findVFA(angularVelAbove, azimuthAbove, responseAbove, elevationAbove,'vfa'));
    
    barT = array2table(bar,'VariableNames',{'CCW_bias','FTV_bias','VFA_bias'}, 'RowNames', {'CCW_session','FTV_session','VFA_session'});
    disp(barT)
    % store the "bar" for mean analysis later
    all(:,:,i) = bar;

    %% plotting the figure, comment this section out if you don't need this
    b = bar3(bar);
    colorbar
    title(sprintf('Day %d strength of biases in each condition',i))
    zlim([0 1])
    caxis([0 1])
    set(gca,'XTickLabel',{'CCW bias';'FTV bias';'VFA bias'}) % for x direction
    set(gca,'YTickLabel',{'CCW session';'FTV session';'VFA session'}) % for y direction
    % for history factors
    %     set(gca,'XTickLabel',{'CCW bias';'FTV bias';'VFA bias'}) % for x direction
    %     set(gca,'YTickLabel',{'CCW session';'FTV session';'VFA session';'preCCW','preFTV','prevVFA'}) % for y direction
    
    % displaying the values of parameters on top of each bar
    textposi = [bar(1,:)+0.1,bar(2,:)+0.16,bar(3,:)+0.16];
    values = strsplit(num2str([bar(1,:),bar(2,:),bar(3,:)]));
    text([1,2,3,1,2,3,1,2,3],[1,1,1,2,2,2,3,3,3],textposi,values)
    %     for ii=1:length(r)
    %         text(r(ii),c(ii),textposi(pos(ii)),num2str(bar(r(ii),c(ii))),'color','red');
    %     end
    % color the sig ones differently.
    %text(r,c,textposi,values,'color','red');
    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
    end
    
    % save the figures as pictures
    saveas(figure1,sprintf('S001_Day%d.jpg', i));
    input('press ENTER to continue')
end

% find the average of the effect for the whole experiments.
avg = mean(all,3);

% plot it
b = bar3(avg);
colorbar
title(sprintf('Mean strength of biases in each condition across days',i))
zlim([0 1])
caxis([0 1])
set(gca,'XTickLabel',{'CCW bias';'FTV bias';'VFA bias'}) % for x direction
set(gca,'YTickLabel',{'CCW session';'FTV session';'VFA session'}) % for y direction

% displaying the values of parameters on top of each bar
textposi = [avg(1,:)+0.1,avg(2,:)+0.16,avg(3,:)+0.16];
values = strsplit(num2str([avg(1,:),avg(2,:),avg(3,:)]));
text([1,2,3,1,2,3,1,2,3],[1,1,1,2,2,2,3,3,3],textposi,values)
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
saveas(figure1,sprintf('P1avgMeanPlot.jpg'));

% find sd values for the 3 biases in each condition (arranged in a matrix
% with the same position as the average matrix)
Sd = [std(all(1,1,:)),std(all(1,2,:)),std(all(1,3,:));
    std(all(2,1,:)),std(all(2,2,:)),std(all(2,3,:));
    std(all(3,1,:)),std(all(3,2,:)),std(all(3,3,:))];

% display the values mean values along with the sd as a table
meanStd = cell(3,3);
for i=1:size(avg,1)
    for j=1:size(avg,2)
        % convert the sd value into string, assign the value to sdAll{i,j}.
        meanStd{i,j} = strcat(num2str(avg(i,j)),' (',num2str(Sd(i,j)),')');
    end
end

T = cell2table(meanStd,'VariableNames',{'CCW_bias','FTV_bias','VFA_bias'}, 'RowNames', {'CCW_session','FTV_session','VFA_session'});;
fprintf('\nAverage bias size and their standard deviation:\n')
disp(T)

% put all p-values of the t-test decisions in a matrix
p = zeros(3,3);

% CCW responses
[h1,p(1,1),ci1,stats1] = ttest(all(1,1,:),all(2,1,:));  % cw and ftv condition
[h2,p(2,1),ci2,stats2] = ttest(all(1,1,:),all(3,1,:));  % cw and vfa condition
[h3,p(3,1),ci3,stats3] = ttest(all(2,1,:),all(3,1,:));  % ftv and vfa condition

% FTV bias
[h4,p(1,2),ci4,stats4] = ttest(all(1,2,:),all(2,2,:));  % cw and ftv condition
[h5,p(2,2),ci5,stats5] = ttest(all(1,2,:),all(3,2,:));  % cw and vfa condition
[h6,p(3,2),ci6,stats6] = ttest(all(2,2,:),all(3,2,:));  % ftv and vfa condition

% VFA bias
[h7,p(1,3),ci7,stats7] = ttest(all(1,3,:),all(2,3,:));  % cw and ftv condition
[h8,p(3,3),ci8,stats8] = ttest(all(1,3,:),all(3,3,:));  % cw and vfa condition
[h9,p(3,3),ci9,stats9] = ttest(all(2,3,:),all(3,3,:));  % ftv and vfa condition

% put all t values in a matrix
tValues = [stats1.tstat,stats4.tstat,stats7.tstat;
			stats2.tstat,stats5.tstat,stats8.tstat;
			stats3.tstat,stats6.tstat,stats9.tstat;];

tTestDecisions = cell(3,3);
for ii=1:size(tTestDecisions,1)
    for jj=1:size(tTestDecisions,2)
        % convert the values into string, assign the value to tTestDecisions{i,j}.
        tTestDecisions{ii,jj} = strcat(num2str(tValues(ii,jj)),' (',num2str(p(ii,jj)),')');
    end
end
T2 = cell2table(tTestDecisions, 'VariableNames', {'CCW_bias','FTV_bias','VFA_bias'},'RowNames', {'CCW and FTV condition','CCW and VFA condition','FTV and VFA condition'});
fprintf('\nT-test results with p-values comparing size of bias between conditions:\n')
disp('(the first number is the t-value, the number in brackets is the p-value)')
disp(T2)

%% find the mean size of each bias for each participants

% plot bar graph for the ranges of the mean values (i.e. how many people in each group)

% get mean ccw bias for each participant
meanCCWeach = mean(all(:,1,:));
meanCCWeach = squeeze(meanCCWeach); % make it into a 1D vector
% do the same things for the FTV and VFA bias
meanFTVeach = squeeze(mean(all(:,2,:)));
meanVFAeach = squeeze(mean(all(:,3,:)));

% find how varied each participant' size of the ccw bias is in 3 conditions
% (find the standard deviation of the same participant, of the same bias)
CCWeach = squeeze(all(:,1,:));
meanStdCWeach = mean(std(CCWeach));
% do the same things for the other biases
meanStdFTVeach = mean(std(squeeze(all(:,2,:))));
meanStdVFAeach = mean(std(squeeze(all(:,3,:))));

