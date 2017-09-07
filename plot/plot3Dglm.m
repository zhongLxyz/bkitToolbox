% create 3D bar graph of model parameter values for each participant. 
% (in this version, both significant and non-significant values will be
% shown on the plot, however, the nonsignificant columns will be plotted in
% a different color)
% (all the relevant variables for the glm in the August experiemnt are stored in 'forPlot3DgraphGLM.mat' in experiment3 folder)

% reading in the datasets
% ** please change the name of the 3 .csv files to correspond to each
% experiment
cw = bkitRead('cwEx.csv');
facing = bkitRead('facingEx.csv');
above = bkitRead('aboveEx.csv');

% finding model parameters and p-values (whether the factors for each participant
% are significant or not) for each condition
[Mcw, ~, ~, indivpValsCW, ~] = bkitModel(cw,'cw');
[Mfacing, ~, ~, indivpValsFacing, ~] = bkitModel(facing,'ftv');
[Mabove, ~, ~, indivpValsAbove, ~] = bkitModel(above,'vfa');

% turn parameter values that are not significant into zero.
% Mcw(indivpValsCW>0.05) = 0;
% Mfacing(indivpValsFacing>0.05) = 0;
% Mabove(indivpValsAbove>0.05) = 0;

figure1 = figure;
for i=1:length(cw)
    
    % for no history factors: order the columns as ccw bias, FTV bias and VFA bias
    bar = zeros(3,3);
    bar(1,:) = [Mcw(i,1),Mcw(i,3),Mcw(i,2)];
    bar(2,:) = [Mfacing(i,1),Mfacing(i,3),Mfacing(i,2)];
    bar(3,:) = [Mabove(i,1),Mabove(i,3),Mabove(i,2)];
    disp(bar)
    
    % whether the values are significant.
    sig = [indivpValsCW(i,1), indivpValsCW(i,3), indivpValsCW(i,2);
           indivpValsFacing(i,1), indivpValsFacing(i,3), indivpValsFacing(i,2);
           indivpValsAbove(i,1), indivpValsAbove(i,3), indivpValsAbove(i,2)];
    disp(i)
    disp(sig<0.05);
       
    % row and column that has sig values
%     [r,c]=find(sig<0.05);
%     r=r';
%     c=c';
%     pos = find(sig<0.5);
       
%     % for history factors
%     bar = zeros(3,6);
%     bar(1,:) = [Mcw(i,1),Mcw(i,3),Mcw(i,2),Mcw(i,4:6)];
%     bar(2,:) = [Mfacing(i,1),Mfacing(i,3),Mfacing(i,2),Mfacing(i,4:6)];
%     bar(3,:) = [Mabove(i,1),Mabove(i,3),Mabove(i,2),Mfacing(i,4:6)];
%%  
    b = bar3(bar);
    colorbar
    title(sprintf('participant %d strength of biases in each condition',i))
    zlim([-3 10])
    caxis([-3 10])
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

    %disp('i')
    
    % save the figures as pictures
    %saveas(figure1,sprintf('glm%d.jpg', i));
    input('press ENTER to continue')
end

% find the average of the effect for the whole experiments.

% calculate the mean parameter values for each condition across
% participants.
meanMcw = mean(Mcw,1);
meanMfacing = mean(Mfacing,1);
meanMabove = mean(Mabove,1);

% arrange them into a 3*3 matrix for plotting
bar = zeros(3,3);
bar(1,:) = [meanMcw(1),meanMcw(3),meanMcw(2)];
bar(2,:) = [meanMfacing(1),meanMfacing(3),meanMfacing(2)];
bar(3,:) = [meanMabove(1),meanMabove(3),meanMabove(2)];

b = bar3(bar);
colorbar
title(sprintf('mean strength of biases in each condition across participants',i))
zlim([-5 20])
caxis([-5 20])
set(gca,'XTickLabel',{'CCW bias';'FTV bias';'VFA bias'}) % for x direction
set(gca,'YTickLabel',{'CCW session';'FTV session';'VFA session'}) % for y direction

% displaying the values of parameters on top of each bar
textposi = [bar(1,:)+0.1,bar(2,:)+0.16,bar(3,:)+0.16];
values = strsplit(num2str([bar(1,:),bar(2,:),bar(3,:)]));
text([1,2,3,1,2,3,1,2,3],[1,1,1,2,2,2,3,3,3],textposi,values)
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

% display the values mean values along with the sd as a table

% find sd values for the 3 biases in each condition (arranged in a matrix 
% with the same position as the average matrix)
Sd = [std(Mcw(:,1)),std(Mcw(:,3)),std(Mcw(:,2)); 
      std(Mfacing(:,1)),std(Mfacing(:,3)),std(Mfacing(:,2)); 
      std(Mabove(:,1)),std(Mabove(:,3)),std(Mabove(:,2))]; 

meanStd = cell(3,3);
for i=1:size(bar,1)
    for j=1:size(bar,2)
        % convert the sd value into string, assign the value to sdAll{i,j}.
        meanStd{i,j} = strcat(num2str(bar(i,j)),' (',num2str(Sd(i,j)),')');
    end
end

T = cell2table(meanStd,'VariableNames',{'CW_bias','FTV_bias','VFA_bias'}, 'RowNames', {'CW_session','FTV_session','VFA_session'});;
disp(T)

% t-test for comparing the difference between VFA responses from the FTV
% and VFA sessions.
h1 = ttest(Mfacing(:,2),Mabove(:,2))
% t-test for comparing the difference between FTV responses from the FTV
% and VFA sessions.
h2 = ttest(Mfacing(:,3),Mabove(:,3))
