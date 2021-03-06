% This script finds whether there is any serial dependency effect in the input dataset.
% It will conduct the test in all 3 variables (cw/ccw, ftv/fa, vfa/vfb)
% It conducts a t-test for the data of each session, to determine if
% repeating answers happen more often than chance would predict.

% things to input before doing the analysis
% Input: N - an integer N, it represents the N-back response from the
% current response.
% type - it specifies what the response column of the input dataset is coding for.
% ('ftv', 'cw' or 'vfa').
% displays: out - for each session, it outputs a row containing the
% probablity values arranged in the following order:
% pA, pAA/(pA*pA),pBA/(pB*pA),pBB/(pB*pB),pAB/(pA*pB);
% (pAA/(pA*pA) is the probabilty of having repeated A responses divided by
% the product of the prob of having an A response multiplied by itself)

% modify this section to specify the dataset:
allData = bkitRead('VFA');
type = 'vfa';  % what the answer column of the dataset is coding for
N = 1;  % N-back response;

disp('serial dependency effect in the CW/CCW variable')
for i=1:length(allData)
    singleSubj = allData{i};
    fprintf('%s:',singleSubj.Session{1});

    response = singleSubj.Response;           % 0: clockwise,  1: ccw
    angularVel = singleSubj.AngularVelocity;     % spinning direction
    azimuth = singleSubj.CameraAzimuth;
    elevation = singleSubj.CameraElevation;             % camera elevation

    dCCW = findCCW(angularVel,azimuth,response,elevation,type); 
    n = length(dCCW);

    pA = sum(dCCW)/n;             % prob of ccw
    pB = 1-pA;                 % prob of cw
	
    % divide by the total number of sequences, so you get the prob of sequence of 1,1 or 1,0 etc.
    pAA = sum(dCCW(N+1:end)==1 & dCCW(1:end-N)==1)/(n-N);   % frequency of 1 then 1
    pBA = sum(dCCW(N+1:end)==1 & dCCW(1:end-N)==0)/(n-N);   % frequency of 0 then 1
    pBB = sum(dCCW(N+1:end)==0 & dCCW(1:end-N)==0)/(n-N);   % frequency of 0 then 0
    pAB = sum(dCCW(N+1:end)==0 & dCCW(1:end-N)==1)/(n-N);   % frequency of 1 then 0    
    
    out(i,:) = [pA, pAA/(pA*pA),pBA/(pB*pA),pBB/(pB*pB),pAB/(pA*pB)];
    disp(out(i,:));
end

AA = out(:,2);
BB = out(:,4);
AB = out(:,5);
BA = out(:,3);


disp('Averaged values across sessions:')
disp('pA, pAA/(pA*pA),pBA/(pB*pA),pBB/(pB*pB),pAB/(pA*pB)')
m = mean(out); % return mean of each column of out (average the quantities across the sessions)
disp(m);

% for each session, for AA and BB series, calculate their mean along
% each row, and compare (do t-test) with the mean in AB and BA series.
[h,p,ci,stats] = ttest(mean(out(:,[2,4]),2),mean(out(:,[3,5]),2)); % p is p-value

mean1 = mean(mean(out(:,[2,4]),2),1); % mean of AA and BB series
mean2 = mean(mean(out(:,[3,5]),2),1); % AB and BA series
std1 = std(mean(out(:,[2,4]),2)); % standard deviation
std2 = std(mean(out(:,[3,5]),2));

fprintf('average strength of repetition((pAA/(pA*pA)+pBB/(pB*pB))/2) across sessions: %d\n', mean1);
fprintf('standard deviation: %d\n', std1);
fprintf('average strength of alternation((pAB/(pA*pB)+pBA/(pB*pA))/2) across sessions: %d\n', mean2);
fprintf('standard deviation: %d\n', std2);
fprintf('t-test results for comparing whether there is significant difference between the two: t-value %d, p-value: %d\n',stats.tstat,p)
eachComp = [mean1; std1; mean2; std2; stats.tstat; p];

disp(' ')
disp('-------------------------------------------------')
disp('Serial dependency effect in the FTV/FA variable')
for i=1:length(allData)
    singleSubj = allData{i};
    fprintf('%s:',singleSubj.Session{1});

    response = singleSubj.Response;           % 0: clockwise,  1: ccw
    angularVel = singleSubj.AngularVelocity;     % spinning direction
    azimuth = singleSubj.CameraAzimuth;
    elevation = singleSubj.CameraElevation;             % camera elevation

    % finding perceived facing Dir 
    d = findFTV(angularVel,azimuth,response,elevation,type); % d is the facing direction perceived. 1: FTV, 0: FA
    n = length(d);

    pA = sum(d)/n;             % prob of FTV
    pB = 1-pA;                 % prob of FA
	
    % divide by the total number of sequences, so you get the prob of sequence of 1,1 or 1,0 etc.
    pAA = sum(d(N+1:end)==1 & d(1:end-N)==1)/(n-N);   % frequency of 1 then 1
    pBA = sum(d(N+1:end)==1 & d(1:end-N)==0)/(n-N);   % frequency of 0 then 1
    pBB = sum(d(N+1:end)==0 & d(1:end-N)==0)/(n-N);   % frequency of 0 then 0
    pAB = sum(d(N+1:end)==0 & d(1:end-N)==1)/(n-N);   % frequency of 1 then 0    
    
    out(i,:) = [pA, pAA/(pA*pA),pBA/(pB*pA),pBB/(pB*pB),pAB/(pA*pB)];
    disp(out(i,:));
end

AA = out(:,2);
BB = out(:,4);
AB = out(:,5);
BA = out(:,3);

disp('Averaged values across sessions:')
disp('pA, pAA/(pA*pA),pBA/(pB*pA),pBB/(pB*pB),pAB/(pA*pB)')
m = mean(out); % return mean of each column of out (average the quantities across the sessions)
disp(m);

% for each session, for AA and BB series, calculate their mean along
% each row, and compare (do t-test) with the mean in AB and BA series.
[h,p,ci,stats] = ttest(mean(out(:,[2,4]),2),mean(out(:,[3,5]),2)); % p is p-value

mean1 = mean(mean(out(:,[2,4]),2),1); % mean of AA and BB series
mean2 = mean(mean(out(:,[3,5]),2),1); % AB and BA series
std1 = std(mean(out(:,[2,4]),2)); % standard deviation
std2 = std(mean(out(:,[3,5]),2));

fprintf('average strength of repetition((pAA/(pA*pA)+pBB/(pB*pB))/2) across sessions: %d\n', mean1);
fprintf('standard deviation: %d\n', std1);
fprintf('average strength of alternation((pAB/(pA*pB)+pBA/(pB*pA))/2) across sessions: %d\n', mean2);
fprintf('standard deviation: %d\n', std2);
fprintf('t-test results for comparing whether there is significant difference between the two: t-value %d, p-value: %d\n',stats.tstat,p)
eachComp = [mean1; std1; mean2; std2; stats.tstat; p];

disp(' ')
disp('-------------------------------------------------')
disp('Serial dependency effect in the VFA/VFB variable')
for i=1:length(allData)
    singleSubj = allData{i};
    fprintf('%s:',singleSubj.Session{1});

    response = singleSubj.Response;           % 0: clockwise,  1: ccw
    angularVel = singleSubj.AngularVelocity;     % actual spinning direction (0:cw and 1:ccw)
    azimuth = singleSubj.CameraAzimuth;
    elevation = singleSubj.CameraElevation;             % camera elevation

    % finding perceived elevation
    d = findVFA(angularVel,azimuth,response,elevation,type);
    n = length(d);
    
    pA = sum(d)/n;             % prob of VFA
    pB = 1-pA;                 % prob of viewed-from-below
	
    % divide by the total number of sequences, so you get the prob of sequence of 1,1 or 1,0 etc.
    pAA = sum(d(N+1:end)==1 & d(1:end-N)==1)/(n-N);   % frequency of 1 then 1
    pBA = sum(d(N+1:end)==1 & d(1:end-N)==0)/(n-N);   % frequency of 0 then 1
    pBB = sum(d(N+1:end)==0 & d(1:end-N)==0)/(n-N);   % frequency of 0 then 0
    pAB = sum(d(N+1:end)==0 & d(1:end-N)==1)/(n-N);   % frequency of 1 then 0    
    
	% handle the cases when there are 100% VFA responses, cuz pB will be zero
	% in the denominator, what I did here is make all entries with pB zero. When
	% looking at the data, we need to see that if the entries appear like 1,1,0,0,0
	% then it is likely that there is ceiling in the VFA.
    if pA ~= 1
        out(i,:) = [pA, pAA/(pA*pA),pBA/(pB*pA),pBB/(pB*pB),pAB/(pA*pB)];
    else
        out(i,:) = [pA, pAA/(pA*pA),0,0,0];
    end
    disp(out(i,:));
end

AA = out(:,2);
BB = out(:,4);
AB = out(:,5);
BA = out(:,3);

disp('Averaged values across sessions:')
disp('pA, pAA/(pA*pA),pBA/(pB*pA),pBB/(pB*pB),pAB/(pA*pB)')
m = mean(out); % return mean of each column of out (average the quantities across the sessions)
disp(m);

% for each session, for AA and BB series, calculate their mean along
% each row, and compare (do t-test) with the mean in AB and BA series.
[h,p,ci,stats] = ttest(mean(out(:,[2,4]),2),mean(out(:,[3,5]),2)); % p is p-value

mean1 = mean(mean(out(:,[2,4]),2),1); % mean of AA and BB series
mean2 = mean(mean(out(:,[3,5]),2),1); % AB and BA series
std1 = std(mean(out(:,[2,4]),2)); % standard deviation
std2 = std(mean(out(:,[3,5]),2));

fprintf('average strength of repetition((pAA/(pA*pA)+pBB/(pB*pB))/2) across sessions: %d\n', mean1);
fprintf('standard deviation: %d\n', std1);
fprintf('average strength of alternation((pAB/(pA*pB)+pBA/(pB*pA))/2) across sessions: %d\n', mean2);
fprintf('standard deviation: %d\n', std2);
fprintf('t-test results for comparing whether there is significant difference between the two: t-value %d, p-value: %d\n',stats.tstat,p)
eachComp = [mean1; std1; mean2; std2; stats.tstat; p];

% comparing whether it is more likely to have serial response of AA or BB.
% meanAA = mean(out(:,2))
% meanBB = mean(out(:,4))
% stdAA = std(out(:,2))
% stdBB = std(out(:,4)) 
% disp('difference between prob of having AA and BB');
% [h2,p2,ci2,stats2] = ttest(out(:,2),out(:,4))
