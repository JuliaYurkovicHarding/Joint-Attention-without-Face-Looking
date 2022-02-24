% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%            CHILDREN WITH ASD ESTABLISH JOINT ATTENTION DURING            
%                 FREE-FLOWING TOY PLAY WITHOUT FACE LOOKS                 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Article Authors:
%   Julia Yurkovic-Harding, Grace Lisandrelli, Rebecca C. Shaffer,
%   Kelli C. Dominick, Ernest V. Pedapati, Craig A. Erickson,
%   Chen Yu, and Daniel P. Kennedy

% Link to Data on OSF:
%   {LINK}

% Accompanying Analyses
%   Script Last Updated: January 2022
%   Script Author: Julia Yurkovic-Harding
%   Matlab Version 2021a
%   See also: jointAttentionWithoutFaceLooks_Analyses.m

% All results are reported in the command window

%% -- PREPARE WORKSPACE ---------------------------------------------------
clc
close all
clear
rng(1)

%% -- LOAD DATA TO WORKSPACE ----------------------------------------------

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Set path to the location of the downloaded .csv files on your computer
% Make sure to add the final folder delineator / at the end of your path
data_dir = '/Users/juliayurkovicharding/Desktop/Joint Attention Data/';
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Load each data file:

LengthOfUsableData = readtable([data_dir,'LengthOfUsableData.csv']);
% File with group name and number, participant number, and the amount of
% usable data for each participant in seconds.

ChildFaceLook = readtable([data_dir,'ChildFaceLook.csv']);
% File with group name and number, participant number, onset and offset of
% each child look to parent face in seconds.

ParentFaceLook = readtable([data_dir,'ParentFaceLook.csv']);
% File with group name and number, participant number, onset and offset of
% each parent look to child face in seconds.

MutualFaceLook = readtable([data_dir,'MutualFaceLook.csv']);
% File with group name and number, participant number, onset and offset of
% each mutual face look in seconds.

JointAttention = readtable([data_dir,'JointAttention.csv']);
% File with group name and number, participant number, onset and offset of
% each joint attention instance in seconds, the ROI (or toy number) of each
% joint attention instance. Ones and zeros in each of the "Preceding"
% columns indicates if the variable of interest was present before joint
% attention or not.

%% -- SET ANALYSIS VARIABLES ----------------------------------------------

% Groups are coded as TD = 1 and ASD = 2
groups = {'TD','ASD'};

% There are 18 TD participants and 19 ASD participants
Participants{1} = transpose(101:118);
Participants{2} = transpose(201:219);

%% -- RESULTS S3a: COUNT AND DURATION OF CHILD LOOK TO PARENT FACE --------

disp('-----------------------------------------')
disp('Child Looks to Parent Face per Minute')
disp(' ')

% Preallocate a cell array for the number of groups
perMinuteChildFaceLook = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    perMinuteChildFaceLook{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the number of times that parent looked to child face
        countChildFaceLook = ...
            sum(ChildFaceLook.Participant==Participants{g}(px));
        
        % Determine the total amount of usable data
        amountUsableData = ...
            LengthOfUsableData.UsableDataDuration...
            (LengthOfUsableData.Participant == Participants{g}(px));
        
        % Convert data in seconds to data in minutes
        amountUsableData = amountUsableData/60;
        
        % Determine the number of parent looks to child face per minute
        perMinuteChildFaceLook{g}(px) = ...
            (countChildFaceLook/amountUsableData);
        
        % Clear temporary variables associated with current participant
        clear countChildFaceLook amountUsableData
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f (%4.1f)\n\n',...
        mean(perMinuteChildFaceLook{g},'omitnan'),...
        std(perMinuteChildFaceLook{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ...
    ttest2(perMinuteChildFaceLook{1},perMinuteChildFaceLook{2});
d = cohensd(perMinuteChildFaceLook{1},perMinuteChildFaceLook{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d


disp('-----------------------------------------')
disp('Median Duration of Child Looks to Parent Face')
disp(' ')

% Preallocate a cell array for the number of groups
durationChildFaceLook = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    durationChildFaceLook{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the duration of parent looks to child face
        %   Offset = timestamp of offset of child look to face (seconds)
        %   Onset = timestamp of onset of child look to face (seconds)
        %   Offset - Onset = duration in seconds of each look to face
        %   Median of (Offset - Onset) = median duration of looks to face
        durationChildFaceLook{g}(px) = ...
            median(ChildFaceLook.Offset(ChildFaceLook.Participant==Participants{g}(px))...
            - ChildFaceLook.Onset(ChildFaceLook.Participant==Participants{g}(px)),'omitnan');
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1fs (%4.1fs)\n\n',...
        mean(durationChildFaceLook{g},'omitnan'),...
        std(durationChildFaceLook{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ...
    ttest2(durationChildFaceLook{1},durationChildFaceLook{2});
d = cohensd(durationChildFaceLook{1},durationChildFaceLook{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

%% -- RESULTS S3b: COUNT AND DURATION OF PARENT LOOK TO CHILD FACE ---------

disp('-----------------------------------------')
disp('Parent Looks to Child Face per Minute')
disp(' ')

% Preallocate a cell array for the number of groups
perMinuteParentFaceLook = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    perMinuteParentFaceLook{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the number of times that parent looked to child face
        countParentFaceLook = ...
            sum(ParentFaceLook.Participant==Participants{g}(px));
        
        % Determine the total amount of usable data
        amountUsableData = ...
            LengthOfUsableData.UsableDataDuration...
            (LengthOfUsableData.Participant == Participants{g}(px));
        
        % Convert data in seconds to data in minutes
        amountUsableData = amountUsableData/60;
        
        % Determine the number of parent looks to child face per minute
        perMinuteParentFaceLook{g}(px) = ...
            (countParentFaceLook/amountUsableData);
        
        % Clear temporary variables associated with current participant
        clear countParentFaceLook amountUsableData
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f (%4.1f)\n\n',...
        mean(perMinuteParentFaceLook{g},'omitnan'),...
        std(perMinuteParentFaceLook{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(perMinuteParentFaceLook{1},perMinuteParentFaceLook{2});
d = cohensd(perMinuteParentFaceLook{1},perMinuteParentFaceLook{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d


disp('-----------------------------------------')
disp('Median Duration of Parent Looks to Child Face')
disp(' ')

% Preallocate a cell array for the number of groups
durationParentFaceLook = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    durationParentFaceLook{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the duration of parent looks to child face
        %   Offset = timestamp of offset of child look to face (seconds)
        %   Onset = timestamp of onset of child look to face (seconds)
        %   Offset - Onset = duration in seconds of each look to face
        %   Median of (Offset - Onset) = median duration of looks to face
        durationParentFaceLook{g}(px) = ...
            median(ParentFaceLook.Offset(ParentFaceLook.Participant==Participants{g}(px))...
            - ParentFaceLook.Onset(ParentFaceLook.Participant==Participants{g}(px)),'omitnan');
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1fs (%4.1fs)\n\n',...
        mean(durationParentFaceLook{g},'omitnan'),...
        std(durationParentFaceLook{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ...
    ttest2(durationParentFaceLook{1},durationParentFaceLook{2});
d = cohensd(durationParentFaceLook{1},durationParentFaceLook{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

%% -- RESULTS S3c: COUNT AND DURATION OF MUTUAL FACE LOOK -----------------

disp('-----------------------------------------')
disp('Mutual Face Looks per Minute')
disp(' ')

% Preallocate a cell array for the number of groups
perMinuteMutualFaceLook = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    perMinuteMutualFaceLook{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the number of times that dyad mutually looked to faces
        countMutualFaceLook = ...
            sum(MutualFaceLook.Participant==Participants{g}(px));
        
        % Determine the total amount of usable data
        amountUsableData = ...
            LengthOfUsableData.UsableDataDuration...
            (LengthOfUsableData.Participant == Participants{g}(px));
        
        % Convert data in seconds to data in minutes
        amountUsableData = amountUsableData/60;
        
        % Determine the number of mutual face looks per minute
        perMinuteMutualFaceLook{g}(px) = ...
            (countMutualFaceLook/amountUsableData);
        
        % Clear temporary variables associated with current participant
        clear countMutualFaceLook amountUsableData
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f (%4.1f)\n\n',...
        mean(perMinuteMutualFaceLook{g},'omitnan'),...
        std(perMinuteMutualFaceLook{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(perMinuteMutualFaceLook{1},perMinuteMutualFaceLook{2});
d = cohensd(perMinuteMutualFaceLook{1},perMinuteMutualFaceLook{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d


disp('-----------------------------------------')
disp('Median Duration of Mutual Face Looks')
disp(' ')

% Preallocate a cell array for the number of groups
durationMutualFaceLook = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    durationMutualFaceLook{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the duration of mutual face looks
        %   Offset = timestamp of offset of mutual face look (seconds)
        %   Onset = timestamp of onset of mutual face look (seconds)
        %   Offset - Onset = duration in seconds of each mutual face look
        %   Median of (Offset - Onset) = median duration of mutual face looks
        durationMutualFaceLook{g}(px) = ...
            median(MutualFaceLook.Offset(MutualFaceLook.Participant==Participants{g}(px))...
            - MutualFaceLook.Onset(MutualFaceLook.Participant==Participants{g}(px)),'omitnan');
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1fs (%4.1fs)\n\n',...
        mean(durationMutualFaceLook{g},'omitnan'),...
        std(durationMutualFaceLook{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ...
    ttest2(durationMutualFaceLook{1},durationMutualFaceLook{2});
d = cohensd(durationMutualFaceLook{1},durationMutualFaceLook{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

%% -- RESULTS S3d: COUNT AND DURATION OF JOINT ATTENTION ------------------

disp('-----------------------------------------')
disp('Joint Attention per Minute')
disp(' ')

% Preallocate a cell array for the number of groups
perMinuteJointAttention = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    perMinuteJointAttention{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine number of times that the dyad was in joint attention
        countJointAttention = ...
            sum(JointAttention.Participant==Participants{g}(px));
        
        % Determine the total amount of usable data
        amountUsableData = ...
            LengthOfUsableData.UsableDataDuration...
            (LengthOfUsableData.Participant == Participants{g}(px));
        
        % Convert data in seconds to data in minutes
        amountUsableData = amountUsableData/60;
        
        % Determine the number of joint attention instances per minute
        perMinuteJointAttention{g}(px) = ...
            (countJointAttention/amountUsableData);
        
        % Clear temporary variables associated with current participant
        clear countJointAttention amountUsableData
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f (%4.1f)\n\n',...
        mean(perMinuteJointAttention{g},'omitnan'),...
        std(perMinuteJointAttention{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(perMinuteJointAttention{1},perMinuteJointAttention{2});
d = cohensd(perMinuteJointAttention{1},perMinuteJointAttention{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d


disp('-----------------------------------------')
disp('Median Duration of Joint Attention')
disp(' ')

% Preallocate a cell array for the number of groups
durationJointAttention = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    durationJointAttention{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the duration of joint attention
        %   Offset = timestamp of offset of joint attention (seconds)
        %   Onset = timestamp of onset of joint attention (seconds)
        %   Offset - Onset = duration in seconds of each joint attention
        %   Median of (Offset - Onset) = median duration of joint attention
        durationJointAttention{g}(px) = ...
            median(JointAttention.Offset(JointAttention.Participant==Participants{g}(px))...
            - JointAttention.Onset(JointAttention.Participant==Participants{g}(px)),'omitnan');
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1fs (%4.1fs)\n\n',...
        mean(durationJointAttention{g},'omitnan'),...
        std(durationJointAttention{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(durationJointAttention{1},durationJointAttention{2});
d = cohensd(durationJointAttention{1},durationJointAttention{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

%% -- COHEN'S D FUNCTION --------------------------------------------------

function [d] = cohensd(data1,data2)

d = ...
    abs((mean(data1,'omitnan') - mean(data2,'omitnan'))/...
    std([data1(:);data2(:)],'omitnan'));

end