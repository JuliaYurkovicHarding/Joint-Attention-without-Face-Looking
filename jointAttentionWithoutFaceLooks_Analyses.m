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
%   See also: jointAttentionWithoutFaceLooks_SupplementalAnalyses.m

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

Naming = readtable([data_dir,'Naming.csv']);
% File with group name and number, participant number, percent of looks
% with joint attention where the parent named the object of the child's
% gaze, and the percent of looks without joint attention where the parent
% named the object of the child's gaze.

%% -- SET ANALYSIS VARIABLES ----------------------------------------------

% Groups are coded as TD = 1 and ASD = 2
groups = {'TD','ASD'};

% There are 18 TD participants and 19 ASD participants
Participants{1} = transpose(101:118);
Participants{2} = transpose(201:219);

%% -- RESULTS 1a: CHILD GAZE TO FACES -------------------------------------

disp('-----------------------------------------')
disp('Child Look to Parent Face')
disp(' ')

% Preallocate a cell array for the number of groups
rateChildFaceLook = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    rateChildFaceLook{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the amount of time that child looked to parent face
        %   Offset = timestamp of offset of child look to face (seconds)
        %   Onset = timestamp of onset of child look to face (seconds)
        %   Offset - Onset = duration in seconds of each look to face
        %   Sum of (Offset - Onset) = total time looking to face
        amountChildFaceLook = ...
            sum(ChildFaceLook.Offset(ChildFaceLook.Participant==Participants{g}(px))...
            - ChildFaceLook.Onset(ChildFaceLook.Participant==Participants{g}(px)));
        
        % Determine the total amount of usable data
        amountUsableData = ...
            LengthOfUsableData.UsableDataDuration...
            (LengthOfUsableData.Participant == Participants{g}(px));
        
        % Determine percent of usable data that child looked to parent face
        rateChildFaceLook{g}(px) = ...
            (amountChildFaceLook/amountUsableData) * 100;
        
        % Clear temporary variables associated with current participant
        clear amountChildFaceLook amountUsableData
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f%% (%4.1f%%)\n\n',...
        mean(rateChildFaceLook{g},'omitnan'),...
        std(rateChildFaceLook{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(rateChildFaceLook{1},rateChildFaceLook{2});
d = cohensd(rateChildFaceLook{1},rateChildFaceLook{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

%% -- Results 1b: Parent Gaze to Faces ------------------------------------

disp('-----------------------------------------')
disp('Parent Look to Child Face')
disp(' ')

% Preallocate a cell array for the number of groups
rateParentFaceLook = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    rateParentFaceLook{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the amount of time that parent looked to child face
        %   Offset = timestamp of offset of parent look to face (seconds)
        %   Onset = timestamp of onset of parent look to face (seconds)
        %   Offset - Onset = duration in seconds of each look to face
        %   Sum of (Offset - Onset) = total time looking to face
        amountParentFaceLook = ...
            sum(ParentFaceLook.Offset(ParentFaceLook.Participant==Participants{g}(px))...
            - ParentFaceLook.Onset(ParentFaceLook.Participant==Participants{g}(px)));
        
        % Determine the total amount of usable data
        amountUsableData = ...
            LengthOfUsableData.UsableDataDuration...
            (LengthOfUsableData.Participant == Participants{g}(px));
        
        % Determine percent of usable data that parent looked to child face
        rateParentFaceLook{g}(px) = ...
            (amountParentFaceLook/amountUsableData) * 100;
        
        % Clear temporary variables associated with current participant
        clear amountParentFaceLook amountUsableData
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f%% (%4.1f%%)\n\n',...
        mean(rateParentFaceLook{g},'omitnan'),...
        std(rateParentFaceLook{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(rateParentFaceLook{1},rateParentFaceLook{2});
d = cohensd(rateParentFaceLook{1},rateParentFaceLook{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

%% -- Results 1c: Child vs. Parent Gaze to Faces --------------------------

disp('-----------------------------------------')
disp('Child vs. Parent Look Faces')
disp(' ')

% Loop through groups
for g = 1:length(groups)
    
    % Run and display t-test and Cohen's d analyses comparing the rate of
    % child and parent looks to faces for each dyad
    disp(groups{g})
    [~,p,~,stats] = ttest(rateChildFaceLook{g},rateParentFaceLook{g});
    d = cohensd(rateChildFaceLook{g},rateParentFaceLook{g});
    fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
    clear p stats d
    
end
clear g

%% -- Results 1d: Mutual Face Looking -------------------------------------

disp('-----------------------------------------')
disp('Mutual Face Looking')
disp(' ')

% Preallocate a cell array for the number of groups
rateMutualFaceLook = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    rateMutualFaceLook{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the amount of time that dyads engaged in mutual face looking
        %   Offset = timestamp of offset of mutual face look (seconds)
        %   Onset = timestamp of onset of mutual face look (seconds)
        %   Offset - Onset = duration in seconds of each mutual face look
        %   Sum of (Offset - Onset) = total time mutually face looking
        amountMutualFaceLook = ...
            sum(MutualFaceLook.Offset(MutualFaceLook.Participant==Participants{g}(px))...
            - MutualFaceLook.Onset(MutualFaceLook.Participant==Participants{g}(px)));
        
        % Determine the total amount of usable data
        amountUsableData = ...
            LengthOfUsableData.UsableDataDuration...
            (LengthOfUsableData.Participant == Participants{g}(px));
        
        % Determine percent of usable data that dyad engaged in mutual face look
        rateMutualFaceLook{g}(px) = ...
            (amountMutualFaceLook/amountUsableData) * 100;
        
        % Clear temporary variables associated with current participant
        clear amountMutualFaceLook amountUsableData
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f%% (%4.1f%%)\n\n',...
        mean(rateMutualFaceLook{g},'omitnan'),...
        std(rateMutualFaceLook{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(rateMutualFaceLook{1},rateMutualFaceLook{2});
d = cohensd(rateMutualFaceLook{1},rateMutualFaceLook{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

%% -- Results 2a: Joint Attention -----------------------------------------

disp('-----------------------------------------')
disp('Joint Attention')
disp(' ')

% Preallocate a cell array for the number of groups
rateJointAttention = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    rateJointAttention{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Determine the amount of time that dyads engaged in joint attention
        %   Offset = timestamp of offset of joint attention (seconds)
        %   Onset = timestamp of onset of joint attention (seconds)
        %   Offset - Onset = duration in seconds of each joint attention
        %   Sum of (Offset - Onset) = total time in joint attention
        amountJointAttention = ...
            sum(JointAttention.Offset(JointAttention.Participant==Participants{g}(px))...
            - JointAttention.Onset(JointAttention.Participant==Participants{g}(px)));
        
        % Determine the total amount of usable data
        amountUsableData = ...
            LengthOfUsableData.UsableDataDuration...
            (LengthOfUsableData.Participant == Participants{g}(px));
        
        % Determine percent of usable data that dyad engaged in joint attention
        rateJointAttention{g}(px) = ...
            (amountJointAttention/amountUsableData) * 100;
        
        % Clear temporary variables associated with current participant
        clear amountJointAttention amountUsableData
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f%% (%4.1f%%)\n\n',...
        mean(rateJointAttention{g},'omitnan'),...
        std(rateJointAttention{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(rateJointAttention{1},rateJointAttention{2});
d = cohensd(rateJointAttention{1},rateJointAttention{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

%% -- Results 2b: Naming during Joint Attention ---------------------------

disp('-----------------------------------------')
% Naming during Looks without Joint Attention
disp('Naming during Looks without Joint Attention: Between Group Comparison')
disp(' ')

% Preallocate a cell array for the number of groups
percentNamingNonJA = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    percentNamingNonJA{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Pull the percent of looks with joint attention where parents
        % named objects
        percentNamingNonJA{g}(px,1) = ...
            Naming.PercentNamingDuringNonJA(Naming.Participant==Participants{g}(px));
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f%% (%4.1f%%)\n\n',...
        mean(percentNamingNonJA{g},'omitnan'),...
        std(percentNamingNonJA{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(percentNamingNonJA{1},percentNamingNonJA{2});
d = cohensd(percentNamingNonJA{1},percentNamingNonJA{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

% Naming during Looks with Joint Attention
disp('Naming during Looks with Joint Attention: Between Group Comparison')
disp(' ')

% Preallocate a cell array for the number of groups
percentNamingJA = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Preallocate a matrix for the number of participants in each group
    percentNamingJA{g} = zeros(size(Participants{g}));
    
    % Loop through participants
    for px = 1:length(Participants{g})
        
        % Pull the percent of looks with joint attention where parents
        % named objects
        percentNamingJA{g}(px,1) = ...
            Naming.PercentNamingDuringJA(Naming.Participant==Participants{g}(px));
        
    end
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f%% (%4.1f%%)\n\n',...
        mean(percentNamingJA{g},'omitnan'),...
        std(percentNamingJA{g},'omitnan'));
    
end
clear g p

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(percentNamingJA{1},percentNamingJA{2});
d = cohensd(percentNamingJA{1},percentNamingJA{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

% Naming during Looks with Joint Attention
disp('Naming during JA and non-JA Looks: Within Group Comparison')
disp(' ')

% Loop through groups
for g = 1:length(groups)
    
    % Run and display t-test and Cohen's d analyses
    disp(groups{g})
    [~,p,~,stats] = ttest(percentNamingNonJA{g},percentNamingJA{g});
    d = cohensd(percentNamingNonJA{g},percentNamingJA{2});
    fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
    clear p stats d
    
end
clear g

% Naming during Looks with Joint Attention
disp('Increase in Likelihood of Naming from JA to non-JA')
disp(' ')

% Divide rate of naming during JA by rate of naming during non-JA to get
% the relative increase in the likelihood of naming
Naming.LikelihoodComparison = ...
    Naming.PercentNamingDuringJA./Naming.PercentNamingDuringNonJA;

% Change Inf values to NaN values for analysis purposes
Naming.LikelihoodComparison(Naming.LikelihoodComparison == Inf) = NaN;

% Preallocate a cell array for the number of groups
namingLikelihoodComparison = cell(length(groups),1);

% Loop through groups
for g = 1:length(groups)
    
    % Pull likelihood comparison for each group
    namingLikelihoodComparison{g} = ...
        Naming.LikelihoodComparison(Naming.GroupNum==g);
    
    % Calculate and display group means and standard deviations
    disp(groups{g})
    fprintf('mean (SD)\n%4.1f%% (%4.1f%%)\n\n',...
        mean(namingLikelihoodComparison{g},'omitnan'),...
        std(namingLikelihoodComparison{g},'omitnan'));
    
end
clear g

% Run and display t-test and Cohen's d analyses
[~,p,~,stats] = ttest2(namingLikelihoodComparison{1},namingLikelihoodComparison{2});
d = cohensd(namingLikelihoodComparison{1},namingLikelihoodComparison{2});
fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
clear p stats d

%% -- Results 3a: Cues for Parents to Follow into Joint Attention ---------

disp('-----------------------------------------')
disp('Cues for Parents to Follow into Joint Attention')
disp(' ')

% Fieldnames for the two cues that parents can follow into joint attention
cues = fieldnames(JointAttention);
cues = cues(contains(cues,'Child'));

% Preallocate a structure
rateCue = struct;

% Loop through cues
for c = 1:length(cues)
    
    disp(cues{c})
    disp(' ')
    
    % Preallocate a cell array for the number of groups
    rateCue.(cues{c}) = cell(length(groups),1);
    
    % Loop through groups
    for g = 1:length(groups)
        
        % Preallocate a matrix for the number of participants in each group
        rateCue.(cues{c}){g} = zeros(size(Participants{g}));
        
        % Loop through participants
        for px = 1:length(Participants{g})
            
            % Determine the sum of joint attention events with the
            % preceding behavioral cue
            jointAttentionWithCue = ...
                sum(JointAttention.(cues{c})...
                (JointAttention.Participant==Participants{g}(px)));
            
            % Determine the total number of joint attention events
            totalJointAttention = sum(JointAttention.Participant==Participants{g}(px));
            
            % Percent of joint attention events with preceding cue
            rateCue.(cues{c}){g}(px) = ...
                (jointAttentionWithCue/totalJointAttention) * 100;
            
            % Clear temporary variables associated with current participant
            clear jointAttentionWithCue totalJointAttention
            
        end
        
        % Calculate and display group means and standard deviations
        disp(groups{g})
        fprintf('mean (SD)\n%4.1f%% (%4.1f%%)\n\n',...
            mean(rateCue.(cues{c}){g},'omitnan'),...
            std(rateCue.(cues{c}){g},'omitnan'));
        
    end
    
    % Run and display t-test and Cohen's d analyses
    [~,p,~,stats] = ttest2(rateCue.(cues{c}){1},...
        rateCue.(cues{c}){2});
    d = cohensd(rateCue.(cues{c}){1},rateCue.(cues{c}){2});
    fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
    clear p stats d
    
end
clear c g px

%% -- Results 3b: Cues for Children to Follow into Joint Attention --------

disp('-----------------------------------------')
disp('Cues for Children to Follow into Joint Attention')
disp(' ')

% Fieldnames for the two cues that parents can follow into joint attention
cues = fieldnames(JointAttention);
cues = cues(contains(cues,'Parent'));

% Loop through cues
for c = 1:length(cues)
    
    disp(cues{c})
    disp(' ')
    
    % Preallocate a cell array for the number of groups
    rateCue.(cues{c}) = cell(length(groups),1);
    
    % Loop through groups
    for g = 1:length(groups)
        
        % Preallocate a matrix for the number of participants in each group
        rateCue.(cues{c}){g} = zeros(size(Participants{g}));
        
        % Loop through participants
        for px = 1:length(Participants{g})
            
            % Determine the sum of joint attention events with the
            % preceding behavioral cue
            jointAttentionWithCue = ...
                sum(JointAttention.(cues{c})...
                (JointAttention.Participant==Participants{g}(px)));
            
            % Determine the total number of joint attention events
            totalJointAttention = sum(JointAttention.Participant==Participants{g}(px));
            
            % Percent of joint attention events with preceding cue
            rateCue.(cues{c}){g}(px) = ...
                (jointAttentionWithCue/totalJointAttention) * 100;
            
            % Clear temporary variables associated with current participant
            clear jointAttentionWithCue totalJointAttention
            
        end
        
        % Calculate and display group means and standard deviations
        disp(groups{g})
        fprintf('mean (SD)\n%4.1f%% (%4.1f%%)\n\n',...
            mean(rateCue.(cues{c}){g},'omitnan'),...
            std(rateCue.(cues{c}){g},'omitnan'));
        
    end
    
    % Run and display t-test and Cohen's d analyses
    [~,p,~,stats] = ttest2(rateCue.(cues{c}){1},...
        rateCue.(cues{c}){2});
    d = cohensd(rateCue.(cues{c}){1},rateCue.(cues{c}){2});
    fprintf('t(%d)=%4.2f, p=%4.2f, d=%4.2f\n\n',stats.df,stats.tstat,p,d);
    clear p stats d
    
end
clear c g px

%% -- Cohen's d Function --------------------------------------------------

function [d] = cohensd(data1,data2)

d = ...
    abs((mean(data1,'omitnan') - mean(data2,'omitnan'))/...
    std([data1(:);data2(:)],'omitnan'));

end