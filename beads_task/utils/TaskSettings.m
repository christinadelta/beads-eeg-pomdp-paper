function set = TaskSettings(taskNb)

% THIS IS A SUBFUNCTION, PART OF THE "LAB-BASED PEARL PROJECT". 

% it takes as input the task number (from the main script) and outputs a
% structure with a few important startup parameters for each
% task

% NOTE:
% All parameters in this function can change as appropreate 

% GLOBAL PARAMETERS:
set.taskNb      = taskNb;   % initialize settings structure
set.EEG         = 0;        % set to 1 when running in the EEGLAB 
set.triggerdur  = 0.002;    % duration of the trigger in sec (3 ms) 
set.fixation    = '+';      % fixation cross 
set.welcomedur  = 2.5;      % welcome screen duration = 2.5 sec
set.jitter      = .4;       % 0.4 sec
set.isi         = .5;       % in seconds


% EXPERIMENTAL SETTINGS
set.infoscreen      = 2.5;  % this screen appears at the beginning of every sequence and informs the participant abou sequence number and probabilities of draws 
set.bead_dur        = 1;  % bead duration in seconds
set.response        = 2.5;  % self-paced or up to 2.5 sec 
set.fix_dur         = 0.5;  % duration of the fixation cross
set.feed_dur        = 1;    % duration of the feedback window self-paced or up to 3 sec

% TASK PARAMETERS 
set.blocks          = 4;    % number of blocks 
set.trials          = 52;                       % total trials
set.blocktrials     = set.trials/set.blocks;    % number of trials per block
set.draws           = 10;                       % draws per sequence/trial
set.conds           = 2;
set.prob            = [0.8 0.6];                % 8:2 or 6:4 proportion of the beads in each of the two urns
set.penalty         = 0.25;                     % every time subject chooses to draw they are panished with a £0.25 loss
set.win             = 10;                       % £10 pounds if the subject gets the urn right
set.loss            = 10;                       % £10 pounds if the subject gets the urn wrong
set.balance         = 0;                        % balance starts from zero
set.conversion      = 0.3;                   % number of points per pence (will be used for conversion)

% Define triggers
if set.EEG == 1

    % 1. start with the sequence related triggers
    set.trigger1    = 1;             % EASY COND - blue urn
    set.trigger2    = 2;             % EASY COND - green urn
    set.trigger3    = 3;             % DIFF COND - blue urn
    set.trigger4    = 4;             % DIFF COND - green urn
    
    set.trigger5    = 5;             % response prompt
    set.trigger6    = 6;             % prediction screen (probability estimate)
    set.trigger7    = 7;             % response trigger
  
    
    % 2. Confidence rating and feedback relatd triggers
    set.trigger9    = 9;             % confidence screen
    set.trigger14   = 14;            % feedback screen (you win!)
    set.trigger15   = 15;            % feedback screen (you lose!)
    set.trigger16   = 16;            % feedback screen (you lose - out of draws)
    set.trigger19   = 17;            % feedback (didn't respond)
    
    % 3. main script triggers 
    set.trigger100  = 100;           % condition trigger -- condition (easy)
    set.trigger101  = 101;           % condition trigger -- condition (difficult) 
    set.trigger102  = 102;           % sequence start 
    set.trigger103  = 103;           % sequence end

end 

return

