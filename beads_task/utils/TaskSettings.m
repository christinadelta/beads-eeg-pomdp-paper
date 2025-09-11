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


if taskNb == 1 % perception task

    % EXPERIMENTAL SETTINGS: Timings
    set.infoscreen  = 2.5;  % this screen appears at the beginning of every sequence and informs the participant abou sequence number and probabilities of draws 
    set.cue_dur     = 0.3;  % cue duration in seconds
    set.prediction  = 1.5;  % prediction time in seconds
    set.delay       = 0.4;  % delay time 
    set.fix_dur     = 0.6;  % duration of the fixation cross
    set.outcome_dur = 0.15; % outcome duration in seconds
    set.feed_dur    = 0.3;  % duration of the feedback window self-paced or up to 3 sec
    
    % TASK PARAMETERS 
    set.blocks      = 4;                        % number of blocks 
    set.trials      = 240;                      % total trials
    set.practice    = 10;                       % Practice trials
    set.numPhases   = 4;                        % probabiltiies switch 4 times in the voaltile condition
    % set.blocktrials = {60,[15,15,15,15]};       % number of trials per block
    set.blocktrials = {60,[15,15,15,15]};       % number of trials per block
    set.nCues       = 2;
    set.combs       = 4;                        % possible cue/outcome combinations
    set.vol         = [1 2];                    % volatility indeces (1 = stable, 2 = volatile)
    set.stoch       = [.9 .1;                   % small stochasticity probabilities
        .8 .2];                                 % large stochasticity probabilities (either 70:30 or 75:3
    set.win         = 2.5;                      % £2.5 pounds if the subject gets >75% correct
    set.balance     = 0;                        % balance starts from zero
    set.threshold   = .75;                      % if accuracy > threshold then subject wilns £2.5 reward
    
    % STIMULUS SETTINGS
    set.stimsize            = 250;              % resize images or not?
    set.stimsize_deg        = 6;                % degrees of visual angle

    if set.EEG == 1

        % EEG TRIGGERS
        set.trigger1        = 1;  % red circle
        set.trigger2        = 2;  % blue circle
        set.trigger3        = 3;  % face
        set.trigger4        = 4;  % house 
        set.trigger5        = 5;  % prediction
        set.trigger14       = 14; % feedback - correct
        set.trigger15       = 15; % feedback incorrect
        
        set.trigger100      = 100; % start of block
        set.trigger101      = 101; % end of block
        set.trigger102      = 102; % stable
        set.trigger103      = 103; % volatile
        set.trigger104      = 104; % small stochasticity
        set.trigger105      = 105; % large stochasticity 
    end

elseif taskNb == 2 % reward task

    % EXPERIMENTAL SETTINGS: Timings
    set.infoscreen      = 2.5;  % this screen appears at the beginning of every sequence and informs the participant abou sequence number and probabilities of draws 
    set.cue_duration    = 1;    % delay time 
    set.fix_dur         = 0.5;  % duration of the fixation cross
    set.feed_dur        = 1.5;    % duration of the feedback window self-paced or up to 3 sec
    set.response        = 1.5;

    % TASK PARAMETERS 
    set.blocks      = 4;                        % number of blocks 
    set.trials      = 240;                      % total trials
    set.practice    = 10;                       % Practice trials
    set. numPhases  = 4;                        % probabiltiies switch 4 times in the voaltile condition
    set.blocktrials = {60,[15,15,15,15]};       % number of trials per block
    set.nCues       = 2;
    set.combs       = 4;                        % possible cue/outcome combinations
    set.vol         = [1 2];                    % volatility indeces (1 = stable, 2 = volatile)
    set.stoch       = [.9 .1;                   % small stochasticity probabilities
        .8 .2];                                 % large stochasticity probabilities (either 70:30 or 75:3
    set.conversion  = 1.1429;                   % to convert point into GBP
    set.balance     = 0;                        % balance starts from zero
    set.threshold   = .75;                      % if accuracy > threshold then subject wilns £2.5 reward

    % define EEG triggers
    if set.EEG == 1

        % start with trial-related triggers
        set.trigger1        = 1; % good option blue, bad option red 
        set.trigger2        = 2; % good option red, bad option blue 
        set.trigger3        = 3; % response onset trigger
        set.trigger4        = 4; % you win! feedback
        set.trigger5        = 5; % you lose.. feedback

        % block related triggers
        set.trigger100      = 100; % start of block
        set.trigger101      = 101; % end of block
        set.trigger102      = 102; % stable
        set.trigger103      = 103; % volatile
        set.trigger104      = 104; % small stochasticity
        set.trigger105      = 105; % large stochasticity 

    end % end of EEG triggers 

elseif taskNb == 3 % beads task

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
end

return

