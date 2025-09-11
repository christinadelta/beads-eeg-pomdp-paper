function set = DefineKeys(taskNb, set)

% This subfunction is part of the PeARL ptojrvt. 
% It runs via the main script of each task using the task number (taskNb)

% it creates a list of keys for each experiment alongside with the "global" 
% keys used in all experiments

% GLOBAL KEYS
% 1. Escape key (allows subject to quit the experiment)
% 2. Space key (allows subject to start the trials after instructions)
KbName('UnifyKeyNames');
set.code20          = KbName('space');
set.code21          = KbName('ESCAPE');

if taskNb == 1 || taskNb == 2
    
    % Task keys
    KbName('UnifyKeyNames');
    set.code1      = KbName('f'); % choose face
    set.code2      = KbName('h'); % choose house

elseif taskNb == 3 % beads task
    
    % TASK 1 KEYS
    KbName('UnifyKeyNames');
    set.code1       = KbName('1!'); % 1! = Blue
    set.code2       = KbName('2@'); % 2@ = Green
    set.code3       = KbName('3#'); % 3# = Draw again

    set.code7       = KbName('a'); % answer a
    set.code8       = KbName('b'); % answer b
    set.code9       = KbName('c'); % answer c
    set.code10      = KbName('d'); % answer d
    

end % end of if statement

end 