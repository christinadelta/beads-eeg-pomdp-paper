function [set] = ShortQuiz(set, scrn, keys)

% this function runs a short quiz related to the general instructions of
% the tasks.

% The input is the task number and the output is the resposes 
% NOTE:
% subject must respond to all the questions correctly in order to complete
% the quiz. If the subject response for a question is incorrect then the
% subject will try again 

% BEADS QUIZ QUESTIONS:
% 1. What are the possible majority colours of the urns?
% 2. On a given trial/sequence, from how many urns can you draw? 
% 3. If the colour split in a given urn is 60:40, what does that mean
% exactly? 
% 4. What is the key-press for choosing the BLUE urn?
% 5. What is the key-press for choosing the GREEN urn?
% 6. What is the key-press for choosing to draw again?
% 7. Up to how many times can you draw before making a decision?
% 8. What is the range of the confidence ratings?
% 9. With how many mouse clicks can you complete a confidence rating?


% POSSIBLE QUESTIONS
% 9. From how many urns can you draw?
% 

% The questions will displayed be randomised 

% FOR SIMPLICITY:
% For odd number questions (1,3,5,7,9), the correct response number has index 2
% For even number questions (2,4,6,8),the correct response number has
% index 4

% UNPACK THE SETTINGS AND SCREEN STRUCTS
window          = scrn.window;       % main window
windrect        = scrn.windrect;
xcenter         = scrn.xcenter;
ycenter         = scrn.ycenter;
ifi             = scrn.ifi;          % frame duration
slack           = scrn.slack;        % slack is ifi/2 (important for timing)
white           = scrn.white;
grey            = scrn.grey;
green           = scrn.green;
red             = scrn.red;
fixsize         = scrn.fixationsize;
textfont        = scrn.textfont;
textsize        = scrn.textsize;

isi             = set.isi;          % interstimulus interval
jitter          = set.jitter; 
fixation        = set.fixation;     % draw fixation
resptime        = 10;               % 2 sec to respond. after that the question will re-appear
fix_dur         = set.fix_dur;      % duration of the fixation

choiceA          = keys.code7;  % subject pressed A
choiceB          = keys.code8;  % subject pressed B
choiceC          = keys.code9;  % subject pressed C
choiceD          = keys.code10; % subject pressed D

% create fixation cross offscreen and paste later (faster)
fixationdisplay = Screen('OpenOffscreenWindow',window);
Screen('FillRect', fixationdisplay, grey);
Screen('TextFont',fixationdisplay, textfont);
Screen('TextSize',fixationdisplay, fixsize);
DrawFormattedText(fixationdisplay, fixation, 'center', ycenter, white);

% fixation -- green (for correct repsonses)
fixation_green = Screen('OpenOffscreenWindow',window);
Screen('FillRect', fixation_green, grey);
Screen('TextFont',fixation_green, textfont);
Screen('TextSize',fixation_green, fixsize);
DrawFormattedText(fixation_green, fixation, 'center', ycenter, green);

% fixation -- red (for correct repsonses)
fixation_red = Screen('OpenOffscreenWindow',window);
Screen('FillRect', fixation_red, grey);
Screen('TextFont',fixation_red, textfont);
Screen('TextSize',fixation_red, fixsize);
DrawFormattedText(fixation_red, fixation, 'center', ycenter, red);

% CREATE WINDOWS FOR FLIPPING 
% Create a background window 
% background_window = Screen('OpenOffscreenWindow',window);
% Screen('TextSize', background_window, textsize);
% Screen('TextFont',background_window, textfont);
% Screen('FillRect', background_window, grey ,windrect);

qs          = 9; % How many questions? 
answ        = 4; % how many potential answers per question?
questions   = 1:qs; % [1:8] list with questions indecies
answers     = 1:answ; % [1:4] list with answer indecies

rand_qs     = questions(randperm(qs)); % randomise the questions


% INIT STRINGS STRUCT. This will contain all the questions and potential answers
strings = [];

% Field 1: Questions  
strings.questions = {'what are the possible majority colours of the urns?', 'On a given trial/sequence, from how many urns can you draw balls?',...
    'If the colour split within a given urn is 60:40, what does that mean exactly?', 'What is the key for choosing the blue urn?',...
    'What is the key for choosing the green urn?', 'What is the key for choosing to darw again?',...
    'Up to how many times can you draw before making a decision?', 'What is the range of the confidence ratings?',...
    'With how many mouse clicks can you complete a confidence rating?'};

% Field 2: Answers 
strings.answers = {'red-green', 'blue-green', 'green-yellow', 'red-blue',
    '2', '4', '3', '1',
    'this is an "easy condition" trial', 'this is a "hard condition" trial', 'nothing, there is no 60:40 colour split in the task', 'The majority of balls are of blue colour',
    '2', '5', '3', '1', 
    '1', '2', '6', '3',
    '2', '6', '4', '3',
    '6', '9', '5', '10',
    '1-10', '0-10', '1-100', '0-100',
    '1', '2', '0', '3'};

% Field 3: Alphabetical order 
strings.order = {'A) ', 'B) ', 'C) ', 'D) '};

Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
fliptime    = Screen('Flip', window); % flip fixation window
trialstart  = fliptime;

% object offset
objectoff   = trialstart + isi + randperm(jitter*1000,1)/1000 - ifi;

i = 1; % counter for while loop

while i <= qs 
    
    % RANDOMISE THE ANSWER INDECIES ON EVERY ITERATION
    rand_answ       = answers(randperm(answ)); 
    thisquestion    = rand_qs(i);
   
    % DISPLAY THE CURRENT QUESTION WITH THE 4 POTENTIAL REPONSES 
    background_window = Screen('OpenOffscreenWindow', window, windrect);
    Screen('TextSize', background_window, textsize);
    Screen('FillRect', background_window, grey ,windrect);
    DrawFormattedText(background_window, strings.questions{thisquestion}, 'center', ycenter-200, white);
    DrawFormattedText(background_window, [strings.order{1}, strings.answers{thisquestion, rand_answ(1)}], 'center', ycenter-100, white); 
    DrawFormattedText(background_window, [strings.order{2}, strings.answers{thisquestion, rand_answ(2)}], 'center', ycenter-50, white);
    DrawFormattedText(background_window, [strings.order{3}, strings.answers{thisquestion, rand_answ(3)}], 'center', ycenter, white);
    DrawFormattedText(background_window, [strings.order{4}, strings.answers{thisquestion, rand_answ(4)}], 'center', ycenter+50, white);
    
    Screen('CopyWindow',background_window, window, windrect, windrect);
    objecton        = Screen('Flip', window, objectoff - slack); 
    startquiz       = objecton;
    
    % INIT RESPONSE 
    responded       = NaN;
    resp_input      = 0;
    
    while resp_input == 0 && (GetSecs - startquiz) < resptime - 2*slack
        
        [~, secs, keycode] = KbCheck; % check for input 
        
        if keycode(1, choiceA)
            resp_input  = choiceA;
            rt          = secs - startquiz; 
            responded   = 1; % answered A
            responset   = secs;
            
        elseif keycode(1, choiceB)
            resp_input  = choiceB;
            rt          = secs - startquiz; 
            responded   = 2; % answered A
            responset   = secs;
            
        elseif keycode(1, choiceC)
            resp_input  = choiceC;
            rt          = secs - startquiz; 
            responded   = 3; % answered A
            responset   = secs;
            
        elseif keycode(1, choiceD)
            resp_input  = choiceD;
            rt          = secs - startquiz; 
            responded   = 4; % answered A
            responset   = secs;
            
        else
            resp_input = 0; 
            rt          = nan; 
            responded   = nan; % answered A
            responset   = GetSecs;
            
            
        end % end of response if 
        
    end % end of response while loop
    
    objectoff   = responset + isi - ifi;                                  % question-response window self paced or on for 2500 ms
    
    % First see if this was an even or odd question
    if mod(thisquestion,2) == 0 % if this is an even question
        pos = find(rand_answ == 4);
    else % if this is an odd question
        pos = find(rand_answ == 2);
    end

    % check if response was correct
    correct = pos == responded;
    
    % if response was correct show green fixation and move to the next
    % question
    if correct
        
        % show green fixation!
        Screen('CopyWindow', fixation_green,window, windrect, windrect) 
        objecton = Screen('Flip', window, objectoff - slack); 
        
        objectoff = objecton + fix_dur - ifi; 
        
        i = i + 1; % move to the next question
        
        set.quizresp{i} = responded;
        set.quizrt{i}   = rt;
        
    else % if response is incorrect 
        
        % show red fixation!
        Screen('CopyWindow', fixation_red,window, windrect, windrect) 
        objecton = Screen('Flip', window, objectoff - slack); 
        
        objectoff = objecton + fix_dur - ifi; 
        
        i = i; % stay in the same question until subject gets it right 
        % i = i + 1; % for testing. while loops are weird

    end % end of correct statement 
    
    % bring white fixation back on
    Screen('CopyWindow', fixationdisplay,window, windrect, windrect)
    objecton = Screen('Flip', window, objectoff - slack);
    
    objectoff = objecton + isi + randperm(jitter*1000,1)/1000 - ifi;
 
end % end of while loop

% show "end of the questions window
Screen('OpenOffscreenWindow', window, windrect);
Screen('TextSize', window, scrn.textsize);
Screen('FillRect', window, scrn.grey ,windrect);
DrawFormattedText(window, 'Well Done! You got all the questions right.', 'center', ycenter-50, white);
DrawFormattedText(window, 'You may now start the experiment', 'center', ycenter, white);
Screen('Flip',window);
WaitSecs(2);
    
    

end