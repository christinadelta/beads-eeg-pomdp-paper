function [set] = RunSlider(scrn, set)

% UNPACK SCREEN RELATED STUFF
window          = scrn.window;      % main window
windrect        = scrn.windrect;
globalrect      = scrn.globalrect;
textsize        = scrn.textsize;
textfont        = scrn.textfont;
grey            = scrn.grey;
white           = scrn.white;
black           = scrn.black;
ycenter         = scrn.ycenter;
ifi             = scrn.ifi;          % frame duration
slack           = scrn.slack;
fixsize         = scrn.fixationsize;

object_offset   = set.object_offset;
fix_dur         = set.fix_dur;      % fixation duration
fixation        = set.fixation;
isi             = set.isi;
jitter          = set.jitter;
EEG             = set.EEG;
taskNb          = set.taskNb;

% UNPACK EEG TRIGGERS
% if EEG == 1
%     responseTrigNotSent     = 1;
%     ioObj                   = set.ioObject;
%     status                  = set.status;
%     triggerdur              = set.triggerdur;
%     address                 = set.address;
%     trgs                    = set.ratetrigger; % confidence rating triggers
% end

% define variables 
anchors         = {'0', '50', '100'};
center          = round([windrect(3) windrect(4)]/2);
line            = 10;
width           = 3;
maxtime         = 10;       % maximum time for dragging the slider
scalelength     = 0.9;      % will change this?
scalepos        = 0.8;      % scale position (0 =top, 1=bottom, and in between)
startpos        = 'left';   % where will be the starting point of the slider?
mousebutton     = 1; 
scalecolour     = white;
slidercolour    = black;

% First define the starting point of the slider
% calculate coordinates of scale line and text bounds
if strcmp(startpos, 'left')
   x = globalrect(3)*(1-scalelength);
elseif strcmp(startpos, 'center')
   x = globalrect(3)/2;
elseif strcmp(startpos, 'right')
   x = globalrect(3)*scalelength;
end

% CREATE WINDOWS FOR FLIPPING
% create fixation cross offscreen and paste later (faster)
fixationdisplay = Screen('OpenOffscreenWindow',window);
Screen('FillRect', fixationdisplay, grey);
Screen('TextFont',fixationdisplay, textfont);
Screen('TextSize',fixationdisplay, fixsize);
DrawFormattedText(fixationdisplay, fixation, 'center', ycenter, white);

% ADD TEXTBOUNDS -- WILL BE USED TO CREATE THE SLIDER 
textbounds          = [Screen('TextBounds', window, sprintf(anchors{1})); Screen('TextBounds', window, sprintf(anchors{3}))];

% calculate coordinates of scale line 
midclick            = [center(1) windrect(4)*scalepos - line - 5 center(1), windrect(4)*scalepos + line + 5];
leftclick           = [windrect(3)*(1-scalelength) windrect(4)*scalepos - line windrect(3)*(1-scalelength) windrect(4)*scalepos + line];
rightclick          = [windrect(3)*(scalelength) windrect(4)*scalepos - line windrect(3)*(scalelength) windrect(4)*scalepos + line];
horzline            = [windrect(3)*scalelength windrect(4)*scalepos windrect(3)*(1-scalelength) windrect(4)*scalepos];

% Calculate the range of the scale, which will be need to calculate the
% position
scalerange          = round(windrect(3)*(1-scalelength)):round(windrect(3)*scalelength); % Calculates the range of the scale (0-100)
scalerangeshifted   = round((scalerange)-mean(scalerange)); % Shift the range of scale so it is symmetrical around zero

% initialise the mouse
SetMouse(round(x), round(windrect(4)*scalepos), window, 1)

t0                         = GetSecs;
respmade                   = 0;

while respmade == 0 
    
    [x,~,buttons,~,~,~] = GetMouse(window, 1);
    
    % Stop at upper and lower bound
    if x > windrect(3)*scalelength
        x = windrect(3)*scalelength;
    elseif x < windrect(3)*(1-scalelength)
        x = windrect(3)*(1-scalelength);
    end

    Screen('TextSize', window, textsize);
    Screen('TextStyle', window, 0)
    Screen('FillRect', window, grey ,windrect);
    DrawFormattedText(window, 'On a scale of 0 to 100, how confident are you for your choice?', 'center', ycenter-150, white);
    DrawFormattedText(window, '0 = Not Confident', 'center', ycenter-50, white);
    DrawFormattedText(window, '100 = Very Confident', 'center', ycenter, white);
    
    % Left, middle and right anchors
    DrawFormattedText(window, anchors{1}, leftclick(1, 1) - textbounds(1, 3)/2,  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Left point
    DrawFormattedText(window, anchors{2}, 'center',  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Middle point
    DrawFormattedText(window, anchors{3}, rightclick(1, 1) - textbounds(2, 3)/2,  windrect(4)*scalepos+40, [],[],[],[],[],[],[]); % Right point

    % Drawing the scale
    Screen('DrawLine', window, scalecolour, midclick(1), midclick(2), midclick(3), midclick(4), width);         % Mid tick
    Screen('DrawLine', window, scalecolour, leftclick(1), leftclick(2), leftclick(3), leftclick(4), width);     % Left tick
    Screen('DrawLine', window, scalecolour, rightclick(1), rightclick(2), rightclick(3), rightclick(4), width); % Right tick
    Screen('DrawLine', window, scalecolour, horzline(1), horzline(2), horzline(3), horzline(4), width);     % Horizontal line
    
    % Draw the slider
    Screen('DrawLine', window, slidercolour, x, windrect(4)*scalepos - line, x, windrect(4)*scalepos  + line, width);
    
    position = round((x)-min(scalerange));                       % Calculates the deviation from 0. 
    position = (position/(max(scalerange)-min(scalerange)))*100; % Converts the value to percentage
    
    DrawFormattedText(window, num2str(round(position)), 'center', windrect(4)*(scalepos - 0.05), white);
    object_onset = Screen('Flip', window, object_offset - slack);    % rating window is on

    % wait for second response 
    secs = GetSecs;
    if buttons(mousebutton) == 1
        respmade = 1;
    end

    rate_rt = (secs - t0);
   
 
end

% release buttons
KbReleaseWait; %Keyboard

object_offset       = object_onset + 0.2 - ifi;

% display fixation
Screen('CopyWindow', fixationdisplay, window, windrect, windrect)
object_onset        = Screen('Flip', window, object_offset - slack); % flip fixation window

object_offset       = object_onset + fix_dur + isi + randperm(jitter*1000,1)/1000 - ifi; 
set.object_offset   = object_offset;
set.rate_rt         = rate_rt;
set.position        = position;

return