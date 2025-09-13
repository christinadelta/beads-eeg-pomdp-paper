function scrn = screenSettings(scrn,taskNb)

% screen settings of the PeARL project
% 1. defines screen settings
% 2. defines object x and y 

%% define screen parameters 

% change these parameters as appropreate 
% screen resolution 
scrn.actscreenRes   = scrn.actscreen;   % get screen's actual resolution
scrn.screenRes      = [1280 800];       % this also the windrect in px
scrn.hz             = 60; 
scrn.distview       = 700;
scrn.width          = scrn.actwidth;
scrn.height         = scrn.actheight;



end