% STARTUP SCRIPT 

% CLEAN UP        
clear;
clc  
close all hidden;

% restoredefaultpath
  
%%

% THE SCRIPT SHOULD RUN IN DIRECTORY : '/lab_study'
% PRESS RUN OR TYPE "startup_script" IN THE COMMAND WINDOW AND ENTER 
% THE TASK CODE NAME. THIS WILL ADD THE CORRECT PATHS OF THE CURRENT
% EXPERIMENT TO YOUR MATLAB PATH, THEN IT WILL RUN THE CORRECT "main_task"
% SCRIPT.

% TASK CODE NAMES:
   
% perception        = perceptual task 
% reward            = reward learning task

% WORKING DIRECTORIES AND INFO:


% working_dir                           = /lab_study/

% core experimental functions           = /working_dir/utils
% experimental stimuli                  = /working_dir/stimuli
% tasks directory                       = /working_dir/tasks
% participant log files                 = /working_dir/results


% DEFINE INITIAL PATHS
startpath           = pwd;
wokingpath          = fullfile(startpath, 'tasks');
% taskpath            = fullfile(workingpath, 'tasks');

% create a user input dialog to gather information
prompt          = {'Enter task name (e.g. perception):','Enter subject number (e.g. 01:'};
dlgtitle        = 'Info window';
dims            = [1 30];
definput        = {'perception','01'}; % this is a default input (this should change)
answer          = inputdlg(prompt,dlgtitle,dims,definput);

startup.answer  = answer;    % participant number and task name
getpath         = answer{1}; % usefull to read the correct task name and run the main script
        

%% 

switch getpath
    
    case 'perception'
        taskpath = wokingpath;
        addpath(taskpath);
        main_perception % run main perception script  
    case 'reward'     
        taskpath = wokingpath;
        addpath(taskpath);
        main_reward % run main reward script  
    case 'beads'
        taskpath = wokingpath;
        addpath(taskpath);
        main_beads % run the main beads script
end

clear startpath workingpath

%%    