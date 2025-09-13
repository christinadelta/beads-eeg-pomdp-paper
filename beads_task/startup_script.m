% STARTUP SCRIPT 

% CLEAN UP        
clear;
clc  
close all hidden;

% restoredefaultpath
  
%%

% THE SCRIPT SHOULD RUN IN DIRECTORY : '/information_sampling_eeg_beads_paper/beads_task'
% PRESS RUN OR TYPE "startup_script" IN THE COMMAND WINDOW AND ENTER 
% THE TASK CODE NAME. THIS WILL ADD THE CORRECT PATHS OF THE CURRENT
% EXPERIMENT TO YOUR MATLAB PATH, THEN IT WILL RUN THE CORRECT "main_task"
% SCRIPT.

% TASK CODE NAMES:
   
% beads                                 = beads task

% WORKING DIRECTORIES AND INFO:


% working_dir                           = pwd

% core experimental functions           = /working_dir/utils
% task directory                        = /working_dir/tasks
% participant log files                 = /working_dir/results


% DEFINE INITIAL PATHS
startpath           = pwd;
wokingpath          = fullfile(startpath, 'tasks');

% create a user input dialog to gather information
prompt          = {'Enter task name (e.g. beads):','Enter subject number (e.g. 01:'};
dlgtitle        = 'Info window';
dims            = [1 30];
definput        = {'beads','01'}; % this is a default input (this should change)
answer          = inputdlg(prompt,dlgtitle,dims,definput);

startup.answer  = answer;    % participant number and task name
getpath         = answer{1}; % usefull to read the correct task name and run the main script
        

%% 

taskpath        = wokingpath;
addpath(taskpath);
main_beads % run the main beads script

clear startpath workingpath

%%  