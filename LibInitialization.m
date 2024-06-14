function LibInitialization()
    %clear all; close all; clc;
    % Add MPH folder to path
    mphFolder = 'Src';  % Replace '/path/to/MPH' with the actual path to your MPH folder
    addpath(genpath(mphFolder)); % genpath adds folders and all subfolders
    
    % Add Src folder to path
    srcFolder = 'MPH';  % Replace '/path/to/Src' with the actual path to your Src folder
    addpath(genpath(srcFolder)); % genpath adds folders and all subfolders
end
