% Main file for running static and dynamic convex optimizations on motor 
% and planetary gear                                                     

clear all; close all; clc;
dataFile;           %Load data       
load('roosload.mat');   %Load roos load data

selection1 = input('Do you want to perform both a static and dynamic optimization? (y/n)', 's');

if(selection1 == 'y') %Dynamic optimization on motor and gearbox
    DynamicMG;
else                  %Static optimization on motor, gearbox or both
    selection2 = input('Do you want to optimize for motor, planetary gear or both? (m/g/b)','s');
    if(selection2 == 'm')
        Motor;
    elseif(selection2 == 'g')
        PlanetaryGear;    
    else
        StaticMG;
    end
end