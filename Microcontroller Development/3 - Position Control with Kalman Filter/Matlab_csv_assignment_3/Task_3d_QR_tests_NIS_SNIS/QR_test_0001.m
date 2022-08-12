clear all
close all
clear global
clc

R = 3*10^-7;
obj = KalmanExperiment.createfromQRC3();
%[probNIS, probSNIS] = analyzeconsistency(obj, 0.95, 5);
analyzeconsistency(obj, 0.95, 5)