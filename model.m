%% Solution 2 (Laplace)
clear
% Select rate constant inputs   "k1,k2,k3,k4,k5,k6,kin" and set time vector
pa = [ 0.81, 0.85, 0.65, 0.70, 0.43, 0.27, 0.16];  tt = 0:1:20;
% Select Drug Input function
func_hydrogel = 4.2138*log(tt)+1.0061;
% Concentration of drug in vitreous over time is ROI
[ROI] = Function(pa,tt,func_hydrogel)
