clear;%clc;%close all

%==========================================================================
% This script contains the user inputs. It creates structures containing
% the volumetric parameters (vol) and the respiratory settings (res) and
% makes the call to Solver.m.
% Lastly, it saves the workspace.
%==========================================================================

%% MANUAL USER INPUTS

% choose morphometry data
data_tog = 1; % 1 = Yeh, 2 = Weibel
%--------------------------------------------------------------------------

% choose constant or varying compliance
comp_tog = 1; % 0 = constant, 1 = varying
%--------------------------------------------------------------------------

% choose to simulate with or without turbulent effects
turb_tog = 0; % 0 = turbulence off, 1 = turbulence on
%--------------------------------------------------------------------------

% choose deposition scenario
mech = 4; % 0 = 'no deposition', 1='impaction only', 2='sedimentation only', 3 = 'diffusion only', 4 = all, 5 = 'no diffusion';
%--------------------------------------------------------------------------

% choose constant flow rates or dynamic
output = 0; %0 = dynamic; 1 = constant
%--------------------------------------------------------------------------

% choose respiratory settings

% Compare to Rissler (2017)
% TLC_liters = 6.7; % average TLC reported by Rissler (2017)
% TV_liters = 0.75;
% RV_liters = 2;
% RR = 10; % enter in breaths/minute
% IE = 1; 
% FRC_liters = 3.4;
% res.Pmusmin = -7.51;% for TV = 0.75 % cmH2O %*pcf; (this parameter changes the amplitude of Pmus)

% Compare to Heyder 
TLC_liters = 6; % average TLC reported by Rissler (2017)
% TV_liters = 1;
TV_liters = 2;
% TV_liters = 0.5;
% TV_liters = 1.5;
RV_liters = 1.86;
% RR = 15; % enter in breaths/minute
RR = 3.75; % enter in breaths/minute
IE = 1; 
FRC_liters = 3;
% res.Pmusmin = -10.65;% TV = 1 for varying compliance
% res.Pmusmin = -8.85;% TV = 1 for constant compliance
res.Pmusmin = -33.7; % TV = 2 for varying compliance
% res.Pmusmin = -17.3; % TV = 2 for constant compliance
% res.Pmusmin = -4.65;% TV = 0.5 for varying compliance
% res.Pmusmin = -4.3;% TV = 0.5 for constant compliance
% res.Pmusmin = -20.2;% TV = 1.5 for varying compliance
% res.Pmusmin = -12.95;% TV = 1.5 for constant compliance

% Compare to Darquenne 
% TLC_liters = 5.6; % average TLC reported by Rissler (2017)
% TV_liters = 0.9;
% RV_liters = 1.86;
% RR = 15; % enter in breaths/minute
% IE = 1; 
% FRC_liters = 3;
% res.Pmusmin = -9.4;% for TV = 0.9 % cmH2O %*pcf; (this parameter changes the amplitude of Pmus)

%% AUTOMATIC SET UP
% Make necessary conversions
% Calculate dependent variables 
% Store volumes in structure called 'vol' 
% Store respiratory parameters in structure called 'res'
% Call the solver

vol.TV = TV_liters;%/1e3; % convert to cubic meters
vol.TLC = TLC_liters;%/1e3; % m^3 - total lung capacity ~ 6 liters in healthy adult male
vol.FRC = FRC_liters;%/1e3;
vol.RV = RV_liters; %/1e3

res.RR = RR; % breaths/min
res.IE = IE; % dimensionless (shifts the transition between inspiration and expiration in Pmus)
res.T = 60/res.RR; % respiratory period in seconds
res.TE = res.T/(1+res.IE);
res.TI = res.TE*res.IE;
res.tau = res.TE/5; 

% MORPHOMETRY

if data_tog==1
    sheet_name1 = 'Yeh1980';
    display('Using Yeh data')
end

if data_tog==2
    sheet_name1 = 'Weibel1963';
    sheet_name2 = 'Yeh1980';
    display('Using Weibel data')
end

Y = readtable('Morphometry','Sheet',sheet_name1,'VariableNamingRule','preserve');

% MAKE CALL TO THE SOLVER--------------------------------------------------
if output==1
    [f,g] = ConstantFlow_Solver(vol,res,Y);
else
    [f,g,t,s] = Solver(vol,res,Y,mech,comp_tog,turb_tog);
    if comp_tog==0
        if turb_tog==0
            save(['Constant_',num2str(res.RR),'_',num2str(vol.TV),'.mat'])
        elseif turb_tog==1
            save(['Constant_',num2str(res.RR),'_',num2str(vol.TV),'_turbulence','.mat'])
        end
    elseif comp_tog==1
        if turb_tog==0
            save(['Varying_',num2str(res.RR),'_',num2str(vol.TV),'.mat'])
        elseif turb_tog==1
            save(['Varying_',num2str(res.RR),'_',num2str(vol.TV),'_turbulence','.mat'])
        end
    end
end
