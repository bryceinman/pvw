classdef pvw
    % Phytoplankton Viscosity Wetsuit Class
    properties
        scale=1;        % Scale cell radius by [um]
        muc='expa';     % Mucous curve choice
        b=1;            % Mucous curve weighting    
        Vf=2;           % Viscosity factor increase from mucous, e.g. 1.1
        Kf=0.975;       % Thermal cond factor decrease from mucous, e.g. 0.9
        q0=0.03;        % Phyto heat flux at R, [pW/um^2 = W/m^2]
        Tcel=20;        % Temperature at infinity, T(inf) [C]
        SA=35;          % Salinity everywhere, S [ppt]
        CR=0;           % Concentration NO3 at cell [umol/L]
        CI=1;           % Concentration NO3 at infinity [umol/L]
        rad=0.000179;   % Radius nutrient molecule, e.g. NO3 [um]
        lunit='um';     % Output length unit [um or m]
        
        % Time Dependent
        rn='00';        % Run Number
        ncm=[1 0 0];    % Parameters: nutrient pulse, constant conc, muc pulse
        dr=1e-3;        % Space increment
        dt=1e-6;        % Time step: a=10 dt=dr^2, a=1 dt=.01*dr^2
        td=1;           % Total run time [s]
        ts=0.1;         % Record data on this time increment [2]
        CP=100;         % Max pulse concentration at outer domain
        pt=0.1;        % Time at which pulse max occurs
        pw=1e-4;        % Duration of pulse (width of exp. curve)
    end
end