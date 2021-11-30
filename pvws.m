function pvws(par)
% PVWS Phtyo Visc Wetsuit Steadystate
%   PVWS(par) simulation of nutrient flux to a phyto 
%    through a mucus layer. par is class pvw.

if nargin~=1 || ~isa(par,'pvw')
    error(['PVW Object Required:'...
        ' pvw (scale muc b Vf Kf q0 Tcel SA CR CI rad)']);
end

%CONSTANTS
    scale=par.scale;    % Scale cell radius by [um]
    muc=par.muc;        % Mucous curve choice
    b=par.b;            % Mucous curve weighting    
    Vf=par.Vf;          % Viscosity factor increase from mucus, e.g. 1.1
    Kf=par.Kf;          % Thermal cond factor decrease from mucus, e.g. 0.9
    q0=par.q0;          % Phyto heat flux at R, [pW/um^2]
    Tcel=par.Tcel;      % Temperature at infinity, T(inf) [C]
    SA=par.SA;          % Salinity everywhere, S [ppt]
    CR=par.CR;          % Concentration NO3 at cell [umol/L]
    CI=par.CI;          % Concentration NO3 at infinity [umol/L]
    rad=par.rad;        % Radius nutrient molecule, e.g. NO3 [um]

    % DOMAIN
    res=0:0.0001:1; curv=res.^2; %curv=(exp(res)-1)/(exp(1)-1);
    dom=1+99*curv;
    inc=dom(2:end)-dom(1:end-1);
    
    % CONSTANT CALCS
    TI=273.15+Tcel;   %Temperature at infinity, T(inf) [K=273.15+C]
    KI=SW_Conductivity(TI,'K',SA,'ppt'); %Thermal cond. SW infinity at TI, Kt(inf) [W/m*K]
    KR=KI*Kf;      %Thermal cond. cell at TI, Kt(R) [W/m*K]
    VI=SW_Viscosity(TI,'K',SA,'ppt');  %Dyn visc SW infinity at TI, [kg/m*s]
    V0R=VI*Vf;       %Dyn visc cell at TI, [kg/m*s = Pa s]
    kB=1.3806488*(10^-23); %Boltzmann Constant [kg*m^2/K*s^2]

    %SCALE CONSTANTS [um]
    dom=dom*scale;
    inc=inc*scale;
    KR=KR*10^-6; %Kt(R) [W/um*K]
    KI=KI*10^-6; %Kt(inf) [W/um*K]
    q0=q0*10^-12; %q0 [W/um^2]
    VI=VI*10^-6; %V(inf) [kg/um*s]
    V0R=V0R*10^-6; %Vo(R) [kg/um*s]
    VratioScale=10^-6; %Arrhenius calc scale [kg/um*s]
    kB=kB*10^12; %kB [kg*um^2/K*s^2]

%MUCOUS (M)
switch muc
    % Exponential decay
        % a=1 M=0.1 @ 2.3 distance, so a=dist/2.3
        % note that b here is the inverse of the manuscript
    case 'expa'
        M=exp(b*(scale-dom)); % absolute scale
    case 'exps'
        M=exp(b*(scale-dom)/scale); % scaled to cell
    % Exponential peak
    case 'expp'
        l=0; % dist from scale, 0.4
        M=exp(-((dom-scale-l).^2)/b);
    % Sawtooth
    case 'sawt'
        st=0.5*sawtooth(pi:pi/2:6*pi,0.5)+.5;
        xst=linspace(dom(1),dom(1)+1/b,length(st));
        st=[st 0]; xst=[xst dom(end)];
        M=interp1(xst,st,dom);
    % Decaying cosine
    case 'cosd'
        xt=0:0.001:5*pi;
        cod=(0.5*cos(xt)+0.5).*exp(-xt/10);
        xcd=linspace(dom(1),dom(1)+1/b,length(cod));
        cod=[cod 0]; xcd=[xcd dom(end)];
        M=interp1(xcd,cod,dom);
    % Decaying cosine shifted
    case 'cos2'
        xt=0:0.001:6*pi;
        sd=(0.5*cos(xt+pi)+0.5).*exp(-xt/10);
        xsd=linspace(dom(1),dom(1)+1/b,length(sd));
        sd=[sd 0]; xsd=[xsd dom(end)];
        M=interp1(xsd,sd,dom);
end

%THERMAL CONDUCTIVITY (Kt)
Kt=KI-(KI-KR)*M; %Kt curve, asssume dT is small compared to dK.

%TEMPERATURE (T)
TR=TI; %Temperature at R is:
for i=1:length(dom)-1 %Integrate (sum) T(r) from R to end of dom
   TR=TR+(inc(i)*q0*scale^2)/(Kt(i)*dom(i)^2); 
end

T=TR*ones(1,length(dom)); %T curve is:
for i=1:length(dom)-1 %Euler T(r) from R to end of dom
    T(i+1)=T(i)-(inc(i)*q0*scale^2)/(Kt(i)*dom(i)^2);
end

%VISCOSITY (V)
V0=VI+(V0R-VI)*M; %Visc curve at TI
V=V0; %Visc curve for T(r) is:
for i=1:length(dom) %V times ratio of SW visc change w/temp
   Vratio=VratioScale*SW_Viscosity(T(i),'K',SA,'ppt')/VI;
   V(i)=V0(i)*Vratio;
end

%DIFFUSIVITY (D)
D=zeros(1,length(dom)); %Diffusivity [um^2/s] curve is:
for i=1:length(dom)
   D(i)=(kB*T(i))/(6*pi*V(i)*rad);
end

%CONCENTRATION (C)
Cint=0;
for i=1:length(dom)-1 %Integrate (sum) C(r)/c1 from R to end of dom
   Cint=Cint+inc(i)/(D(i)*dom(i)^2); 
end
c1=(CI-CR)/Cint;

C=CR*ones(1,length(dom)); %C curve is:
for i=1:length(dom)-1 %Euler C(r) from R to end of dom
    C(i+1)=C(i)+(inc(i)*c1)/(D(i)*dom(i)^2);
end

%CONSTANT COMPARISONS
Dc=D(end); % Concentration with constant D
CintC=0;
for i=1:length(dom)-1 %Integrate (sum) C(r)/c1 from R to end of dom
   CintC=CintC+inc(i)/(Dc*dom(i)^2); 
end
cc1=(CI-CR)/CintC;

Cc=CR*ones(1,length(dom)); %C curve is:
for i=1:length(dom)-1 %Euler C(r) from R to end of dom
    Cc(i+1)=Cc(i)+(inc(i)*cc1)/(Dc*dom(i)^2);
end

TcR=TI; %Temperature with constant Kt at R is:
for i=1:length(dom)-1 %Integrate (sum) T(r) from R to end of dom
   TcR=TcR+(inc(i)*q0*scale^2)/(KI*dom(i)^2); 
end

Tc=TcR*ones(1,length(dom)); %T curve with constant Kt is:
for i=1:length(dom)-1 %Euler Tc(r) from R to end of dom
    Tc(i+1)=Tc(i)-(inc(i)*q0*scale^2)/(KI*dom(i)^2);
end

Gc=(Cc(2)-Cc(1))/inc(1);
Gm=(C(2)-C(1))/inc(1);
Jc=-Dc*Gc; %Constant D Flux
Jm=-D(1)*Gm; %Variable D Flux

% CHANGE LENGTH UNITS?
if par.lunit=='m'
    uc=10^-6;       % convert um to m
    scale=scale*uc; % um        ->  m
    q0=q0*uc^-2;    % W/um^2    ->  W/m^2
    rad=rad*uc;     % um        ->  m
    Kt=Kt/uc;       % W/um*K    ->  W/m*K
    V=V/uc;         % kg/um*s   ->  kg/m*s
    V0=V0/uc;
    D=D*uc^2;       % um^2/s    ->  m^2/s
    Dc=Dc*uc^2;
    dom=dom*uc;     % um        ->  m
    inc=inc*uc;
    Gc=Gc/uc;       % 1/um      ->  1/m
    Gm=Gm/uc;   
    Jc=Jc*uc;       % um/s      ->  m/s
    Jm=Jm*uc;
end

% ASSIGN DATA TO BASE WORKSPACE VARIABLES
    assignin('base','scale',scale);
    assignin('base','muc',muc);
    assignin('base','b',b);
    assignin('base','Vf',Vf);
    assignin('base','Kf',Kf);
    assignin('base','q0',q0);
    assignin('base','TI',TI);
    assignin('base','SA',SA);
    assignin('base','CR',CR);
    assignin('base','CI',CI);
    assignin('base','rad',rad);
    assignin('base','M',M);
    assignin('base','Kt',Kt);
    assignin('base','T',T);
    assignin('base','V',V);
    assignin('base','V0',V0);
    assignin('base','D',D);
    assignin('base','C',C);
    assignin('base','Tc',Tc);
    assignin('base','Cc',Cc);
    assignin('base','dom',dom);
    assignin('base','inc',inc);
    assignin('base','Dc',Dc);
    assignin('base','Gc',Gc);
    assignin('base','Gm',Gm);
    assignin('base','Jc',Jc);
    assignin('base','Jm',Jm);

end

function k = SW_Conductivity(T,uT,S,uS)
    % SW_Conductivity    Thermal conductivity of seawater
    %=========================================================================
    % USAGE:  k = SW_Conductivity(T,uT,S,uS)
    %
    % DESCRIPTION:
    %   Thermal conductivity of seawater at 0.1 MPa given by [1]
    %   Values at temperature higher than the normal boiling temperature
    %   are calculated at the saturation pressure.
    %
    % INPUT:
    %   T  = temperature
    %   uT = temperature unit
    %        'C'  : [degree Celsius] (ITS-90)
    %        'K'  : [Kelvin]
    %        'F'  : [degree Fahrenheit]
    %        'R'  : [Rankine]
    %   S  = salinity
    %   uS = salinity unit
    %        'ppt': [g/kg] (reference-composition salinity)
    %        'ppm': [mg/kg] (in parts per million)
    %        'w'  : [kg/kg] (mass fraction)
    %        '%'  : [kg/kg] (in parts per hundred)
    %
    %   Note: T and S must have the same dimensions
    %
    % OUTPUT:
    %   k = thermal conductivity [W/m K]
    %
    %   Note: k will have the same dimensions as T and S
    %
    % VALIDITY: 0 < T < 180 C; 0 < S < 160 g/kg
    %
    % ACCURACY: 3.0%
    %
    % REVISION HISTORY:
    %   2009-12-18: Mostafa H. Sharqawy (mhamed@mit.edu), MIT
    %               - Initial version
    %   2012-06-06: Karan H. Mistry (mistry@alum.mit.edu), MIT
    %               - Allow T,S input in various units
    %               - Allow T,S to be matrices of any size
    %
    % DISCLAIMER:
    %   This software is provided "as is" without warranty of any kind.
    %   See the file sw_copy.m for conditions of use and licence.
    %
    % REFERENCES:
    %  [1] D. T. Jamieson, and J. S. Tudhope, Desalination, 8, 393-401, 1970.
    %=========================================================================

    % CHECK INPUT ARGUMENTS

    % CHECK THAT S&T HAVE SAME SHAPE
    if ~isequal(size(S),size(T))
        error('check_stp: S & T must have same dimensions');
    end

    % CONVERT TEMPERATURE INPUT TO Â°C
    switch lower(uT)
        case 'c'
        case 'k'
            T = T - 273.15;
        case 'f'
            T = 5/9*(T-32);
        case 'r'
            T = 5/9*(T-491.67);
        otherwise
            error('Not a recognized temperature unit. Please use ''C'', ''K'', ''F'', or ''R''');
    end

    % CONVERT SALINITY TO PPT
    switch lower(uS)
        case 'ppt'
        case 'ppm'
            S = S/1000;
        case 'w'
            S = S*1000;
        case '%'
            S = S*10;
        otherwise
            error('Not a recognized salinity unit. Please use ''ppt'', ''ppm'', ''w'', or ''%''');
    end

    % CHECK THAT S & T ARE WITHIN THE FUNCTION RANGE
    if ~isequal((T<0)+(T>180),zeros(size(T)))
        warning('Temperature is out of range for thermal conductivity function 0<T<180 C');
    end

    if ~isequal((S<0)+(S>160),zeros(size(S)))
        warning('Salinity is out of range for thermal conductivity function 0<S<160 g/kg');
    end

    % BEGIN

    T = 1.00024*T;      %convert from T_90 to T_68
    S = S / 1.00472;    %convert from S to S_P
    k = 10.^(log10(240+0.0002*S)+0.434*(2.3-(343.5+0.037*S)./(T+273.15)).*(1-(T+273.15)./(647.3+0.03*S)).^(1/3)-3);

end

function mu = SW_Viscosity(T,uT,S,uS)
    % SW_Viscosity    Dynamic viscosity of seawater
    %=========================================================================
    % USAGE:  mu = SW_Viscosity(T,uT,S,uS)
    %
    % DESCRIPTION:
    %   Dynamic viscosity of seawater at atmospheric pressure (0.1 MPa) using
    %   Eq. (22) given in [1] which best fit the data of [2], [3] and [4].
    %   The pure water viscosity equation is a best fit to the data of [5].
    %   Values at temperature higher than the normal boiling temperature
    %   are calculated at the saturation pressure.
    %
    % INPUT:
    %   T  = temperature
    %   uT = temperature unit
    %        'C'  : [degree Celsius] (ITS-90)
    %        'K'  : [Kelvin]
    %        'F'  : [degree Fahrenheit]
    %        'R'  : [Rankine]
    %   S  = salinity
    %   uS = salinity unit
    %        'ppt': [g/kg]  (reference-composition salinity)
    %        'ppm': [mg/kg] (in parts per million)
    %        'w'  : [kg/kg] (mass fraction)
    %        '%'  : [kg/kg] (in parts per hundred)
    %
    %   Note: T and S must have the same dimensions
    %
    % OUTPUT:
    %   mu = dynamic viscosity [kg/m-s]
    %
    %   Note: mu will have the same dimensions as T and S
    %
    % VALIDITY: 0 < T < 180 C and 0 < S < 150 g/kg;
    %
    % ACCURACY: 1.5%
    %
    % REVISION HISTORY:
    %   2009-12-18: Mostafa H. Sharqawy (mhamed@mit.edu), MIT
    %               - Initial version
    %   2012-06-06: Karan H. Mistry (mistry@alum.mit.edu), MIT
    %               - Allow T,S input in various units
    %               - Allow T,S to be matrices of any size
    %
    % DISCLAIMER:
    %   This software is provided "as is" without warranty of any kind.
    %   See the file sw_copy.m for conditions of use and licence.
    %
    % REFERENCES:
    %   [1] M. H. Sharqawy, J. H. Lienhard V, and S. M. Zubair, Desalination
    %       and Water Treatment, 16, 354-380, 2010. (http://web.mit.edu/seawater/)
    %   [2] B. M. Fabuss, A. Korosi, and D. F. Othmer, J., Chem. Eng. Data 14(2), 192, 1969.
    %   [3] J. D. Isdale, C. M. Spence, and J. S. Tudhope, Desalination, 10(4), 319 - 328, 1972
    %   [4] F. J. Millero, The Sea, Vol. 5, 3 – 80, John Wiley, New York, 1974
    %   [5] IAPWS release on the viscosity of ordinary water substance 2008
    %=========================================================================

    % CHECK INPUT ARGUMENTS

    % CHECK THAT S&T HAVE SAME SHAPE
    if ~isequal(size(S),size(T))
        error('check_stp: S & T must have same dimensions');
    end

    % CONVERT TEMPERATURE INPUT TO °C
    switch lower(uT)
        case 'c'
        case 'k'
            T = T - 273.15;
        case 'f'
            T = 5/9*(T-32);
        case 'r'
            T = 5/9*(T-491.67);
        otherwise
            error('Not a recognized temperature unit. Please use ''C'', ''K'', ''F'', or ''R''');
    end

    % CONVERT SALINITY TO PPT
    switch lower(uS)
        case 'ppt'
        case 'ppm'
            S = S/1000;
        case 'w'
            S = S*1000;
        case '%'
            S = S*10;
        otherwise
            error('Not a recognized salinity unit. Please use ''ppt'', ''ppm'', ''w'', or ''%''');
    end

    % CHECK THAT S & T ARE WITHIN THE FUNCTION RANGE
    if ~isequal((T<0)+(T>180),zeros(size(T)))
        warning('Temperature is out of range for Viscosity function 0<T<180 C');
    end

    if ~isequal((S<0)+(S>150),zeros(size(S)))
        warning('Salinity is out of range for Viscosity function 0<S<150 g/kg');
    end

    % BEGIN

    S = S/1000;

    a = [
        1.5700386464E-01
        6.4992620050E+01
       -9.1296496657E+01
        4.2844324477E-05
        1.5409136040E+00
        1.9981117208E-02
       -9.5203865864E-05
        7.9739318223E+00
       -7.5614568881E-02
        4.7237011074E-04
    ];

    mu_w = a(4) + 1./(a(1)*(T+a(2)).^2+a(3));


    A  = a(5) + a(6) * T + a(7) * T.^2;
    B  = a(8) + a(9) * T + a(10)* T.^2;
    mu = mu_w.*(1 + A.*S + B.*S.^2);

end
