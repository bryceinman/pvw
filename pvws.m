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
