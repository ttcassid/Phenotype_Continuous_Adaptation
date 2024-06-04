% Matlab 2020b Code to reproduce the Pressley et al. example in the pre-print:
% Reducing phenotype-structured PDE models of cancer evolution to systems
% of ODEs: a generalised moment dynamics approach (arXiv:2406.01505)
% T Cassidy, June 2024.
close all
clear all
%% Load model parameters from Pressley et al. (Table 2)
PA.rMax = 0.45 ; % Intrinsic growth rate
PA.kMax = 1e4; % Carrying capacity
PA.k = 2; % innate cell response
PA.b = 10; % Resistance benefit
PA.delta = 0.7*PA.kMax;  % Threshold for disease progression  (Not used here).
PA.g = 0.1;  % Cost of resistance in growth rate
PA.d = 0.01;   % Intrinsic death rate
PA.m = 1; % Drug effect

FigureSwitch = 2; % 1 to reproduce figures A & B, 2 to reproduce figures C & D.
if FigureSwitch == 1
    PA.alpha = 0.0025; % Speed of adaptation
    tf = 65;
elseif FigureSwitch == 2
    PA.alpha = 0.0001; % Speed of adaptation
    tf = 500;
end
%% Simulation time
t0 = 0;
TotalTime = [t0 tf];
%% Initial condtions
PIC = 6e3;
XIC = 0;
SigmaIC = 0;
%% Treatment options 
PA.AdaptiveTherapySwitch = 1; % turn on adaptive therapy vs MTD
NTest = 1; % A variable to pass the current treatment state 

%Adaptive therapy requires switch from MTD to no dosing  if tumour size is less than 0.5*PIC;
PA.AdaptiveTherapySwitchConcOff = 0.5.*PIC; % Lower bound for treatment off
PA.AdaptiveTherapySwitchConcOn = 1.0.*PIC;  % Upper bound for treatment on
PA.TreatmentIndicator = 1;  % Indicator if treatment is applied or not

%% Solve the ODE systems
ICAdaptive = [PIC,XIC,0];
%% Solve the Adaptive Dynamics ODE system for the Pressley et al. model
PA.AdaptiveCase = 1;
[solMean] =  PressleyAdaptiveModel(ICAdaptive,TotalTime,PA);

%% Solve the Adaptive Dynamics ODE with heterogeneity system for beta = 2 alpha
ICVariance = [PIC,XIC,SigmaIC,0];
PA.beta = PA.alpha/0.5;
PA.AdaptiveCase = 1;
[solVariance1] =  PressleyAdaptiveWithVarianceModel(ICVariance,TotalTime,PA);

% Solve the Adaptive Dynamics ODE with heterogeneity system for beta = alpha
PA.beta = PA.alpha/1;
PA.AdaptiveCase = 1;
[solVariance2] =  PressleyAdaptiveWithVarianceModel(ICVariance,TotalTime,PA);

% Solve the Adaptive Dynamics ODE with heterogeneity system for beta = alpha/2
PA.beta = PA.alpha/2;
PA.AdaptiveCase = 1;
[solVariance3] =  PressleyAdaptiveWithVarianceModel(ICVariance,TotalTime,PA);

% Solve the Adaptive Dynamics ODE with heterogeneity system for beta = alpha/10
PA.beta = PA.alpha/10;
PA.AdaptiveCase = 1;
[solVariance4] =  PressleyAdaptiveWithVarianceModel(ICVariance,TotalTime,PA);

%% Relative change figures
NColor(1,:) = [224,130,20]/256;
NColor(2,:) = [84,39,136]/256;
NColor(3,:) = [227,26,28]/256;
NColor(4,:) = [27,120,55]./256;  

% TransparencyValue = 0;
CurveSize = 1.25;
VarianceSize = 0.75;

Fig = figure(1);
g3 =  plot(solMean.x(:), solMean.y(1,:),'LineWidth',CurveSize,'Color', [67,147,195]/256,'LineStyle','-'); %grey
hold on
g4 =  plot(solVariance1.x(:), solVariance1.y(1,:),'LineWidth',CurveSize,'Color', NColor(1,:),'LineStyle','-');
hold on
g4 =  plot(solVariance2.x(:), solVariance2.y(1,:),'LineWidth',CurveSize,'Color', NColor(2,:),'LineStyle','-');
hold on
g4 =  plot(solVariance3.x(:), solVariance3.y(1,:),'LineWidth',CurveSize,'Color', NColor(3,:),'LineStyle','-');
hold on
g4 =  plot(solVariance4.x(:), solVariance4.y(1,:),'LineWidth',CurveSize,'Color', NColor(4,:),'LineStyle','-');
hold on
ax = gca(Fig);
ax.FontSize = 12; 
ylim([0 1e4]);
yticks([0 2e3 4e3 6e3 8e3 10e3])
if FigureSwitch == 1
    xlim([-1 min(tf,300)])
    xticks([0 20 40 60])
elseif FigureSwitch == 2
    xlim([175 425])
    xticks([0 200 400])
end

%% figure for the population phenotype 
CurveSize = 1.0;
Fig = figure(2);
g3 =  plot(solMean.x(:), solMean.y(2,:),'LineWidth',CurveSize,'Color', [67,147,195]/256,'LineStyle','-'); %grey
hold on

g4 =  plot(solVariance1.x(:), solVariance1.y(2,:),'LineWidth',CurveSize,'Color',NColor(1,:),'LineStyle','-'); %grey
hold on
g5 =  plot(solVariance1.x(:), solVariance1.y(2,:)-solVariance1.y(3,:) ,'LineWidth',VarianceSize,'Color',NColor(1,:),'LineStyle','--'); %grey
hold on
g6 =  plot(solVariance1.x(:), solVariance1.y(2,:)+solVariance1.y(3,:) ,'LineWidth',VarianceSize,'Color', NColor(1,:),'LineStyle','--'); %grey
hold on

g4 =  plot(solVariance2.x(:), solVariance2.y(2,:),'LineWidth',CurveSize,'Color',NColor(2,:),'LineStyle','-'); %grey
hold on
g5 =  plot(solVariance2.x(:), solVariance2.y(2,:)-solVariance2.y(3,:) ,'LineWidth',VarianceSize,'Color',NColor(2,:),'LineStyle','--'); %grey
hold on
g6 =  plot(solVariance2.x(:), solVariance2.y(2,:)+solVariance2.y(3,:) ,'LineWidth',VarianceSize,'Color', NColor(2,:),'LineStyle','--'); %grey
hold on

g4 =  plot(solVariance3.x(:), solVariance3.y(2,:),'LineWidth',CurveSize,'Color',NColor(3,:),'LineStyle','-'); %grey
hold on
g5 =  plot(solVariance3.x(:), solVariance3.y(2,:)-solVariance3.y(3,:) ,'LineWidth',VarianceSize,'Color',NColor(3,:),'LineStyle','--'); %grey
hold on
g6 =  plot(solVariance3.x(:), solVariance3.y(2,:)+solVariance3.y(3,:) ,'LineWidth',VarianceSize,'Color', NColor(3,:),'LineStyle','--'); %grey
hold on

g4 =  plot(solVariance4.x(:), solVariance4.y(2,:),'LineWidth',CurveSize,'Color',NColor(4,:),'LineStyle','-'); %grey
hold on
g5 =  plot(solVariance4.x(:), solVariance4.y(2,:)-solVariance4.y(3,:) ,'LineWidth',VarianceSize,'Color',NColor(4,:),'LineStyle','--'); %grey
hold on
g6 =  plot(solVariance4.x(:), solVariance4.y(2,:)+solVariance4.y(3,:) ,'LineWidth',VarianceSize,'Color', NColor(4,:),'LineStyle','--'); %grey
hold on
ax = gca(Fig);
ax.FontSize = 12; 

if FigureSwitch == 1
    xlim([-1 min(tf,300)])
    xticks([0 20 40 60])
    ylim([-0.03 1.25])
    yticks([0 0.5 1])
elseif FigureSwitch == 2
    yticks([0 0.2 0.4 0.6 0.8])
    xlim([-1 min(tf,700)])
    xticks([0 200 400])
end

%% Figure to see the treatment administered; 
% Fig = figure(3);  
%  g3 =  plot(solMean.x(:), solMean.y(end,:),'LineWidth',CurveSize,'Color', [67,147,195]/256,'LineStyle','-'); %grey
%  hold on
%  g4 =  plot(solVariance1.x(:), solVariance1.y(end,:),'LineWidth',CurveSize,'Color',  NColor(1,:),'LineStyle','-'); %grey
%  hold on
%   g4 =  plot(solVariance2.x(:), solVariance2.y(end,:),'LineWidth',CurveSize,'Color',  NColor(2,:),'LineStyle','-'); %grey
%  hold on
%   g4 =  plot(solVariance3.x(:), solVariance3.y(end,:),'LineWidth',CurveSize,'Color',  NColor(3,:),'LineStyle','-'); %grey
%  hold on
%   g4 =  plot(solVariance4.x(:), solVariance4.y(end,:),'LineWidth',CurveSize,'Color',  NColor(4,:),'LineStyle','-'); %grey
%  hold on



function [Output,NTest] = TreatmentTwoOutput(PA,P);  

if PA.AdaptiveTherapySwitch == 1 % If we are simulating adaptive therapy treatment.

    if P > PA.AdaptiveTherapySwitchConcOn -1e-6 % If tumour size is above upper threshold
        OutputVec =  PA.m; % Output treatment
        PA.TreatmentIndicator = 1; % Indicate that treatment is being administered 
        NTest = 1; % PA.AdaptiveCase;
    elseif P < PA.AdaptiveTherapySwitchConcOff % If tumour size is below lower threshold
        OutputVec = 0; % Output no treatment
        PA.TreatmentIndicator = 0; % Indicate that treatment is not being administered 
        NTest = 0; % PA.AdaptiveCase;
    else % If not above upper threshold, or below lower threshold, continue: applying treatment if treatment is being applied, or not, if not. 
        NTest = PA.AdaptiveCase ;  % NTest+1
        PA.TreatmentIndicator = PA.AdaptiveCase ;
    end

    Output = PA.TreatmentIndicator.*PA.m;
else % if not adaptive therapy , output treatment
    Output = PA.m;
    NTest = 0;
end

end

function [AdaptiveCase] = CaseCheck(PA,P)
[~, NTest] = TreatmentTwoOutput(PA,P);
AdaptiveCase = NTest;
end


function [Output] = Treatment(PA,P)
[OutputTest, ~] = TreatmentTwoOutput(PA,P);
Output = OutputTest;
end

function [Output] = GrowthRate(PA,x) % cost of resistance
Output = PA.rMax.*exp(-PA.g.*x);
end

function [Output] = FitnessFunction(PA,P,x) % The G function
Output = GrowthRate(PA,x).*(1-P/PA.kMax) - PA.d - Treatment(PA,P).*(1/(PA.k+PA.b.*x) ) ;
end

function [Output] = AdaptationFunction(PA,P,x) % First derivative of the fitness function G (also the drift function)
Output = (-PA.g).*GrowthRate(PA,x).*(1-P/PA.kMax)  + PA.b.* Treatment(PA,P).*(1/(PA.k+PA.b.*x)).^2;
end

function [Output] = SecondGDerivative(PA,P,x) % Second derivative of the fitness function G
Output =   ( (-PA.g).^2.*GrowthRate(PA,x).*(1-P/PA.kMax)  - 2.*PA.b.^2.* Treatment(PA,P).*(1/(PA.k+PA.b.*x)).^3 );
end

function [Output] = ThirdGDerivative(PA,P,x) % Third derivative of the fitness function G
Output =   (  (-PA.g).^3.*GrowthRate(PA,x).*(1-P/PA.kMax)  + 6.*PA.b.^3.* Treatment(PA,P).*(1/(PA.k+PA.b.*x)).^4 );
end

function [sol] = PressleyAdaptiveModel(totaltime,IC,PA) %ODE model with therapy
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
sol = ode45(@PressleyAdaptiveDynamics,IC,totaltime,opts);
    function dydt = PressleyAdaptiveDynamics(t,y);
        PA.AdaptiveCase =  CaseCheck(PA,y(1)); %Check if treatment is being applied at the current time step. To be passed to the treatment function.
        dydt(1) = FitnessFunction(PA,y(1),y(2))*y(1) ; % ODE for Tumour size
        dydt(2) = PA.alpha*AdaptationFunction(PA,y(1),y(2)); % ODE for mean phenotype
        dydt(3) = Treatment(PA,y(1)); % keep track of total treatment
        dydt = dydt';
    end
end

function [sol] = PressleyAdaptiveWithVarianceModel(totaltime,IC,PA) %ODE model with therapy
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',1e-2);
sol = ode45(@PressleyAdaptiveDynamics,IC,totaltime,opts);
    function dydt = PressleyAdaptiveDynamics(t,y);
        PA.AdaptiveCase =  CaseCheck(PA,y(1));
        dydt(1) = (FitnessFunction(PA,y(1),y(2)) + SecondGDerivative(PA,y(1),y(2))*y(3)/2 )*y(1) ; % ODE for Tumour size
        dydt(2) = PA.alpha*AdaptationFunction(PA,y(1),y(2)) + ( PA.alpha*ThirdGDerivative(PA,y(1),y(2))/2 + AdaptationFunction(PA,y(1),y(2)) -3*y(2)*SecondGDerivative(PA,y(1),y(2))/2 )*y(3);  % ODE for mean phenotype
        dydt(3) = 2*PA.beta  + 2.*PA.alpha*AdaptationFunction(PA,y(1),y(2)).*y(2)  ...
            + 2*( PA.alpha*SecondGDerivative(PA,y(1),y(2)) - ( 2*PA.alpha*ThirdGDerivative(PA,y(1),y(2))/2 + AdaptationFunction(PA,y(1),y(2))/2 )*y(2) + SecondGDerivative(PA,y(1),y(2))/2*(y(2)^2+y(3)/2) ).*y(3); % ODE for variance 
        dydt(4) = Treatment(PA,y(1)); % keep track of total treatment
        dydt = dydt';
    end
end
 

