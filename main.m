%% This code accompanies the manuscript entitled "A compartmental framework
%% for transitioning from patient-level to population-scale epidemiological
%% dynamics" by Hart et al. For further information about the paper or this
%% code, please email william.hart@keble.ox.ac.uk

%% We request that users cite the original publication when referring to
%% this code or any results generated from it.

%% This code reproduces the panels in Figure 2A-B of our paper.

clear all; close all; clc;


%% Patient-level dynamics

% Load patient-level data (the vector V_mean_vector contains the mean viral
% load, calculated over 10,000 within-host realisations, at the times since
% infection contained in x_vector).

load('Data/patient_level_data.mat','x_vector','V_mean_vector')


% Calculate expected infectiousness, beta_vector, at times since infection
% x_vector, assuming infectiousness is proportional to viral load.

R0 = 1.5; %Basic reproduction number
N = 1000; %Population size
beta_vector = R0*V_mean_vector/N;


% Plot expected infectiousness curve in Figure 1

figure(1); hold on;
plot(x_vector,beta_vector,'k','linewidth',3)


%% Parameters for population-scale dynamics

% Initial conditions.

I0 = 1; %Initial number of infected individuals
S0 = N - I0; %Initial number of susceptibles


% Times at which population-scale dynamics are to be calculated.

tmax = 80; %Maximum time
dt = 0.001; %Time step
t_vector = 0:dt:tmax; %Time grid


%% Population-scale dynamics using compartmental method

% Parameters for compartmental framework.

n_vector = [10,20,50]; %Values of the number of compartments, n
T = 7; %Expected infectiousness very small for greater times since infection 


% Use the compartmental method to calculate the population-scale dynamics,
% for the values of the number of compartments, n, specified by n_vector.

figure(2); hold on; %Plot results in Figure 2

for n = n_vector
    [~,dS_dt_vector_compartmental] = compartmental_solution(x_vector,beta_vector,n,T,S0,I0,t_vector);
    plot(t_vector,-dS_dt_vector_compartmental,'linewidth',3)
end


%% Population-scale dynamics using IDE method

[~,dS_dt_vector_IDE] = IDE_solution(x_vector,beta_vector,S0,I0,tmax,dt);
plot(t_vector,-dS_dt_vector_IDE,'k:','linewidth',3)


%% Format figures

figure(1);
set(gcf,'Position',[360 278 560 560])
ax1 = gca;
ax1.FontSize = 24;
ax1.TitleFontSizeMultiplier = 1;
ax1.LabelFontSizeMultiplier = 1;
ax1.FontWeight = 'bold';
ax1.LineWidth = 1.5;
axis square
xlabel('Time since infection (days)');
ylabel('Expected infectiousness (day^{-1})');
xlim([0,7])
xticks(0:7);
ylim([0,1.02e-3])

figure(2);
set(gcf,'Position',[360 278 560 560])
ax1 = gca;
ax1.FontSize = 24;
ax1.TitleFontSizeMultiplier = 1;
ax1.LabelFontSizeMultiplier = 1;
ax1.FontWeight = 'bold';
ax1.LineWidth = 1.5;
axis square
xlim([0,80])
xticks(0:20:80)
xlabel('Time since start of outbreak (days)');
ylabel('Rate of new cases (day^{-1})');

legendstr = [];
for n = n_vector
    legendstr = [legendstr,strcat("{\itn=}",num2str(n)," compartments")];
end
legendstr = [legendstr,'IDE model'];
legend(legendstr,'Location','northeast')