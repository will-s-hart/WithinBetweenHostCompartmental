%% This code accompanies the manuscript entitled "A theoretical framework
%% for transitioning from patient-level to population-scale epidemiological
%% dynamics: influenza A as a case study" by Hart et al. For further information about the paper or this
%% code, please email william.hart@keble.ox.ac.uk

%% We request that users cite the original publication when referring to
%% this code or any results generated from it.

%% This code reproduces the panels in Figure 2 of our paper.

clear all; close all; clc;


%% Patient-level dynamics

% Load patient-level data (the vector V_mean_vector contains the mean viral
% load, calculated over 10,000 within-host realisations, at the times since
% infection contained in x_vector).

load('patient_level_data.mat','x_vector','V_mean_vector')


% Calculate expected infectiousness, beta_vector, at times since infection
% x_vector, assuming infectiousness is proportional to viral load.

R0 = 1.5; %Basic reproduction number
N = 1000; %Population size
beta_vector = R0*V_mean_vector/N;


% Plot expected infectiousness curve

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

n = 1000; %Values of the number of compartments, n
T = 7; %Expected infectiousness very small for greater times since infection 


% Use the compartmental method to calculate the population-scale dynamics,
% for the values of the number of compartments, n, specified by n_vector.

[~,dS_dt_vector_compartmental] = compartmental_solution(x_vector,beta_vector,n,T,S0,I0,t_vector);

% Plot results

figure(2); hold on;
plot(t_vector,-dS_dt_vector_compartmental,'color',[0,0.5,1],'linewidth',3)


%% Population-scale dynamics using IDE method

[~,dS_dt_vector_IDE] = IDE_solution(x_vector,beta_vector,S0,I0,tmax,dt);

% Plot results

plot(t_vector,-dS_dt_vector_IDE,'k--','linewidth',3)

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
ylim([0,26])
yticks(0:5:25)
xlabel('Time since start of outbreak (days)');
ylabel('Rate of new cases (day^{-1})');

l=legend({'Compartmental method','IDE method'},'Location','northeast');
l.Position=[0.3929    0.8750    0.5125    0.1027];