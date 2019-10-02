function [S_vector,dS_dt_vector] = IDE_solution(x_vector,beta_vector,S0,I0,tmax,dt)
% IDE_solution  Uses an IDE method to transition from patient-level to
%   population-scale epidemiological dynamics.
%   [S_vector,dS_dt_vector] =
%       IDE_solution(x_vector,beta_vector,S0,I0,tmax,dt)
%   returns vectors S_vector and dS_dt_vector, containing the number of
%   susceptible individuals and the rate of change of the number of
%   susceptible individuals, respectively, on a time grid with step dt up
%   to maximum time tmax. The input vector beta_vector contains the expected
%   infectiousness at the times since infection contained in x_vector; S0
%   and I0 give the intial numbers of susceptible and infected individuals,
%   respectively.
%
%   See also compartmental_solution.


% Setup time grid.

t_vector = 0:dt:tmax; %Time grid
n = length(t_vector); %Number of time-points


% Evaluate the expected infectiousness on the time grid.

beta_vector = interp1(x_vector,beta_vector,t_vector,'linear',0);


% Use an Euler method to solve the K&M IDE model.

S_vector = zeros(1,n); %Initialise vector containing number of susceptibles

S_vector(1) = S0;
S_vector(2) = S_vector(1)*(1-dt*I0*beta_vector(1));

for i = 2:(n-1)
    S_vector(i+1) = S_vector(i)*(1+dt*(beta_vector(2:i)*(S_vector(i:(-1):2)-S_vector((i-1):(-1):1))'-I0*beta_vector(i)));
end


% Calculate the rate of change of the number of susceptibles.

dS_dt_vector = gradient(S_vector,dt);
end