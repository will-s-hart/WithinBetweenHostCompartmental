function [S_vector,dS_dt_vector] = compartmental_solution(x_vector,beta_vector,n,T,S0,I0,t_vector)
% compartmental_solution  Uses a compartmental method to transition from
% patient-level to population-scale epidemiological dynamics.
%   [S_vector,dS_dt_vector]
%       = compartmental_solution(x_vector,beta_vector,n,T,S0,I0,t_vector)
%   returns vectors S_vector and dS_dt_vector, containing the number of
%   susceptible individuals and the rate of change of the number of
%   susceptible individuals, respectively, at the times contained in
%   t_vector. The input vector beta_vector contains the expected
%   infectiousness at the times since infection contained in x_vector; n
%   is the number of infected compartments; T is a time since infection
%   beyond which the expected infectiousness is very small; S0 and I0 give
%   the intial numbers of susceptible and infected individuals, respectively.
%
%   See also IDE_solution.


% Calculate transmission rates in SI_{n}R model.

int_beta_vector = cumtrapz(x_vector,beta_vector);
int_beta_eval = interp1(x_vector,int_beta_vector,[0:(T/n):((n-1)*T/n),max(x_vector)]','linear',0);
beta_compartments = (n/T)*diff(int_beta_eval);


% Rates of transition between infected compartments.

mu_compartments = (n/T)*ones(n,1);


% Initial conditions

I0_compartments = zeros(n,1);
I0_compartments(1) = I0;


% Use ode45 to solve compartmental model. Note y(1) = S, y(2) = I_{1},
% y(3) = I_{2),..., y(n+1) = I_{n}

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[~,y] = ode45(@(t,y) SInR_RHS(t,y,n,beta_compartments,mu_compartments),t_vector,[S0;I0_compartments],options);


% Number of susceptibles, and rate of change of number of susceptibles.

S_vector = y(:,1)';
dS_dt_vector = -(y(:,2:end)*beta_compartments)'.*S_vector;


function ydot = SInR_RHS(t,y,n,beta_compartments,mu_compartments)
% SInR_RHS  Produces the right hand side of the system of
% differential equations in the SI_{n}R model.
%   ydot = SInR_RHS(t,y,n,beta_compartments,mu_compartments) returns the
%   rate of change of y at time t, where y(1) = S, y(2) = I_{1},
%   y(3) = I_{2),..., y(n+1) = I_{n}; n is the number of infected
%   compartments (this function requires n > 1); beta_compartments are the
%   transmission rate parameters in the SI_{n}R model; mu_compartments are
%   rates of transition through the infected compartments.


ydot(1) = -beta_compartments'*y(1)*y(2:n+1);
ydot(2) = -ydot(1) - mu_compartments(1)*y(2);
ydot(3:(n+1)) = mu_compartments(1:n-1).*y(2:n) - mu_compartments(2:n).*y(3:n+1);
ydot = ydot';
end
end