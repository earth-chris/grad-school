% Matlab(R) code for 
% "day-to-day model" 
% as presented in 
% "A statistically predictive model for future monsoon failure in India"
% by J. Schewe and A. Levermann
%
% Author:	Jacob Schewe 
% Address:  Potsdam Institute for Climate Impact Research, Potsdam, Germany
% email:	jacob.schewe@pik-potsdam.de
% Date:		03 Jun 2012
%
% Please refer to the manuscript for a discussion of the model. 
%
clear all

% --- Setup ---
% --- Parameters
nyears = 6030; % Number of years (seasons) to run
lseason = 135; % Length of season in days
Pstrong  = 9.; % Precipitation in strong (wet) state, mm/day
Pweak    = 0.; % Precipitation in weak (dry) state, mm/day
tau 	   = 17; % Memory length in time steps (days) 
prmax   = 0.8; % Maximum probability of either state
p_init = 0.75; % Initial probability of strong state - p_init <= prmax
% --- Random numbers
% % Model results will differ from run to run because of the stochastic 
% % nature of the model. For robust results, average over many runs. To
% % obtain reproducable results, however, a constant stream of random 
% % numbers must be used - for this purpose, uncomment the following 
% % lines: 
% s = RandStream.create('mt19937ar','seed',5489);
% RandStream.setDefaultStream(s);

% --- Compute ---
Psum = zeros(1,nyears);
% --- Loop over years (seasons)
for i = [1:nyears]
    P = zeros(1,lseason);
    n=1;
	% --- Loop over days of season
    while n <= lseason
        pr = rand; % draw random number 0 <= pr <= 1
		% --- Determine probability
        if n>tau
            p = (sum(P(n-tau:n-1))/tau - Pweak)/(Pstrong-Pweak); % ...from memory effect
        else
            p = p_init; % ...constant in the beginning of the season
		end
		% --- Limit high and low probabilities
		if p > prmax
			p = prmax;
		elseif p < (1-prmax)
			p = (1-prmax);
		end
		% --- Assign precipitation
        if pr < p
            P(n) = Pstrong; 
        else
            P(n) = Pweak; 
        end
        n=n+1;
    end
	% --- Seasonal mean
    Psum(i) = sum(P)/lseason;
end

% --- Plot ---
bins = 0:0.2:10; % Intervals for histogram bins
figure;
hist(Psum,bins)
xlim([-2,12])
ylim([0 500])
title('day-to-day model')
xlabel('mm/day')

% --- END ---
