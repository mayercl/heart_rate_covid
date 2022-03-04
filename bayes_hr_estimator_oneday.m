function steps_all = bayes_hr_estimator_Snyder(hr_data, steps_data, bin_size, id)

%% Close figures

close all

%% Inputs: 

% 1. hr_data - a two column array, the first column lists the epoch
% time date in days (e.g., 738157 is 01-Jan-2021), the second column lists
% the heart rate value at that time
% 2. steps_data - a two column array like the heart rate data, but the
% second column is the steps value at the respective time
% 3. bin_size - the number of minutes to bin the data (5 is suggested)

%% Outputs: 

% A results.csv file that lists the parameter estimates for each day there
% is data along with the uncertainties in the estimates and the number of
% data points for that day and the daily step count for that day


%% Global parameters

% These parameters are used in the get_likelihood_ar1 function down below
% the main script

global hr_avg step_avg jumps N bin_size1

bin_size1 = bin_size;

%% Initialize parameters

% Customizeable parameters for TMCMC
% Number of samples per estimate, burn-in ratio (usual MCMC burn-in)
num_samples = 100000;
burnin_ratio = 0.5;

% Credible interval should contain this fraction of samples
confidence = 0.8;

% nbars = 48; % Number of bars per day (so 48 is half-hour activity bars) 
% thresh = 250; % Number of steps in bin for bar to be full
% minbar = 30; % Bins with fewer than this number of steps aren't shown
% nytick = 15; % Tick every this number of days on the Y Axis


% Optimal parameter values and uncertainties are stored in these lists

all_fits = [];
all_stds = [];
phases = [];

%% Open the results file

% Open results.csv file and add a header line

results_file = fopen(['results_snyder/' id '_results.csv'], 'w');

fprintf(results_file, 'Date (dd-mmm-yyyy), Constant Vertical Shift in HR, Amplitude of Circadian Oscillation, Horizontal Shift (Phase Estimate-Raw), Phase Estimate (Time), Heart Rate Increase per Step, Uncorrelated noise (standard deviation), Error correlation decay rate, STD of constant shift parameter, STD of amplitude, Uncertainty in phase (80 confidence), STD of HRpS, STD of uncorrelated noise, STD of correlated noise, Number of data points, MSE, Daily step count \n');

%% Process data

hr_dates = hr_data(:,1);

steps_dates = steps_data(:,1);

hr_vals = hr_data(:,2);
steps_vals = steps_data(:,2);

% Save all the days in the interval of data given to loop over

hr_days = floor(hr_dates);

dmy_hr = unique(hr_days, 'rows', 'stable');

dmy_hr = (dmy_hr(1):(dmy_hr(end)+1));

num_days = length(dmy_hr)-1;

%% Remove hr_times and hr_values where steps are zero consecutively for more than 2 hours

[hr_dates, hr_vals] = constant_steps_filter2(hr_dates, hr_vals, steps_dates, steps_vals,120);

%% Main loop

l = 1;
steps_all = [];
for i = 1:num_days
    
    % Load steps information for that day
    
    steps_times = steps_dates(steps_dates>=dmy_hr(i) & steps_dates < dmy_hr(i+1));
    
    steps_values = steps_vals(steps_dates>=dmy_hr(i) & steps_dates < dmy_hr(i+1));
        
    % Load heart rate data for the specific day
    
    hr_times = hr_dates(hr_dates < dmy_hr(i+1) & hr_dates >= dmy_hr(i));
    
    hr_values = hr_vals(hr_dates < dmy_hr(i+1) & hr_dates >= dmy_hr(i));

    % Convert the times data to minutes

    hr_times = hr_times*24*60;
    steps_times = steps_times*24*60;

    daily_step = sum(steps_values);
    
    % Ask whether the hr data is empty for that specific day. If so, make
    % the fit values NaN and continue in the loop
    
    if(isempty(hr_times))
        
        all_fits = [all_fits; NaN NaN NaN NaN NaN NaN];
        all_stds = [all_stds; 0];
        if(i==1)
            
            phases = [];
        else
            
            phases = [phases; phase_est];
            
        end
        
        last_ten(2:10) = last_ten(1:9);
        last_ten_var(2:10) = last_ten_var(2:10);
        last_ten(1) = last_ten(2);
        last_ten_var(1) = last_ten_var(2);
        
        continue;
        
    else

        raw_steps_data = [steps_times steps_values];
         
        raw_hr_data = [hr_times hr_values];
        
        % Find all unique binned times for HR data
        [hr_avg, ~, idx] = unique(floor(raw_hr_data(:, 1) / bin_size), 'stable');

        % Average HR measurements in same bin
        val = accumarray(idx, raw_hr_data(:, 2), [], @mean);

        % Combine unique times and averages
        hr_avg = [hr_avg, val];

        try

            % Find leftmost and rightmost bin
            left_min = floor(min(raw_hr_data(1, 1), raw_steps_data(1, 1)) / bin_size);
            right_max = floor(max(raw_hr_data(end, 1), raw_steps_data(end, 1)) / bin_size);
            clearvars raw_hr_data;

        catch

            left_min = floor(raw_hr_data(1, 1) / bin_size);
            right_max = floor(raw_hr_data(end, 1) / bin_size);
            clearvars raw_hr_data;

        end

        % Total length of interval containing HR data
        period_offset = [left_min, right_max - left_min + 1] * bin_size;

        % Precomputed quantities for efficiency: gaps in consecutive measurements
        jumps = [hr_avg(1, 1) + 0.1; hr_avg(2:end, 1) - hr_avg(1:(end - 1), 1)];

        % " " " ": number of measurements after averaging same times
        N = size(hr_avg, 1);

        % Fill in any gaps in step data with zeros
        step_int = int32(raw_steps_data(:, 1));
        steps_new = zeros(period_offset(2), 2);
        steps_new(:, 1) = period_offset(1) + (0:(period_offset(2) - 1));
        steps_new(step_int - period_offset(1) + 1, 2) = raw_steps_data(:, 2);
        clearvars raw_steps_data;

        % Average steps data into bins (including zeros)
        [step_avg, ~, idx] = unique(floor(steps_new(:, 1) / bin_size) , 'stable');
        val = accumarray(idx, steps_new(:, 2), [], @mean);
        val2 = accumarray(idx, steps_new(:, 2), [], @sum);
        
        step_avg = [step_avg, val];
        val2_all = [step_avg(:,1), val2];
        % We only need to keep the step bins corresponding to HR data
        step_avg = step_avg(hr_avg(:, 1) - step_avg(1, 1) + 1, :);
        val2_all = val2_all(hr_avg(:,1) - step_avg(1,1) + 1,:);
        
        if(length(hr_avg(:,1))<20)
            all_fits = [all_fits; NaN NaN NaN NaN NaN NaN];
            all_stds = [all_stds; 0];

            if(l)
                phases = [];
            else
                phases = [phases; phase_est];
            end
            
            if(~l)
                last_ten(2:10) = last_ten(1:9);
                last_ten_var(2:10) = last_ten_var(2:10);
                last_ten(1) = last_ten(2);
                last_ten_var(1) = last_ten_var(2);
            end

            continue;            
        end

        % On the first day, fit data from scratch
        if l == 1

            % On the first day, estimate the initial parameters using the
            % particleswarm global optimization procedure on the likelihood
            
            fun = @(x) -0.5 * get_likelihood_ar1(x);

            optfit = particleswarm(@(x) get_likelihood_ar1(x), 6, [min(hr_avg(:,2)) 0 0 0 0 0], [max(hr_avg(:,2)) 20 1440/bin_size 2 15 1]);

            % Tracking the most recent fitted phases. These are only used to guess
            % starting points for the following day, so that convergence ideally
            % happens somewhat faster -- should not affect final results
            last_ten = ones(10,1) * optfit(3);
            last_ten_var = ones(10,1) * 40000;
            
            l = 0;


        % On future days, likelihood should incorporate previous day (plus
        % Gaussian noise, sd 1 hr) as a prior
        else
            fun = @(x) -0.5 * get_likelihood_ar1(x) - ((mod(x(3) - optfit(3) + 720 / bin_size, 1440 / bin_size) - 720 / bin_size) ^ 2) / (2 * ((stdfit / norminv(0.5 * (1 + confidence))) ^ 2 + (60 / bin_size) ^ 2));
        end

        % Start by guessing a repeat of the previous day's fit
        params = optfit;

        % Specifically for phase (parameter of interest), refine the guess
        % according to the last 10 phases weighted by their certainty
        params(3) = circ_mean(last_ten*2*pi/(24*60/bin_size), 1./last_ten_var)*(24*60/bin_size)/(2*pi);

        % Make sure phase is in valid range
        while params(3) > 2880 / bin_size
            params(3) = params(3) - 1440 / bin_size;
        end
        while params(3) < 1440 / bin_size
            params(3) = params(3) + 1440 / bin_size;
        end
        lb = params(3) - 720 / bin_size;
        rb = params(3) + 720 / bin_size;

        % Reasonable minimum values to prevent divergence of the likelihood
        params([1 2 4 5 6]) = max(params([1, 2, 4, 5, 6]), [50, 0.5, 0.05, 0.5, 0.05]);


        % TMCMC parameters. Number of chains to run in parallel:
        num_walkers = 500;

        % The prior should ensure only physical values of the parameters
        fun_pri = @(x) logical(prod(x([1, 2, 4, 5, 6]) > 0)) && (x(6) < 1);

        % Start walkers off as random perturbations around the starting guess
        start_locs = repmat(params', 1, num_walkers) .* (1 + 0.1 * randn(6, num_walkers));

        % Ensure that starting locations are in valid range
        start_locs(6, :) = min(start_locs(6, :), 0.999);
        start_locs(3, :) = max(min(start_locs(3, :), rb - 0.1) , lb + 0.1);

        % Run GWMCMC to sample from the likelihood
        [samples, ~] = gwmcmc_periodic(start_locs, {fun_pri fun}, num_samples, 'BurnIn', burnin_ratio, 'BinSize', bin_size);

        % Best fit is posterior mean of all samples
        optfit= mean(samples, [2 3])';
        optfit_std = std(samples, 1, [2 3])';

        % To avoid periodic boundary effects, phase was allowed to wander the
        % whole real line. At the end, map back to 1 day interval and take
        % the circular mean to get the estimate

        ss = samples(3,:,end);
        optfit(3) = mod(circ_mean(ss'*2*pi*bin_size/1440),2*pi)*24*60/(2*pi*bin_size);

        %% Compute the uncertainty in the phase parameter

        % Finding the uncertainty: start by assuming an uncertainty of 6 hours
        % in either direction, then use binary search (10 iterations) to find interval that
        % contains ``confidence'' fraction of all samples
        stdfit = 6 * 60 / bin_size;
        std_jump = 0.5 * stdfit;
        search_space = mod(samples(3, :, :), 24 * 60 / bin_size);
        search_space = abs(search_space-optfit(3));
        search_space_angles = search_space*2*pi*bin_size/(24*60);
        for j = 1:10
            if mean(search_space_angles<(stdfit*bin_size*2*pi/(24*60)), 'all') < confidence
                stdfit = stdfit + std_jump;
            else
                stdfit = stdfit - std_jump;
            end
            std_jump = 0.5 * std_jump;
        end 

        % Add fits and uncertainties to list
        
        phase_est = mod(optfit(3) /(60 / bin_size) - 12, 24);
        phases = [phases; phase_est];
        phase_est_str = datestr(phase_est/24, 'HH:MM');

        all_fits = [all_fits; optfit];
        all_stds = [all_stds; stdfit];

        %% Write day's results to the results.csv file

        num_points = length(hr_avg(:,1));

        mse = mean((hr_avg(:, 2) - (optfit(1) + optfit(2) * cos(pi * (hr_avg(:, 1) - optfit(3)) / (24 * 30 / bin_size1)) + optfit(4) * step_avg(:,2))).^2);
        fprintf(results_file, '%s, %f, %f, %f, %s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %i, %f, %i\n', datestr(dmy_hr(i)), optfit(1:3), phase_est_str, optfit(4:6), optfit_std(1:2), stdfit/12, optfit_std(4:6), num_points, mse, daily_step);
        
        steps_all = [steps_all; val2_all(:,2)];
        
        figure(1)
        set(gcf, 'Position',  [100, 100, 800, 600])
        set(gcf,'color','w');
        set(gca,'Visible','off')
        set(0,'DefaultAxesTitleFontWeight','normal');
        scatter(hr_avg(:,1),hr_avg(:,2),'filled')
        hold on 
        plot(hr_avg(:,1), optfit(1) + optfit(2) * cos(pi * (hr_avg(:, 1) - optfit(3)) / (24 * 30 / bin_size1)) + optfit(4) * step_avg(:,2),'r')
        hold off 
        
        last_ten(2:10) = last_ten(1:9);
        last_ten_var(2:10) = last_ten_var(2:10);
        last_ten(1) = optfit(3);
        last_ten_var(1) = stdfit ^ 2;

    end
    
end


end

%% Likelihood function

function likelihood = get_likelihood_ar1(params)

global hr_avg step_avg jumps N bin_size1

steps = step_avg(:,2);

% Difference between observed data and prediction without noise
r = hr_avg(:, 2) - (params(1) + params(2) * cos(pi * (hr_avg(:, 1) - params(3)) / (24 * 30 / bin_size1)) + params(4) * steps);

% When there are gaps in HR, noise is the convolution of several
% individual steps of noise. Prepower computes these up to a max gap
% size of 10 (for computational efficiency). Terms decay fairly quickly
% so in practice 10 was enough to lose correlation completely
prepower = ones(N, 1);
b_max = min(max(jumps(2:end)) - 1, 10);
for i = 1:b_max
    prepower(jumps > i) = prepower(jumps > i) + params(6) ^ (2 * i);
end
prepower(1) = 1;

% Contribution to likelihood due to error at each step
likelihood = sum((r(2:end) - (params(6) .^ jumps(2:end)) .* r(1:(end - 1))) .^ 2 ./ (prepower(2:end))) / (params(5) * params(5));

% Contribution to likelihood due to relation between consecutive errors
% (here using AR(1) error model)
likelihood = likelihood + sum(log(params(5) * params(5) * prepower)) - log(1 - params(6) * params(6)) + r(1) * r(1) * (1 - params(6) * params(6)) / (params(5) * params(5));

end