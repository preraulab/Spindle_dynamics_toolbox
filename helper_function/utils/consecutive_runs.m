function [run_lengths, run_inds, filtered_vector] = consecutive_runs(data, min_len, max_len, value)
%CONSECUTIVE_RUNS  Extract consecutive runs of a specified value from a binary vector
%
%   Usage:
%       [run_lengths, run_inds, filtered_vector] = consecutive_runs(data, min_len, max_len, value)
%
%   Inputs:
%       data    : 1xN double/logical - binary vector -- required
%       min_len : integer - minimum run length (default: 1)
%       max_len : integer - maximum run length (default: inf)
%       value   : scalar - value to extract runs of (default: extract runs of ones)
%
%   Outputs:
%       run_lengths     : 1xR integer - length of each selected run
%       run_inds        : 1xR cell - start/end indices of each selected run
%       filtered_vector : 1xN double - binary vector with ones at kept-run positions
%
%   Notes:
%       Runs with length in [min_len, max_len] (inclusive) are returned. If
%       value is omitted, runs of ones are extracted.
%
%   Example:
%       lengths = consecutive_runs([1 0 1 1 0 1 1 1 0 0]);
%       % returns [1 3]
%       [lengths, run_inds] = consecutive_runs([0 0 1 1 0 1 1 1 0 0 1], 2, 4, 1);
%       % returns lengths = [2 3] and run_inds = {[3 4], [6 8]}.
%
%   See also: get_chunks, find, histcounts
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%Default lengths impose no filtering
if nargin <2
    min_len=1;
end

if nargin <3
    max_len=inf;
end

%Check for valid lengths
if min_len>=max_len
    error('Min size must be less than max size');
end

%Select a single value if necessary
if nargin == 4
    data = data==value;
end

% Determine length of data
len = length(data);

% Compute maximum number of runs based on min_len
max_runs = floor(len / min_len);

% Initialize output variables
run_inds = cell(1, max_runs); % Preallocate cell array for maximum possible number of runs
run_lengths = zeros(1, max_runs); % Preallocate vector for maximum possible number of runs

% Initialize run variables
cur_run = 0;
cur_len = 0;

% Iterate through binary vector
for ii = 1:len
    if data(ii)
        % Increment current run
        cur_len = cur_len + 1;
    elseif ii>1 && data(ii-1)
        if cur_len >= min_len && cur_len <= max_len
            % End current run
            cur_run = cur_run + 1;
            runs_start = (ii - cur_len);
            runs_end = (ii - 1);
            run_inds{cur_run} = runs_start:runs_end;
            run_lengths(cur_run) = cur_len;
        end
        % Reset current run
        cur_len = 0;
    end
end

if data(len) && cur_len >= min_len && cur_len <= max_len
    % End current run
    cur_run = cur_run + 1;
    runs_start = (len - cur_len + 1);
    runs_end = (len);
    run_inds{cur_run} = runs_start:runs_end;
    run_lengths(cur_run) = cur_len;
end

% Trim output vectors to only include filled cells
run_inds = run_inds(1:cur_run);
run_lengths = run_lengths(1:cur_run);

%Create the filtered binary vector if needed
if nargout == 3
    filtered_vector = zeros(size(data));
    filtered_vector([run_inds{:}]) = 1;
end

end

