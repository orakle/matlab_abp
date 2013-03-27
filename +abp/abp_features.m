function [m, s, c] = abp_features(varargin)
% ABP_FEATURES - Mean and stdev of ABP features across conditions
%
% [m, s, c] = abp_features(file);
% abp_features(files);
%
% Where
%
% FILE is an EDF+ file containing two data channels. The first channel
% contains the ABP time series, and the second contains an annotation
% signal indicating the transitions between experimental conditions. 
%
% C is a Kx1 vector with condition codes.
%
% M and S are is a Kx9 numeric matrices. The former contains the mean value
% of 9 different ABP features across the K conditions. The latter contains
% the standard deviations. The features in M and S are:
%
%       Col 1: Systolic BP [mmHg]
%           2: Diastolic BP [mmHg]
%           3: Pulse Pressure [mmHg]
%           4: Mean Pressure [mmHg]
%           5: mean_dyneg
%           6: Area under systole 0.3*sqrt(RR)  method
%           7: Area under systole 1st min-slope method
%           8: Heart rate
%           9: Liljestrand and Zander's cardiac output estimate
%
% If no output arguments are provided, then abp_features will generate a
% .csv file (abp_features_XXXXXX.csv) that contains the computed ABP
% features.
%
% See also: cardiac_output

import io.edfplus.read;
import cardiac_output.*;
import datahash.DataHash;
import abp.abp_features;

STD_PSYS  = 120;
STD_PDIAS = 80;
MIN_COND_DUR = 5000; % In samples
DISCARD = 1000;     % Samples to discard from the annotations
NB_COND = 5;

if nargin > 1,
    m = cell(1, nargin);
    s = cell(1, nargin);
    c = cell(1, nargin);
    for i = 1:nargin
        [thisM, thisS, thisC] = abp_features(varargin{i});
        m{i} = thisM; 
        s{i} = thisS;
        c{i} = thisC;
    end
    if nargout < 1,
        save2file(varargin, m,s,c);
    end
    return;
end

file = varargin{1};

[hdr, dat] = read(file, 'verbose', false);

bpIdx = find(ismember(hdr.channel_type, 'BP'));
if numel(bpIdx) ~= 1,
    error('There must be exactly one BP signal: %d found', numel(bpIdx));
end

% Must use 100 Hz for James Sun's functions
abp = resample(dat(bpIdx,:)', 125, hdr.sr(1));
ann = resample(dat(end,:)', 125, hdr.sr(1));
ann = medfilt1(ann, 20);

% Determine onsets of each experimental condition
% This is very fragile. Try something smarter? 
idx = find(abs(diff(flipud(ann))) > 0.25);
idx(idx < DISCARD) = [];
condOnset = [];
prevOnset = -Inf;
for i = 1:numel(idx),
    if idx(i)-prevOnset > MIN_COND_DUR,
        condOnset = [condOnset; idx(i)]; %#ok<AGROW>
        prevOnset = idx(i);
    end
end
condOnset = numel(ann)-condOnset;
condOnset = flipud(condOnset(1:NB_COND+1));
% Approximately in mmHg
ssf = abp(1:1000);
abp = (STD_PSYS-STD_PDIAS)*(abp - mean(ssf))/range(ssf)+STD_PDIAS;

% Detect beat onsets and extract features for each beat
r = wabp(abp);
out = abpfeature(abp, r);

% Estimate heart rate
t1 = [0;diff(out(:, 1))];
t2 = [0;diff(out(:, 3))];
t = (t1 + t2)/2;
t = medfilt1(t, 5);
hr = 60./(t/100);

% Select relevant columns and remove outliers
out = out(:, [2 4:6 8 10 12]);
inRange = r > idx(1) & r < idx(end);
out(inRange, :) = medfilt1(out(inRange, :), 4);

% HR and CO
out = [out hr (out(:,3)./(out(:,1)+out(:,2))).*hr];

nCond = numel(condOnset)-1;
m = nan(nCond, size(out,2));
s = nan(nCond, size(out,2));
c = nan(nCond, 1);
t = 1:numel(ann);
for i = 1:nCond
    inRange = r > condOnset(i) & r < condOnset(i+1);
    c(i) = round(mean(ann(t > condOnset(i) & t < condOnset(i+1))));
    this = out(inRange, :);
    m(i, :) = mean(this);
    s(i, :) = std(this);                    
end

% First and last are always condition 0
c(1)   = 0;
c(end) = 0;

if nargout < 1,
    save2file({file}, {m}, {s}, {c});
end

end


function save2file(file, m, s, c)
import datahash.DataHash;
import safefid.safefid;
import somsds.parse_filename;
import mperl.join;

filesHash = DataHash(file);
fName = ['abp_features_' filesHash(1:6) '.csv'];

fid = safefid.fopen(fName, 'w');

% print file header
tags = parse_filename(file{1});
tagNames = fieldnames(tags);
tags = join(',', tagNames);
fprintf(fid, '%s,counter,%s,%s,%s,%s,%s,%s,%s,%s\n', tags, ...
    'DiastolicBP',       ...
    'PulsePressure',     ...
    'MeanPressure',      ...
    'MeanDyneg',         ...
    'AreaUnderSystole1', ...
    'AreaUnderSystole2', ...
    'HeartRate',         ...
    'CO');

% Print features to file
for i = 1:numel(file)
    
    tags = parse_filename(file{1});
    tagVals = cellfun(@(x) tags.(x), tagNames, 'UniformOutput', false);
    tagVals = cellfun(@(x) join(',', x), tagVals, 'UniformOutput', false);
    for k = 1:size(m{i}, 1)
        fprintf(fid, '%s,%d,%s,%s\n', join(',', tagVals), c{i}(k), ...
            join(',', m{i}(k,:)), ...
            join(',', s{i}(k,:)));
    end
    
end

end
