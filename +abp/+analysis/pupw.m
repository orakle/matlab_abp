function pupw(varargin)
% PUPW - Extract ABP features from the pupilator physiology files
%
% feat = abp.analysis.pupilator()
% abp.analysis.pupilator()
% abp.analysis.pupilator('file1', 'file2', ...)
%
% Where
%
% FILE1, FILE2, ... are the names of the raw data files from which the
% features will be extracted (using abp.abp_features). If not provided,
% then all physilogy files of the pupw recording will be used.
%
% FEAT is a cell array having as columns individual features and as rows
% individual data records.
%
% If no output argument is specified the features will be stored in a .csv
% file named abp_features_XXXXX.csv where XXXXX is an alphanumeric hash
% generated from the file names that were used to generate the features.
%
% See also: abp

import somsds.link2rec;
import datahash.DataHash;

if nargin < 1,
    
    folder = DataHash(randn(1,100));
    while (exist(folder, 'dir'))
        folder = DataHash(randn(1,100));
    end
    
    files = link2rec('pupw', 'modality', 'physiology', 'condition', ...
        {'morning-supine', 'morning-sitting', 'afternoon-supine', ...
        'afternoon-sitting'}, 'folder', folder);
    
else
    
    files = varargin;
    
end

try
    
    abp.features(files{:});
    
catch ME
    
    if nargin < 1,
        rmdir(folder, 's');
    end
    rethrow(ME);
    
end

if nargout < 1,
    rmdir(folder, 's');
end


end