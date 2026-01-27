function featid = get_featid(data, param, mod)
% get_featid
%
% # Syntax
%   featid = get_featid(data, param, mod)
%
%_______________________________________________________________________
% Copyright (C) 2022 University College London

% Written by Agoston Mihalik (cca-pls-toolkit@cs.ucl.ac.uk)
% $Id$

% This file is part of CCA/PLS Toolkit.
%
% CCA/PLS Toolkit is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CCA/PLS Toolkit is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CCA/PLS Toolkit. If not, see <https://www.gnu.org/licenses/>.

if isfield(param, (['VAR' lower(mod)]))
    if isfield(data, ['G' mod])
        % Reduce dimensionality by number of explained variance per group
        idstart = cumsum([0 data.(['G' mod])]);
        for i=1:numel(idstart)-1
            featid(idstart(i)+1:idstart(i+1)) = cumsum(data.(['L' mod])(idstart(i)+1:idstart(i+1))) / ...
                sum(data.(['L' mod])(idstart(i)+1:idstart(i+1))) <= param.(['VAR' lower(mod)]);
        end
    else
        % Reduce dimensionality by explained variance
        featid = cumsum(data.(['L' mod])) / sum(data.(['L' mod])) <= param.(['VAR' lower(mod)]);
    end
else
    % Keep all variance
    featid = true(size(data.(['L' mod])));
end

if isfield(param, (['PCA' lower(mod)]))
    if isfield(data, ['G' mod])
        if any(data.(['G' mod]) < param.(['PCA' lower(mod)]))
            warning('Dimensionality in one of the groups of data %s is lower than number of PCA components', mod);
        else
            % Reduce dimensionality by number of principal components per group
            idstart = cumsum([0 data.(['G' mod])]);
            for i=1:numel(idstart)-1
                featid(idstart(i)+param.(['PCA' lower(mod)])+1:idstart(i+1)) = false;
            end
        end
    else
        if sum(featid) < param.(['PCA' lower(mod)])
            warning('Dimensionality of data %s is lower than number of PCA components', mod);
        else
            % Reduce dimensionality by number of principal components
            featid(param.(['PCA' lower(mod)])+1:end) = false;
        end
    end
end