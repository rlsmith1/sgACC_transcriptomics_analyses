function [trdata, tedata] = permute_data(res, trdata, tedata, split, iperm)
% permute_data
%
% # Syntax
%   [trdata, tedata] = permute_data(res, trdata, tedata, split, iperm)
%
%_______________________________________________________________________
% Copyright (C) 2022 University College London
%
% Written by Agoston Mihalik (cca-pls-toolkit@cs.ucl.ac.uk)
% $Id$
%
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
%
% Modified by Agoston Mihalik (am3022@cantab.ac.uk)
%   refactor code a bit
%   update for GRCCA analysis

% Load cfg
cfg = loadmat(res, fullfile(res.dir.frwork, 'cfg.mat'), 'cfg');

% Load training and test indexes
[otrid, oteid] = loadmat(res, fullfile(res.dir.frwork, 'outmat.mat'), 'otrid', 'oteid');
trid = otrid(:,split);
teid = oteid(:,split);

% Subjects permutations index
permid = loadmat(res, fullfile(res.dir.frwork, 'perm', sprintf('permmat_%d.mat', ...
    res.stat.nperm)), 'permid');
sid = permid{res.frwork.split.all==split}(:,iperm);

if strfind(cfg.stat.perm, 'train')
    % Subject orders within training set
    if strcmp(cfg.stat.perm, 'train+test')
        [~, trsid] = sort(sid(trid));
        [~, trsid] = sort(trsid);
    else
        trsid = sid;
    end

    % Permute training subjects
    switch cfg.machine.name
        case {'pls' 'spls'}
            trdata.Y = trdata.Y(trsid,:);

        case {'cca' 'rcca' 'grcca'}
            trdata.RY = trdata.RY(trsid,:);
    end
end

if strfind(cfg.stat.perm, 'test')
    % Subject orders within test set
    [~, tesid] = sort(sid(teid));
    [~, tesid] = sort(tesid);
    
    % Permute test subjects
    switch cfg.machine.name
        case {'pls' 'spls'}
            tedata.Y = tedata.Y(tesid,:);
            
        case {'cca' 'rcca' 'grcca'}
            tedata.RY = tedata.RY(tesid,:);
    end
end
