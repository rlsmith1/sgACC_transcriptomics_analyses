function data = transform_data(cfg, transform, data, mod, varargin)
% transform_data
%
% # Syntax:
%   data = transform_data(cfg, transform, data, mod, varargin)
% 
%_______________________________________________________________________
% Copyright (C) 2023 MQ: Transforming Mental Health
%
% Written by Agoston Mihalik (am3022@cantab.ac.uk)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Modified by Agoston Mihalik (am3022@cantab.ac.uk)
%   update for GRCCA analysis

p = parse_input([], varargin{:});

% Load data modality group
if isfield(cfg.data.(mod), 'group')
    S = load(cfg.data.(mod).group.fname);
    if isfield(S, 'group')
        group = S.group;
    else
        error('group vector not available in %s.', m, cfg.data.(mod).group.fname);
    end
    
    % Renumber group if needed 
    [~, ~, group] = unique(group);

    % Calculate number of groups
    ngroups = numel(unique(group));
end

if contains(transform, 'svd')
    %----- RCCA transform to solve RCCA problem in principal components basis

    % Transform features
    if isfield(p, 'V')
        if strcmp(transform, 'svd-inverse-transform') && ~isfield(data, mod)
            % Inverse transform weight vector(s) into principal component basis
            data = p.V' * data;
        elseif strcmp(transform, 'svd-left-transform') && ~isfield(data, mod)
            % Transform weight vector(s) from principal component basis into feature space
            data = p.V * data;
        elseif strcmp(transform, 'svd-right-transform') && isfield(data, mod)
            % Transform data into principal component basis
            data.(['R' mod]) = data.(mod) * p.V;

            % Define fields to be saved
            if isfield(p, 'fname') && strncmp(p.fname, 'te_', 3)
                fields = {['R' mod]};
            end
        else
            error('SVD transformation is not as expected.')
        end

    elseif isfield(p, 'fname') && (strncmp(p.fname, 'tr_', 3) || strncmp(p.fname, 'boot_', 5))
        if strcmp(transform, 'svd-right-transform') && isfield(data, mod)
            % Calculate PCA transformation
            [data.(['V' mod]), data.(['R' mod]), data.(['L' mod])] = fastsvd(data.(mod), ...
                0, cfg.machine.svd.tol, cfg.machine.svd.(['var' lower(mod)]), 'V', 'R', 'L');

            % Define fields to be saved
            fields = {['V' mod] ['R' mod] ['L' mod]};
        else
            error('SVD transformation is not as expected.')
        end
    end

    % Remove original data to free up RAM
    if isfield(data, mod)
        data = rmfield(data, mod);
    end

    % Save SVD results
    if exist('fields', 'var') && (cfg.env.save.svd || (strcmp(cfg.machine.name, 'grcca') && ~exist('group', 'var') ...
            && isfield(p, 'fname') && ~exist_file(cfg, fullfile(cfg.dir.load, 'svd', p.fname))))
        name_value = parse_struct(data, 1, fields);
        savemat(cfg, fullfile(cfg.dir.load, 'svd', p.fname), name_value{:});
    end

elseif contains(transform, 'grcca') && exist('group', 'var')
    %----- GRCCA transform to reduce GRCCA problem to RCCA problem (see Tuzhilina et al. 2021)

    % Check inputs
    if any(~isfield(p, {'lambda' 'mu'}))
        error('Missing input arguments for GRCCA transformation')
    end

    for g=1:ngroups
        nel = sum(group==g);

        if strcmp(transform, 'grcca-inverse-transform')
            % Calculate block-wise inverse transformation matrix
            invA = [eye(nel)-repmat(1/nel, nel) repmat(sqrt(p.mu / p.lambda / nel), nel, 1)]';

            % Inverse transform weight vector(s)
            if ~isfield(data, mod)
                data(end+1,:) = invA(end,:) * data(group==g,:);
                data(group==g,:) = invA(1:end-1,:) * data(group==g,:);
            else
                error('GRCCA inverse transformation is not as expected.')
                % data.(mod)(:,group==g) = data.(mod)(:,end-ngroups+g) * invA(end,:) + ...
                %     data.(mod)(:,group==g) * invA(1:end-1,:);
                % data.(mod)(:,end-ngroups+g) = [];
            end
        else

            % Calculate block-wise transformation matrix
            A = [eye(nel)-repmat(1/nel, nel) repmat(sqrt(p.lambda / p.mu / nel), nel, 1)];

            if strcmp(transform, 'grcca-left-transform') && ~isfield(data, mod)
                % Transform weight vector(s)
                data(group==g,:) = A * data([group==g; true],:);
                data(numel(group)+1,:) = [];
            elseif strcmp(transform, 'grcca-right-transform') && isfield(data, mod)
                % Transform data
                data.(mod)(:,end+1) = data.(mod)(:,group==g) * A(:,end);
                data.(mod)(:,group==g) = data.(mod)(:,group==g) * A(:,1:end-1);
            else
                error('GRCCA transformation is not as expected.')
            end
        end
    end
end
