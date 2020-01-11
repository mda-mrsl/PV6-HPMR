function struct = readBrukerHeader(filename)
%READBRUKERHEADER reads Bruker JCAMP header file as MATLAB structure
%
%   Usage: struct = readBrukerHeader(filename)
%
%   See also READBRUKERSTUDY, READBRUKERFID
%
%   ??/200?, Dustin Ragan
%   06/2019, Keith Michel

if ~nargin, help(mfilename), end

%%
fid = fopen(filename, 'rt');
rawtext = fread(fid, inf, 'char=>char')';
fclose(fid);

% Remove lines starting in $$
rawtext = regexprep(rawtext, '\n\${2}[^\n]+', ' ');
% Find locations of parameter names as words between ##$ and =
% and corresponding values as everything between = and subsequent ##
params = regexp(rawtext, '\#{2}\$*(\w+)=(.*?)(?=\n\#{2})', 'tokens');
headercell = cell(numel(params),2);

%%
for ind = 1:numel(params)
    headercell{ind,1} = params{ind}{1};
    value = params{ind}{2};
    if isempty( regexp(value, '\n', 'once'))
    	% When there is no newline we have either a numeric or string
        value = str2double(value);
        if isnan(value)
            value = params{ind}{2};
        end
    else
        % When there is a newline it is usually preceded by integers 
        % enclosed in parentheses that indicate the parameter shape
        shape = regexp(value, '^\(([\d, ]+)\)', 'tokens');
        if isempty(shape)
            % Parentheses contained something other than numerals,
            % commas and spaces. Return everything as single line string
            value = regexprep(value, '\n', ' ');
        else
            shape = squeeze( str2num(shape{1}{1}));
            value = regexprep(value, '\n', ' ');
            value = regexp(value, '\) (.*)$', 'tokens');
            value = [value{1}{1}, ' '];
            if isempty( str2num(value))
                % This is not a numeric array, capture it as a cell
                value = regexprep(value, '\n', ' ');
                if strcmp(value(1), '(')
                    value = regexp(value, '(\(.*?\))(?= )', 'tokens');
                elseif strcmp(value(1), '<')
                    value = regexp(value, '(\<.*\>)(?= )', 'tokens');
                else
                    value = regexp(value, '(\S+)', 'tokens');
                end
                value = vertcat(value{:});
                if numel(value) == 1
                    value = value{1};
                end
            else
                value = str2num(value);
                if numel(shape)>1
                    value = reshape(value, shape);
                end
            end
        end
    end
    
    if ischar(value)
        value = regexprep(value, '\s+$', ''); % remove trailing spaces
        value = regexprep(value, '^<(.*)>$', '$1'); % remove angle braces
    end
    headercell{ind,2} = value;
end

struct = cell2struct(headercell(:,2), headercell(:,1));

end