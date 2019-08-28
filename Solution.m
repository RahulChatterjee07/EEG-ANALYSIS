%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\User\Downloads\JD032918.csv


%% Initialize variables.
filename = 'C:\Users\User\Downloads\JD032918.csv';
delimiter = ',';

%% Read columns of data as text:

formatSpec = '%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,4,5,6,8,10]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,4,5,6,8,10]);
rawStringColumns = string(raw(:, [2,3,7,9]));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
JD1 = table;
JD1.Id = cell2mat(rawNumericColumns(:, 1));
JD1.Reference = rawStringColumns(:, 1);
JD1.Type = rawStringColumns(:, 2);
JD1.xcm = cell2mat(rawNumericColumns(:, 2));
JD1.ycm = cell2mat(rawNumericColumns(:, 3));
JD1.zcm = cell2mat(rawNumericColumns(:, 4));
JD1.SpecialSensor = rawStringColumns(:, 3);
JD1.Date = cell2mat(rawNumericColumns(:, 5));
JD1.T1423123080000500 = rawStringColumns(:, 4);
JD1.VarName10 = cell2mat(rawNumericColumns(:, 6));

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R;
%% Convert Table to String array
L=table2array(JD1);
%% Convert String array to Cell array
M=cellstr(L);
%% Sorting the Cell array
R=natsortrows(M,2);
%% Convert sorted cell array to Table
P=cell2table(R,...
    'VariableNames',{'Id' 'Reference' 'Type' 'x_cm' 'y_cm' 'z_cm' 'SpecialSensor' 'Date' 'T1423123080000500' 'VarName10'})
%% natsortrows function 
function [X,ndx,dbg] = natsortrows(X,col,varargin)

% Alphanumeric / Natural-Order sort the rows of a cell array of strings (1xN char).
%
% (c) 2012 Stephen Cobeldick
%
% Alphanumeric sort of the rows of a cell array of strings: sorts by both
% character order and also the values of any numbers that occur within the
% strings. The cell of strings must be a matrix. SORTROWS input <col> is
% also supported, so NATSORTROWS is a drop-in replacement for SORTROWS.
%


%% Input Wrangling %%
%
assert(iscell(X),'First input <X> must be a cell array.')
tmp = cellfun('isclass',X,'char') & cellfun('size',X,1)<2 & cellfun('ndims',X)<3;
assert(all(tmp(:)),'First input <X> must be a cell array of char row vectors (1xN char).')
assert(ismatrix(X),'First input <X> must be a matrix (size RxC).')
%
%% Select Columns to Sort %%
%
[m,n] = size(X);
dbg{n} = [];
ndx = 1:m;
drn = {'descend','ascend'};
isn = false;
%
if nargin<2 || isnumeric(col)&&isempty(col)
	vec = n:-1:1;
elseif ischar(col)&&isrow(col)&&any(strcmpi(col,drn))
	% Sort all columns descending/ascending.
	vec = n:-1:1;
	varargin{max(2,end+1)} = col;
elseif isnumeric(col)
	% Sort columns according to the provided indices.
	assert(isreal(col)&&isvector(col),'Second input <col> must be a real numeric vector.')
	assert(all(fix(col)==col)&&all(abs(col)<=n)&&all(col),...
		'Second input <col> must be a vector of column indices into the first input <X>.')
	vec = reshape(col(end:-1:1),1,[]);
	varargin{max(2,end+1)} = [];
	isn = true;
else
	error('Second input <col> must be a numeric vector of indices, or ''ascend''/''descend''.')
end
%
%% Sort Columns %%
%
for k = vec
	if isn
		varargin(end) = drn((3+sign(k))/2);
		k = abs(k); %#ok<FXSET>
	end
	if nargout<3 % faster:
		[~,ids] = natsort(X(ndx,k),varargin{:});
	else % for debugging:
		[~,ids,tmp] = natsort(X(ndx,k),varargin{:});
		[~,idd] = sort(ndx);
		dbg{k} = tmp(idd,:);
	end
	ndx = ndx(ids);
end
%
ndx = ndx(:);
X = X(ndx,:);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsortrows
%% natsort function
function [X,ndx,dbg] = natsort(X,xpr,varargin) %#ok<*SPERR>
% Alphanumeric / Natural-Order sort the strings in a cell array of strings (1xN char).
%
% (c) 2012 Stephen Cobeldick
%
% Alphanumeric sort of a cell array of strings: sorts by character order
% and also by the values of any numbers that are within the strings. The
% default is case-insensitive ascending with integer number substrings:
% optional inputs control the sort direction, case sensitivity, and number
% matching (see the section "Number Substrings" below).
%


%% Input Wrangling %%
%
assert(iscell(X),'First input <X> must be a cell array.')
tmp = cellfun('isclass',X,'char') & cellfun('size',X,1)<2 & cellfun('ndims',X)<3;
assert(all(tmp(:)),'First input <X> must be a cell array of char row vectors (1xN char).')
%
% Regular expression:
if nargin<2 || isnumeric(xpr)&&isempty(xpr)
	xpr = '\d+';
else
	assert(ischar(xpr)&&isrow(xpr),'Second input <xpr> must be a regular expression (char row vector).')
end
%
% Optional arguments:
tmp = cellfun('isclass',varargin,'char') & 1==cellfun('size',varargin,1) & 2==cellfun('ndims',varargin);
assert(all(tmp(:)),'All optional arguments must be char row vectors (1xN char).')
% Character case matching:
ChrM = strcmpi(varargin,'matchcase');
ChrX = strcmpi(varargin,'ignorecase')|ChrM;
% Sort direction:
DrnD = strcmpi(varargin,'descend');
DrnX = strcmpi(varargin,'ascend')|DrnD;
% Relative sort-order of numbers compared to characters:
RsoB = strcmpi(varargin,'beforechar');
RsoA = strcmpi(varargin,'afterchar');
RsoX = strcmpi(varargin,'asdigit')|RsoB|RsoA;
% SSCANF conversion format:
FmtX = ~(ChrX|DrnX|RsoX);
%
if nnz(FmtX)>1
	tmp = sprintf(', ''%s''',varargin{FmtX});
	error('Overspecified optional arguments:%s.',tmp(2:end))
end
if nnz(DrnX)>1
	tmp = sprintf(', ''%s''',varargin{DrnX});
	error('Sort direction is overspecified:%s.',tmp(2:end))
end
if nnz(RsoX)>1
	tmp = sprintf(', ''%s''',varargin{RsoX});
	error('Relative sort-order is overspecified:%s.',tmp(2:end))
end
%
%% Split Strings %%
%
% Split strings into number and remaining substrings:
[MtS,MtE,MtC,SpC] = regexpi(X(:),xpr,'start','end','match','split',varargin{ChrX});
%
% Determine lengths:
MtcD = cellfun(@minus,MtE,MtS,'UniformOutput',false);
LenZ = cellfun('length',X(:))-cellfun(@sum,MtcD);
LenY = max(LenZ);
LenX = numel(MtC);
%
dbg = cell(LenX,LenY);
NuI = false(LenX,LenY);
ChI = false(LenX,LenY);
ChA = char(double(ChI));
%
ndx = 1:LenX;
for k = ndx(LenZ>0)
	% Determine indices of numbers and characters:
	ChI(k,1:LenZ(k)) = true;
	if ~isempty(MtS{k})
		tmp = MtE{k} - cumsum(MtcD{k});
		dbg(k,tmp) = MtC{k};
		NuI(k,tmp) = true;
		ChI(k,tmp) = false;
	end
	% Transfer characters into char array:
	if any(ChI(k,:))
		tmp = SpC{k};
		ChA(k,ChI(k,:)) = [tmp{:}];
	end
end
%
%% Convert Number Substrings %%
%
if nnz(FmtX) % One format specifier
	fmt = varargin{FmtX};
	err = ['The supplied format results in an empty output from sscanf: ''',fmt,''''];
	pct = '(?<!%)(%%)*%'; % match an odd number of % characters.
	[T,S] = regexp(fmt,[pct,'(\d*)([bdiuoxfeg]|l[diuox])'],'tokens','split');
	assert(isscalar(T),'Unsupported optional argument: ''%s''',fmt)
	assert(isempty(T{1}{2}),'Format specifier cannot include field-width: ''%s''',fmt)
	switch T{1}{3}(1)
		case 'b' % binary
			fmt = regexprep(fmt,[pct,'(\*?)b'],'$1%$2[01]');
			val = dbg(NuI);
			if numel(S{1})<2 || ~strcmpi('0B',S{1}(end-1:end))
				% Remove '0B' if not specified in the format string:
				val = regexprep(val,'(0B)?([01]+)','$2','ignorecase');
			end
			val = cellfun(@(s)sscanf(s,fmt),val, 'UniformOutput',false);
			assert(~any(cellfun('isempty',val)),err)
			NuA(NuI) = cellfun(@(s)sum(pow2(s-'0',numel(s)-1:-1:0)),val);
		case 'l' % 64-bit
			NuA(NuI) = cellfun(@(s)sscanf(s,fmt),dbg(NuI)); %slow!
		otherwise % double
			NuA(NuI) = sscanf(sprintf('%s\v',dbg{NuI}),[fmt,'\v']); % fast!
	end
else % No format specifier -> double
	NuA(NuI) = sscanf(sprintf('%s\v',dbg{NuI}),'%f\v');
end
% Note: NuA's class is determined by SSCANF or the custom binary parser.
NuA(~NuI) = 0;
NuA = reshape(NuA,LenX,LenY);
%
%% Debugging Array %%
%
if nargout>2
	dbg(:) = {''};
	for k = reshape(find(NuI),1,[])
		dbg{k} = NuA(k);
	end
	for k = reshape(find(ChI),1,[])
		dbg{k} = ChA(k);
	end
end
%
%% Sort Columns %%
%
if ~any(ChrM) % ignorecase
	ChA = upper(ChA);
end
%
ide = ndx.';
% From the last column to the first...
for n = LenY:-1:1
	% ...sort the characters and number values:
	[C,idc] = sort(ChA(ndx,n),1,varargin{DrnX});
	[~,idn] = sort(NuA(ndx,n),1,varargin{DrnX});
	% ...keep only relevant indices:
	jdc = ChI(ndx(idc),n); % character
	jdn = NuI(ndx(idn),n); % number
	jde = ~ChI(ndx,n)&~NuI(ndx,n); % empty
	% ...define the sort-order of numbers and characters:
	jdo = any(RsoA)|(~any(RsoB)&C<'0');
	% ...then combine these indices in the requested direction:
	if any(DrnD) % descending
		idx = [idc(jdc&~jdo);idn(jdn);idc(jdc&jdo);ide(jde)];
	else % ascending
		idx = [ide(jde);idc(jdc&jdo);idn(jdn);idc(jdc&~jdo)];
	end
	ndx = ndx(idx);
end
%
ndx  = reshape(ndx,size(X));
X = X(ndx);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsort
