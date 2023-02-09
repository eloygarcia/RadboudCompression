function [X, meta] = nrrdRead(filename)
%NRRDREAD  Import NRRD imagery and metadata.
%   [X, META] = NRRDREAD(FILENAME) reads the image volume and associated
%   metadata from the NRRD-format file specified by FILENAME.
%
%   Current limitations/caveats:
%   * "Block" datatype is not supported.
%   * Only tested with "gzip" and "raw" file encodings.
%   * Very limited testing on actual files.
%   * I only spent a couple minutes reading the NRRD spec.
%
%   See the format specification online:
%   http://teem.sourceforge.net/nrrd/format.html

% Copyright 2012 The MathWorks, Inc.


% Open file.
fid = fopen(filename, 'rb');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

% Magic line.
theLine = fgetl(fid);
assert(numel(theLine) >= 4, 'Bad signature in file.')
assert(isequal(theLine(1:4), 'NRRD'), 'Bad signature in file.')

% The general format of a NRRD file (with attached header) is:
% 
%     NRRD000X
%     <field>: <desc>
%     <field>: <desc>
%     # <comment>
%     ...
%     <field>: <desc>
%     <key>:=<value>
%     <key>:=<value>
%     <key>:=<value>
%     # <comment>
% 
%     <data><data><data><data><data><data>...

meta = struct([]);

% Parse the file a line at a time.
while (true)

  theLine = fgetl(fid);
  
  if (isempty(theLine) || feof(fid))
    % End of the header.
    break;
  end
  
  if (isequal(theLine(1), '#'))
      % Comment line.
      continue;
  end
  
  % "fieldname:= value" or "fieldname: value" or "fieldname:value"
  parsedLine = regexp(theLine, ':=?\s*', 'split','once');
  
  assert(numel(parsedLine) == 2, 'Parsing error')
  
  field = lower(parsedLine{1});
  value = parsedLine{2};
  
  field(isspace(field)) = '';
  meta(1).(field) = value;
  
end

datatype = getDatatype(meta.type);

if(~isfield(meta, 'endian'))
    meta.endian='little';
end


% Get the size of the data.
assert(isfield(meta, 'sizes') && ...
       isfield(meta, 'dimension') && ...
       isfield(meta, 'encoding') && ...
       isfield(meta, 'endian'), ...
       'Missing required metadata fields.')

dims = sscanf(meta.sizes, '%d');

ndims = sscanf(meta.dimension, '%d');
meta.dimension = ndims;

spacedims = ndims;
if(isfield(meta, 'spacedimension'))
    spacedims = sscanf(meta.spacedimension, '%d');
    meta.spacedimension = spacedims;
end


% assert(numel(dims) == spacedims);

data = readData(fid, meta, datatype);
data = adjustEndian(data, meta);

% Reshape and get into MATLAB's order.
% Reshape the spaceorigin and sizes into the meta data too.
X = reshape(data, dims');
if(spacedims >3)
    X = permute(X, [2 1 3 4]);
%     sp = sscanf(meta.spacedirections,'(%f,%f,%f,%f) (%f,%f,%f,%f) (%f,%f,%f,%f) (%f,%f,%f,%f)');            
%     x = sp(5:8);  y = sp(1:4); z =sp(9:12); t=sp(13:16);
%     x= x([2 1 3:end]); y= y([2 1 3:end]); z= z([2 1 3:end]); t= t([2 1 3:end]);
%     meta.spacedirections = [x y z t];
    orig = sscanf(meta.spaceorigin,'(%f,%f,%f,%f)');
    meta.spaceorigin = orig([2 1 3:end]);
    sizes = sscanf(meta.sizes,'%f %f %f %f');
    meta.sizes = sizes([2 1 3:end]);
    meta.space = 'posterior-left-superior-time';
else
    if(spacedims >2)
        X = permute(X, [2 1 3]);

    %     sp = sscanf(meta.spacedirections,'(%f,%f,%f) (%f,%f,%f) (%f,%f,%f)');            
    %     x = sp(4:6);  y = sp(1:3); z =sp(7:9);
    %     x= x([2 1 3]); y= y([2 1 3]); z= z([2 1 3]);
    %     meta.spacedirections = [x y z];
        orig = sscanf(meta.spaceorigin,'(%f,%f,%f)');
        meta.spaceorigin = orig([2 1 3]);
        sizes = sscanf(meta.sizes,'%f %f %f');
        meta.sizes = sizes([2 1 3]);
        meta.space = 'posterior-left-superior';
    else
        if(ndims == 3)% vector type
            X = permute(X, [3 2 1]);
            sizes = sscanf(meta.sizes,'%f %f %f');
            meta.sizes = sizes([3 2 1]);
        else
            X = permute(X, [2 1]);
            sizes = sscanf(meta.sizes,'%f %f');
            meta.sizes = sizes([2 1]);
        end
        
        orig = sscanf(meta.spaceorigin,'(%f,%f)');
        meta.spaceorigin = orig([2 1]);
        
%         meta.spacedimension=spacedims;
    end
        
end









function datatype = getDatatype(metaType)

% Determine the datatype
switch (metaType)
 case {'signed char', 'int8', 'int8_t'}
  datatype = 'int8';
  
 case {'uchar', 'unsigned char', 'uint8', 'uint8_t'}
  datatype = 'uint8';

 case {'short', 'short int', 'signed short', 'signed short int', ...
       'int16', 'int16_t'}
  datatype = 'int16';
  
 case {'ushort', 'unsigned short', 'unsigned short int', 'uint16', ...
       'uint16_t'}
  datatype = 'uint16';
  
 case {'int', 'signed int', 'int32', 'int32_t'}
  datatype = 'int32';
  
 case {'uint', 'unsigned int', 'uint32', 'uint32_t'}
  datatype = 'uint32';
  
 case {'longlong', 'long long', 'long long int', 'signed long long', ...
       'signed long long int', 'int64', 'int64_t'}
  datatype = 'int64';
  
 case {'ulonglong', 'unsigned long long', 'unsigned long long int', ...
       'uint64', 'uint64_t'}
  datatype = 'uint64';
  
 case {'float'}
  datatype = 'single';
  
 case {'double'}
  datatype = 'double';
  
 otherwise
  assert(false, 'Unknown datatype')
end



function data = readData(fidIn, meta, datatype)

switch (meta.encoding)
 case {'raw'}
  
  data = fread(fidIn, inf, [datatype '=>' datatype]);
  
 case {'gzip', 'gz'}

  tmpBase = tempname();
  tmpFile = [tmpBase '.gz'];
  fidTmp = fopen(tmpFile, 'wb');
  assert(fidTmp > 3, 'Could not open temporary file for GZIP decompression')
  
  tmp = fread(fidIn, inf, 'uint8=>uint8');
  fwrite(fidTmp, tmp, 'uint8');
  fclose(fidTmp);
  
  gunzip(tmpFile)
  
  fidTmp = fopen(tmpBase, 'rb');
  cleaner = onCleanup(@() fclose(fidTmp));
  
  meta.encoding = 'raw';
  data = readData(fidTmp, meta, datatype);
  
 case {'txt', 'text', 'ascii'}
  
  data = fscanf(fidIn, '%f');
  data = cast(data, datatype);
  
 otherwise
  assert(false, 'Unsupported encoding')
end



function data = adjustEndian(data, meta)

[~,~,endian] = computer();

needToSwap = (isequal(endian, 'B') && isequal(lower(meta.endian), 'little')) || ...
             (isequal(endian, 'L') && isequal(lower(meta.endian), 'big'));
         
if (needToSwap)
    data = swapbytes(data);
end