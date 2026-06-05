function writeHash(matFilepath)
%%  writeHash - compute and store MD5 hash of a .mat file as a sidecar .hash.json
%   call immediately after save(...) to stamp the output.
%   the sidecar lives in the same folder as matFilepath. ::
    [folder, name] = fileparts(matFilepath);
    jsonPath = fullfile(folder, [name '.hash.json']);
    s = struct('filepath', matFilepath, ...
               'hash',     computeFileHash(matFilepath), ...
               'saved',    datestr(now, 'yyyy-mm-ddTHH:MM:SS'));
    fid = fopen(jsonPath, 'w');
    fprintf(fid, '%s\n', jsonencode(s));
    fclose(fid);
end
