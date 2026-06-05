function checkHash(matFilepath)
%%  checkHash - verify a .mat file against its sidecar .hash.json
%   call before load(...) to detect stale intermediates.
%   behavior:
%     no sidecar  → warning  (first run; upstream stage not yet stamped)
%     hash match  → silent pass
%     hash mismatch → error  (file changed since it was last saved; re-run upstream stage) ::
    [folder, name] = fileparts(matFilepath);
    jsonPath = fullfile(folder, [name '.hash.json']);
    if ~isfile(jsonPath)
        warning('checkHash: no hash sidecar for %s — run the upstream stage to stamp it.', matFilepath);
        return
    end
    fid = fopen(jsonPath, 'r');
    s   = jsondecode(fread(fid, '*char')');
    fclose(fid);
    if ~strcmp(computeFileHash(matFilepath), s.hash)
        error('checkHash: %s has changed since it was stamped (hash mismatch). Re-run the upstream stage.', matFilepath);
    end
end
