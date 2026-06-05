function hash = computeFileHash(filepath)
%%  computeFileHash - MD5 hash of a file via Java MessageDigest
%   reads in 8 kB chunks so large .mat files do not load into MATLAB memory.
%   returns a 32-character lowercase hex string. ::
    import java.security.MessageDigest
    md  = MessageDigest.getInstance('MD5');
    fid = fopen(filepath, 'rb');
    while ~feof(fid)
        chunk = fread(fid, 8192, '*uint8');
        if ~isempty(chunk)
            md.update(chunk);
        end
    end
    fclose(fid);
    hash = lower(sprintf('%02x', mod(double(md.digest()), 256)));
end
