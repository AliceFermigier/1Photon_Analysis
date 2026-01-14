function red_path = green2red_path(green_path)
    % Split path
    parts = strsplit(green_path, filesep);

    % Replace only last two folders
    parts{end-1} = strrep(parts{end-1}, 'G', 'R');
    parts{end}   = strrep(parts{end},   'G', 'R');

    % Reassemble
    red_path = strjoin(parts, filesep);
end
