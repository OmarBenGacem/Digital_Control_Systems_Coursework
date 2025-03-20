function toOverleaf(eq, fileName, ismatrix, varargin)


    % Check input arguments
    if nargin < 2 || nargin > 3
        error('Function accepts either 2 or 3 arguments.');
    end

    if isstring(fileName)
        fileName = char(fileName);
    end
    targetFolder = fullfile(pwd, '..', 'equations');
  

    
    fullFileName = fullfile(targetFolder, [fileName, '.tex']);
    fid = fopen(fullFileName, 'w');

    if fid == -1
        error('Could not open or create file %s for writing.', fullFileName);
    end


    if ismatrix
        disp(mattex(eq))
        fprintf(fid, mattex(eq));

    else

        eqLatex = latex(eq);
        eqLatex = strrep(eqLatex, '\,', '');

        if nargin == 3
        
            thirdArg = varargin{1};
            if isstring(thirdArg)
                thirdArg = char(thirdArg);
            end
            % fprintf(fid, '\\begin{equation}\n%s = %s\n\\end{equation}\n', thirdArg, eqLatex);
            fprintf(fid, '\\n%s = %s\n\\n', thirdArg, eqLatex);
            % fprintf(fid, 'third_arg = %s\n\n', thirdArg);
        else
            
            % fprintf(fid, '\\begin{equation}\n%s\n\\end{equation}\n', eqLatex);
            fprintf(fid, '\\n%s\n\\n', eqLatex);
    
        end
    end

    fclose(fid);


end

function s = mattex(M)
    %% Convert MATLAB matrix to LaTeX bmatrix format
    [m, n] = size(M);  % Get matrix dimensions

    % Create first line
    s = sprintf('\\begin{bmatrix}\n');  

    % Add matrix content
    for k = 1:m
        for l = 1:n
            s = sprintf('%s %6.3f', s, M(k, l)); % Print 3 decimal places
            if l < n
                s = sprintf('%s &', s); % Add '&' separator between columns
            end
        end
        if k < m
            s = sprintf('%s \\\\\n', s);  % Add double backslash for new row
        else
            % s = sprintf('%s\n', s);
        end
    end

    % Add last line
    s = sprintf('%s\\end{bmatrix}\n', s);
end
