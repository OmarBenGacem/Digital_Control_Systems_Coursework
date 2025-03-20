function toOverleaf(eq, fileName, varargin)


    % Check input arguments
    if nargin < 2 || nargin > 3
        error('Function accepts either 2 or 3 arguments.');
    end

    if isstring(fileName)
        fileName = char(fileName);
    end
    targetFolder = fullfile(pwd, '..', 'equations');
  
    eqLatex = latex(eq);
    eqLatex = strrep(eqLatex, '\,', '');
    
    fullFileName = fullfile(targetFolder, [fileName, '.tex']);
    fid = fopen(fullFileName, 'w');
    if fid == -1
        error('Could not open or create file %s for writing.', fullFileName);
    end
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

    fclose(fid);


end


