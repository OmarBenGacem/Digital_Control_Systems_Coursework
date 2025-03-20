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
        mattex(eq, fid)
        % fprintf(fid, mattex(eq));

    else

        eqLatex = latex(eq);
        eqLatex = strrep(eqLatex, '\,', '');

        fprintf(fid, '\\n%s\n\\n', eqLatex);
    
    end

    fclose(fid);


end

function [] = mattex(M, fid, frmt)
    
    %% https://uk.mathworks.com/matlabcentral/fileexchange/135251-matrix-to-latex

    % MAT2LATEX(M) displays the latex code for the numerical 2D matrix M in the
    % command window.
    %
    % MAT2LATEX(M, frmt) displays the latex code for the numerical 2D matrix M 
    % in the command window, using frmt, e.g. '%0.2E', as a format specifier 
    % for the matrix content. Default value of frmt is '%0.4f'.
    %
    % Created on Sep. 2023, by Ahmed Mahfouz, in Luxembourg.
    if ~exist("frmt", 'var')
        frmt = '%0.4f';
    end

    % Validate format specifier
    matlabFormatPattern = '^%[-\+]?[0-9]*\.?[0-9]*[bcdeEfFgGosuxX]$';
    if isempty(regexp(frmt, matlabFormatPattern, 'once'))
        error('Invalid format specifier. Please use a valid MATLAB format specifier.\n');
    end
    
    % if class(M) == 'symfun'
    if isa(M, 'symfun')
        M = formula(M);
    end

    
    [nrows, ncols] = size(M);
    fprintf(fid, '\\begin{bmatrix}\n');

    for k = 1:nrows
        for l = 1:ncols
            if isa(M(k,1), 'sym')
                fprintf(fid, char(M(k,1)));
            else
                fprintf(fid, frmt, M(k, l)); 
            end
            if l < ncols
                fprintf(fid, ' & '); 
            end
        end
        if k < nrows
            fprintf(fid, ' \\\\\n');
        else
            fprintf(fid, '\n');
        end
    end

    fprintf(fid, '\\end{bmatrix}\n');

end

