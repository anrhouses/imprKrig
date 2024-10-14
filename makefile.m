% DESCRIPTION
% - MAKEFILE that compiles all the MEX function in the util folder
%
% INPUT
% - TYPE:
% 'all': builds all the files in the folder
% 'precise': builds only precise kriging
% 'clean': remove all mex compiled files from the folder
%
function makefile( TYPE )

WIND = 1;
UNI = 2;
MAC = 3;

switch computer
    case 'PCWIN'
        mySys = WIND ;
        arch = 32 ;
    case 'PCWIN64'
        mySys = WIND ;
        arch = 64 ;
    case 'GLNX86'
        mySys = UNI ;
        arch = 32 ;
    case 'GLNXA64'
        mySys = UNI ;
        arch = 64 ;
    otherwise
        mySys = MAC ;
        if (strfind(computer,64)>0)
            arch = 64 ;
        else
            arch = 32 ;
        end
end

% if nothing is specified make all
if ~exist( 'TYPE','var' )
    TYPE = 'all';
end

function compile(TYPE)
    if (mySys == WIND)
        if strcmp( TYPE, 'vargeo' )
            if(arch==64)
                mex -c -v -L.\lib64\win -lgsl -lcblas -I.\include .\src\variogram.cpp .\src\geostat.cpp ;
            else
                mex -c -v -L.\lib32\win -lgsl -lcblas -I.\include .\src\variogram.cpp .\src\geostat.cpp ;
            end
        elseif strcmp( TYPE, 'precise' )
            if(arch==64)
                mex -v -L.\lib64\win -lgsl -lcblas -I.\include .\src\precise_kriging.cpp variogram.obj geostat.obj;
            else
                mex -v -L.\lib32\win -lgsl -lcblas -I.\include .\src\precise_kriging.cpp variogram.obj geostat.obj;
            end
        elseif strcmp( TYPE, 'intervallist' )
            if(arch==64)
                mex -v -L.\lib64\win -lgsl -lcblas -I.\include .\src\intervallist_kriging.cpp variogram.obj geostat.obj;
            else
                mex -v -L.\lib32\win -lgsl -lcblas -I.\include .\src\intervallist_kriging.cpp variogram.obj geostat.obj;
            end
        end
    elseif (mySys == UNI)
      % TODO FOR LINUX
    else
      % TODO FOR MAC
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE ALL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp( TYPE, 'all' ) 
   compile('vargeo');
   compile('precise');
   compile('intervallist');
   delete('.\*.obj');
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE CLEAN %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp( TYPE, 'clean')   
    % clear screen
    clc
    disp('removing all compiled obj mex files.');
    
    % TODO better handle of the specific architectures...
    delete('*.mexw64');
    delete('*.mexw32');
    delete('*.mex');
    delete('.\*.obj');
    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE PRECISE OR INTERVALLIST %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    compile('vargeo');
    compile(TYPE);
    delete('.\*.obj');
end


disp('done!');

end