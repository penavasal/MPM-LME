function makeVtkMP(mpD,mpFileName)

%VTK output file generation: material point data
%--------------------------------------------------------------------------
% Author: Mian Xie(after William Coombs)
% Date:   12/01/2022
% Description:
% Function to generate a VTK file containing the material point data.
%
%--------------------------------------------------------------------------
% MAKEVTKMP(mpD,mpFileName)
%--------------------------------------------------------------------------
% Input(s):
% mpD.x      - material point coordinates (nmp,nD)
% mpD.sig_raw- material point stresses (raw) before mitigating kinematic locking
% mpD.sig    - material point stresses (processed) after mitigating kinematic locking
% mpFileName - VTK file name, for example 'mpData1.vtk'  
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

[nmp,nD]=size(mpD.x);

fid=fopen(mpFileName,'wt');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'MATLAB generated vtk file, MX\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i double\n',nmp);
if nD==3
    for i=1:nmp
        fprintf(fid,'%f %f %f \n',mpD.x(i,:));
    end
elseif nD==2
    for i=1:nmp
        fprintf(fid,'%f %f %f\n',mpD.x(i,:),0);
    end
elseif nD==1
    for i=1:nmp
        fprintf(fid,'%f %f %f\n',mpD.x(i),0,0);
    end
end
fprintf(fid,'\n');

fprintf(fid,'POINT_DATA %i\n',nmp);

fprintf(fid,'SCALARS sigma_xx FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',mpD.sig(i,1));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_yy FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',mpD.sig(i,2));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_zz FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',mpD.sig(i,4));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_xy FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',mpD.sig(i,3));
end
fprintf(fid,'\n');


fprintf(fid,'SCALARS sigma_xx_raw FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',mpD.sig_raw(i,1));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_yy_raw FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',mpD.sig_raw(i,2));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_zz_raw FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',mpD.sig_raw(i,4));
end
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_xy_raw FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:nmp
    fprintf(fid,'%f\n',mpD.sig_raw(i,3));
end
fprintf(fid,'\n');



fclose('all');
end