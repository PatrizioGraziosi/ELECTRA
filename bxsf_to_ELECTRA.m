% Written and developed by                                                %
% Patrizio Graziosi, patrizio.graziosi@gmail.com, during the              %  
% Marie Curie - Individual Fellowships  GENESIS - project ID 788465       %
% Generic transport simulator for new generation thermoelectric materials %
% ----------------------------------------------------------------------- %
% This file is distributed under the terms of the GNU                     %
% General Public License. See the file `LICENSE' in  the root directory   %
% of the present distribution.                                            %
% ----------------------------------------------------------------------- %
%                                                                         %
% Please cite the code source when publishing results obtained            %
% using the present code                                                  %
%                                                                         %
% ----------------------------------------------------------------------- %

function bxsf_to_ELECTRA(fileName,material_name,alat,reduce_flag) %codegen
if strcmp(reduce_flag,'y')
    reduce_SOC = 'SOC_not_magn';
else
    reduce_SOC = 'no';
end
% if nk_new == 0
%     bands_interpolation =   'no'; % if you wants numerical bands interpolation, it uses griddata for sparse points
% else
%     bands_interpolation =   'yes'; % if you wants numerical bands interpolation, it uses griddata for sparse points
% end
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


% extraction codes run

% this sub-function takes the data from the bxsf file and compose a 4D matrix
[points_in_axis_kx, points_in_axis_ky , points_in_axis_kz, num_of_bands,...
BandMatrix, Fermi, a, b, c] = Taking_data_from_bxsf(fileName);

% now we compose the matrixes of the coordinates in k-space

Ek = zeros(points_in_axis_kx,points_in_axis_ky,points_in_axis_kz,num_of_bands);
for id_band = 1:num_of_bands %#ok<*FXUP>
    id_k=1;
    for id_x=1:points_in_axis_kx
        for id_y=1:points_in_axis_ky
            for id_z=1:points_in_axis_kz
                Ek(id_x,id_y,id_z,id_band) = BandMatrix(id_k,id_band)-Fermi;                
                id_k=id_k+1;
            end
        end
    end
end

% alat shall be inputted      
    blat = alat;
    clat = alat;


B_matrix = [ a*2*pi/(alat*1e-9) ; b*2*pi/(blat*1e-9); c*2*pi/(clat*1e-9) ] ;

kx_matrix = zeros(points_in_axis_kx, points_in_axis_ky, points_in_axis_kz);
ky_matrix = kx_matrix; kz_matrix = kx_matrix;
    for id_x = (points_in_axis_kx - 1) : -1 : 0
        for id_y = (points_in_axis_ky - 1) : -1 : 0 % 0 : points_in_axis_ky - 1
            for id_z = (points_in_axis_kz - 1) : -1 : 0
                
                k_vector_not_norm = [id_x id_y id_z]*B_matrix; 
                
                kx_matrix(id_x+1,id_y+1,id_z+1) = 1/points_in_axis_kx * k_vector_not_norm(1);
                ky_matrix(id_x+1,id_y+1,id_z+1) = 1/points_in_axis_ky * k_vector_not_norm(2);
                kz_matrix(id_x+1,id_y+1,id_z+1) = 1/points_in_axis_kz * k_vector_not_norm(3);
            end
        end
    end
    
% clearvars -except Ek kx_matrix ky_matrix kz_matrix a b c kx_array ky_array kz_array bands_interpolation bands_centering nk_new material_name alat


% additional features
if strcmp(reduce_SOC,'SOC_not_magn')
    Ek_full = Ek;
    Ek=zeros(points_in_axis_kx,points_in_axis_ky,points_in_axis_kz,num_of_bands/2);
    for i = 1:num_of_bands/2
        Ek(:,:,:,i) = Ek_full(:,:,:,2*i-1);
    end    
end
% if strcmp(bands_interpolation,'yes')
%     [ Ek, kx_matrix, ky_matrix, kz_matrix ] = Interpolation_3D(Ek,kx_matrix,ky_matrix,kz_matrix,a,b,c,alat,nk_new) ;
% end


save_filename=['Ek_',material_name,'.mat'];   
save(save_filename, 'Ek', 'kx_matrix', 'ky_matrix', 'kz_matrix', 'a', 'b', 'c', 'alat', 'save_filename') 


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [points_in_axis_kx, points_in_axis_ky , points_in_axis_kz, num_of_bands,...
BandMatrix, Fermi, a, b, c] = Taking_data_from_bxsf(fileName)

inputFile = strcat(fileName);
fid = fopen(inputFile);
% searching for the point where the bands start
for i = 1:60
    temp = fgetl(fid);
    if size(temp,2) > 5 && strcmp(temp(1:5),'BAND:')
        num_lines_to_skip = i;
    end
end
temp = fgetl(fid);
first_eigenvalues = str2num(temp);
num_elements_each_row = size(first_eigenvalues,2) ;

% clearvars -except inputFile num_lines_to_skip num_elements_each_row

fid = fopen(inputFile);
for i=1:(num_lines_to_skip-12)
    fgetl(fid);
end
temp = fgetl(fid);
numstr = regexp(temp,'(-)?\d+(\.\d+)?(e(-|+)\d+)?','match') ;
Fermi = str2double(numstr) ; 

% skip the lines
for i=1:4
    fgetl(fid);
end

% number of bands in the file
temp = fgetl(fid);
num_of_bands = str2num(temp) ; 

% number of points in the k space, that corresponds to the number of values in each band 
temp = fgetl(fid);
points_array = str2num(temp) ; 
points_in_axis_kx = points_array(1);
points_in_axis_ky = points_array(2);
points_in_axis_kz = points_array(3);
num_of_points = points_in_axis_kx * points_in_axis_ky * points_in_axis_kz  ; 


% the a, b, c vectors that ate in from 3 to 1 lines above the number of
% line to skip, which indiactes the origin of the grid, usually Gamma
fgetl(fid);
% vectors
temp = fgetl(fid); 
a = str2num(temp); %#ok<*NASGU>
temp = fgetl(fid); 
b = str2num(temp); 
temp = fgetl(fid); 
c = str2num(temp); 


fgetl(fid); % thus we arrive to the num_line_to_skip, extraction starts

if rem(num_of_points,num_elements_each_row) ~= 0
    num_lines_to_read = floor(num_of_points/num_elements_each_row)+1;
elseif rem(num_of_points,num_elements_each_row) == 0
    num_lines_to_read = floor(num_of_points/num_elements_each_row);
end

for id_band = 1:num_of_bands
    if id_band == 1
        starting_line = num_lines_to_skip+1;
    else
        starting_line = num_lines_to_skip+(num_lines_to_read+1)*(id_band-1)+1;
    end
    end_line = starting_line+num_lines_to_read-1;

    
    Band_temp=[];
    for id_line = starting_line:end_line
        
        temp = fgetl(fid);
        
        Band_temp = [Band_temp temp]; %#ok<*AGROW>

    end
    BandMatrix(:,id_band) = str2num(Band_temp); %#ok<*ST2NM>
    num_lines_inbetween = 1; % number of line in between two bands, 1 is the typical values in bsxf files
    % skip the lines
    for i=1:num_lines_inbetween
        temp=fgetl(fid);
    end
    
end

%these are only the axes, regardless of the angles between them
% kx_array = ( -cell_vector_width:(2*cell_vector_width/round(num_of_points_kx-1)):cell_vector_width ) * 2*pi/(alat*1e-9) ; 
% ky_array = ( -cell_vector_width:(2*cell_vector_width/round(num_of_points_ky-1)):cell_vector_width ) * 2*pi/(alat*1e-9) ; 
% kz_array = ( -cell_vector_width:(2*cell_vector_width/round(num_of_points_kz-1)):cell_vector_width ) * 2*pi/(alat*1e-9) ; 
%

% save E_temp % it contains the BandMatrix, one column each band
end

end