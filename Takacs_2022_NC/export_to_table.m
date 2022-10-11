%% This function exports analysis results from the cellList into a table
function res = export_to_table(cellList)

% Get the total number of cells
n = sum(cellfun(@length,cellList.meshData));

% Allocate space for all the results
cell_length_in_um = zeros(n,1);
average_fluor_che = zeros(n,1);
average_fluor_gfp = zeros(n,1);
number_of_spots_che = zeros(n,1);
number_of_spots_gfp = zeros(n,1);
spot_density_per_10um_che = zeros(n,1);
spot_density_per_10um_gfp = zeros(n,1);
mean_spot_intensity_che = zeros(n,1);
mean_spot_intensity_gfp = zeros(n,1);
cv_spot_intensity_che = zeros(n,1);
cv_spot_intensity_gfp = zeros(n,1);
mean_spot_fluor_density_che = zeros(n,1);
mean_spot_fluor_density_gfp = zeros(n,1);
mean_nonspot_fluor_density_che = zeros(n,1);
mean_nonspot_fluor_density_gfp = zeros(n,1);
mean_spot_nonspot_fluor_density_diff_che = zeros(n,1);
mean_spot_nonspot_fluor_density_diff_gfp = zeros(n,1);

% Enumeration index
ii = 1;
for frame = 1:length(cellList.meshData)
    for cc = 1:length(cellList.meshData{frame})
        cs = cellList.meshData{frame}{cc};
        if isfield(cs,'cell_length_in_um')
            cell_length_in_um(ii) = cs.cell_length_in_um;
        end
        if isfield(cs,'che')
            average_fluor_che(ii) = cs.che.average_fluor;
            if isfield(cs.che,'spotPosition')
                number_of_spots_che(ii) = cs.che.number_of_spots;
                spot_density_per_10um_che(ii) = cs.che.spot_density_per_10um;
                mean_spot_intensity_che(ii) = cs.che.mean_spot_intensity;
                cv_spot_intensity_che(ii) = cs.che.cv_spot_intensity;
                mean_spot_fluor_density_che(ii) = cs.che.mean_spot_fluor_density;
                mean_nonspot_fluor_density_che(ii) = cs.che.mean_nonspot_fluor_density;
                mean_spot_nonspot_fluor_density_diff_che(ii) = cs.che.mean_spot_nonspot_fluor_density_diff;
            end
        end
        if isfield(cs,'gfp')
            average_fluor_gfp(ii) = cs.gfp.average_fluor;
            if isfield(cs.gfp,'spotPosition')
                number_of_spots_gfp(ii) = cs.gfp.number_of_spots;
                spot_density_per_10um_gfp(ii) = cs.gfp.spot_density_per_10um;
                mean_spot_intensity_gfp(ii) = cs.gfp.mean_spot_intensity;
                cv_spot_intensity_gfp(ii) = cs.gfp.cv_spot_intensity;
                mean_spot_fluor_density_gfp(ii) = cs.gfp.mean_spot_fluor_density;
                mean_nonspot_fluor_density_gfp(ii) = cs.gfp.mean_nonspot_fluor_density;
                mean_spot_nonspot_fluor_density_diff_gfp(ii) = cs.gfp.mean_spot_nonspot_fluor_density_diff;
            end
        end
        ii = ii + 1;
    end
    
end

% Put all the results into a table
res = table(cell_length_in_um,...
    average_fluor_che,...
    average_fluor_gfp,...
    number_of_spots_che,...
    number_of_spots_gfp,...
    spot_density_per_10um_che,...
    spot_density_per_10um_gfp,...
    mean_spot_intensity_che,...
    mean_spot_intensity_gfp,...
    cv_spot_intensity_che,...
    cv_spot_intensity_gfp,...
    mean_spot_fluor_density_che,...
    mean_spot_fluor_density_gfp,...
    mean_nonspot_fluor_density_che,...
    mean_nonspot_fluor_density_gfp,...
    mean_spot_nonspot_fluor_density_diff_che,...
    mean_spot_nonspot_fluor_density_diff_gfp);

% Remove rows with all zeros
m = table2array(res);
res(all(m == 0,2),:) = [];
end