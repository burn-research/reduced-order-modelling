function write(obj, filename, varargin)
%% Description
% Write data to HDF5 file.

%% Notes
% >> fid = H5F.create('myfile.h5');
% >> plist = 'H5P_DEFAULT';
% >> gid = H5G.create(fid,'my_group',plist,plist,plist);
% >> gid2 = H5G.create(fid,'my_group2',plist,plist,plist);
% >> gid1 = H5G.create(gid,'my_group_nested',plist,plist,plist);
% >> H5G.close(gid1);
% >> H5G.close(gid2);
% >> H5G.close(gid);
% >> H5F.close(fid);
% >> h5disp('myfile.h5')
% HDF5 myfile.h5 
% Group '/' 
%   Group '/my_group' 
%       Group '/my_group/my_group_nested' 
%   Group '/my_group2'

% Write a 5-by-5 data set of uint8 values to the root group.
% 
% >> hdf5write('myfile.h5', '/dataset1', uint8(magic(5)))
% Write a 2-by-2 data set of text entries into a subgroup.
% 
% >> dataset = {'north', 'south'; 'east', 'west'};
% >> hdf5write('myfile2.h5', '/group1/dataset1.1', dataset);
% Write a data set and attribute to an existing group.
% 
% >> dset = single(rand(10,10));
% >> dset_details.Location = '/group1/dataset1.2';
% >> dset_details.Name = 'Random';
% 
% >> attr = 'Some random data';
% >> attr_details.Name = 'Description';
% >> attr_details.AttachedTo = '/group1/dataset1.2/Random';
% >> attr_details.AttachType = 'dataset';
% 
% >> hdf5write('myfile2.h5', dset_details, dset, ...
% >>            attr_details, attr, 'WriteMode', 'append');
% Write a data set using objects.
% 
% >> dset = hdf5.h5array(magic(5));
% >> hdf5write('myfile3.h5', '/g1/objects', dset);

% fid = H5F.create('myfile.h5');
%         type_id = H5T.copy('H5T_NATIVE_DOUBLE');
%         dims = [10 5];
%         h5_dims = fliplr(dims);
%         h5_maxdims = h5_dims;
%         space_id = H5S.create_simple(2,h5_dims,h5_maxdims);
%         dcpl = 'H5P_DEFAULT';
%         dset_id = H5D.create(fid,'DS',type_id,space_id,dcpl);
%         H5S.close(space_id);
%         H5T.close(type_id);
%         H5D.close(dset_id);
%         H5F.close(fid);
%         h5disp('myfile.h5');

%% Function
% Create file
filename = [filename, '.h5'];
fid = H5F.create(filename);
fileattrib(filename,'+w');
type_id = H5T.copy('H5T_NATIVE_DOUBLE');

% Add groups
plist = 'H5P_DEFAULT';
gid = cell(length(obj.variable_names) + 1, 1);
for ii = 1 : length(gid) - 1
    % Create group corresponding to the variable
    gid{ii}.root = H5G.create(fid, obj.variable_names{ii}, plist, plist, plist);
    
    % Add data set
    add_dataset_all(obj, true, ii, gid{ii}.root, type_id, 'data', 'TRAINING');
    add_dataset_all(obj, true, ii, gid{ii}.root, type_id, 'pca', 'PCA_RECOVERED');
    add_dataset_all(obj, true, ii, gid{ii}.root, type_id, 'cpca', 'CPCA_RECOVERED');
    add_dataset_all(obj, true, ii, gid{ii}.root, type_id, 'lpca', 'LPCA_RECOVERED');
    add_dataset_all(obj, true, ii, gid{ii}.root, type_id, 'lcpca', 'LCPCA_RECOVERED');
    add_dataset_all(obj, false, ii, gid{ii}.root, type_id, 'data', 'VALIDATION');
    add_dataset_all(obj, false, ii, gid{ii}.root, type_id, 'pca', 'PCA_PREDICTED');
    add_dataset_all(obj, false, ii, gid{ii}.root, type_id, 'cpca', 'CPCA_PREDICTED');
    add_dataset_all(obj, false, ii, gid{ii}.root, type_id, 'lpca', 'LPCA_PREDICTED');
    add_dataset_all(obj, false, ii, gid{ii}.root, type_id, 'lcpca', 'LCPCA_PREDICTED');
    
    % Write to the data set
    write_dataset_all(fid, obj, true, ii, 'data', '/TRAINING');
    write_dataset_all(fid, obj, true, ii, 'pca', '/PCA_RECOVERED');
    write_dataset_all(fid, obj, true, ii, 'cpca', '/CPCA_RECOVERED');
    write_dataset_all(fid, obj, true, ii, 'lpca', '/LPCA_RECOVERED');
    write_dataset_all(fid, obj, true, ii, 'lcpca', '/LCPCA_RECOVERED');
    write_dataset_all(fid, obj, false, ii, 'data', '/VALIDATION');
    write_dataset_all(fid, obj, false, ii, 'pca', '/PCA_PREDICTED');
    write_dataset_all(fid, obj, false, ii, 'cpca', '/CPCA_PREDICTED');
    write_dataset_all(fid, obj, false, ii, 'lpca', '/LPCA_PREDICTED');
    write_dataset_all(fid, obj, false, ii, 'lcpca', '/LCPCA_PREDICTED');
end

% Create a group for info about the ROM
% gid{end}.root = H5G.create(fid, 'ROM', plist, plist, plist);
% add_dataset(obj.pca_approximation_order, gid{end}.root, type_id, 'APPROXIMATION_ORDER');
% add_dataset(obj.number_of_clusters, gid{end}.root, type_id, 'N_CLUSTERS');
% write_dataset(fid, obj.pca_approximation_order, '/APPROXIMATION_ORDER');
% write_dataset(fid, obj.number_of_clusters, '/N_CLUSTERS');

% Close file
H5F.close(fid);


end


% Local functions
function testdata = get_testdata(obj, l, ii, str1)
if l
    testdata = obj.get_variable(obj.training_points, obj.variable_names{ii}, str1);
else
    testdata = obj.get_variable(obj.prediction_points, obj.variable_names{ii}, str1);
end
1;
end

function gid_out = add_dataset(testdata, gid, type_id, str2)    

plist = 'H5P_DEFAULT';
dims = size(testdata);  h5_dims = fliplr(dims); h5_maxdims = h5_dims;
space_id = H5S.create_simple(2, h5_dims, h5_maxdims);
gid_out = H5D.create(gid, str2, type_id, space_id, plist);

end

function gid_out = add_dataset_all(obj, l, ii, gid, type_id, str1, str2)    

plist = 'H5P_DEFAULT';
testdata = get_testdata(obj, l, ii, str1);
dims = size(testdata);  h5_dims = fliplr(dims); h5_maxdims = h5_dims;
space_id = H5S.create_simple(2, h5_dims, h5_maxdims);
gid_out = H5D.create(gid, str2, type_id, space_id, plist);

end

function write_dataset(fid, testdata, dataset_name)

dset_id = H5D.open(fid, dataset_name);
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',testdata);
H5D.close(dset_id);

end

function write_dataset_all(fid, obj, l, ii, str1, str2)

if l
    testdata = obj.get_variable(obj.training_points, obj.variable_names{ii}, str1);
else
    testdata = obj.get_variable(obj.prediction_points, obj.variable_names{ii}, str1);
end
write_dataset(fid, testdata, ['/' obj.variable_names{ii} str2]);

end
