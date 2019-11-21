function plot_pca_reconstruction_error(obj, varargin)

y1 = obj.local_pca_reconstruction_error_observations;
y2 = obj.local_pca_reconstruction_error_variables;

obj.plot_reconstruction_error(y1, y2, 'Local PCA reconstruction error');

end