function plot_pca_reconstruction_error(obj, varargin)

y1 = obj.pca_reconstruction_error_observations;
y2 = obj.pca_reconstruction_error_variables;

obj.plot_reconstruction_error(y1, y2, 'PCA reconstruction error');

end