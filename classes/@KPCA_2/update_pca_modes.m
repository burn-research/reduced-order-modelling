function update_pca_modes(obj)
        obj.modes = obj.modes_init(:, 1:obj.k);
        fprintf('The PCA modes were updated.\n');
end