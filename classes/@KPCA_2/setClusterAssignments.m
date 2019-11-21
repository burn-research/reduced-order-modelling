function setClusterAssignments(obj, nc)

i = 1;
while(obj.number_of_clusters_stored(i) ~= nc)
    i = i + 1; 
end

obj.local_idx = obj.local_idx_stored{i};

obj.runLpca();
obj.runLpcaKriging();

end





