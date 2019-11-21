function setPrivateProperty(obj, property_name, value)
evaluation_string = ['obj.', property_name, ' = value;'];
eval(evaluation_string);
end

