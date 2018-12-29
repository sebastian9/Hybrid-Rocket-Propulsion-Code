function [ normalizedVector ] = normalize( vector )
    % normalize vector
    if (norm(vector) == 0 )
        normalizedVector = vector;
    else
        normalizedVector = vector/norm(vector);
    end
end