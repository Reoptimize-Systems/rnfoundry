function matpoly = polymodel2matpoly(polymodel)
% converts a polymodel structure as produced by polyfitn to a single matrix
% format in which each row corresponds to a term in the model and the first
% column is the coefficient of the term with the following columns the
% power to which the variable is raised in that term

    matpoly = [polymodel.Coefficients' polymodel.ModelTerms];

end