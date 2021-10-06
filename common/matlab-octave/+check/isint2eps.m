function result = isint2eps(X)
% isint2eps determines if the numbers in a matrix are integers to the limit
% of the machine floating-point relative accuracy for those values

    theMod = mod(X,1);
    
    result = theMod <= eps(theMod);

end