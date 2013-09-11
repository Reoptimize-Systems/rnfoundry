function tf = isodd(X)

     tf = bitget(abs(X),1)~=0;

end