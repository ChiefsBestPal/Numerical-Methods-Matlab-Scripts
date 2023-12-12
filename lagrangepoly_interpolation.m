function f = lagrangepoly_interpolation(xi, yi)

  n = length(xi);
  f = zeros(1, n);
  Ilog = logical(eye(n));
  for i = 1:n
    Pi = poly(xi(~Ilog(i,:)));
    Pi = Pi ./ polyval(Pi, xi(Ilog(i,:)));
    Pi = Pi .* yi(Ilog(i,:));
    f = f + Pi;
   endfor
  endfunction
