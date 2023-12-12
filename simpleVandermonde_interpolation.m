function a_hat = simpleVandermonde_interpolation (xi,yi)
  # A is vandermonde coefficient matrix
  A = flip(vander(xi'),2);

  a_hat = A \ yi';
endfunction
