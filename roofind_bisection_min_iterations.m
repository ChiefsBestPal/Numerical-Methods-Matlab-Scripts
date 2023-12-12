function n = roofind_bisection_min_iterations (a0,b0,epsilon)
  % epsilon is max absolute error
  min_n = (log(2.5 - 0.5) - log(epsilon))/log(2) - 1;
  n = ceil(min_n);
endfunction
