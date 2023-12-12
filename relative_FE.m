function rel_fe = relative_FE (A,b,xr)
  r = A \ b;
  fe = norm(r - xr,inf);
  rel_fe = fe / norm(r,inf);
endfunction
