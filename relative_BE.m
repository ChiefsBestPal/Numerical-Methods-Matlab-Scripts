function rel_be = relative_BE (A,b,xr)
  be = norm(A*xr - b,inf);
  rel_be = be / norm(b,inf);
endfunction
