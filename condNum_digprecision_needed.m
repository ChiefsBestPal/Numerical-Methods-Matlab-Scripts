function sol_num_dig = condNum_digprecision_needed (A,sol_dig_prec_wanted)
  _cond = norm(A,inf)* norm(inv(A),inf); # >= max(Magnifying factor)
  dig_lost = round(log10(_cond));
  sol_num_dig = dig_lost + sol_dig_prec_wanted;
endfunction
