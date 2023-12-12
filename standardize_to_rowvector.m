function u = standardize_to_rowvector(v)
  sizedim = size(v);
  if (!any(sizedim == 1))
    error ("v is not a column or row vector");
   elseif (sizedim(2) == 1)
    u = v';
  else
    u = v;
   endif
endfunction
