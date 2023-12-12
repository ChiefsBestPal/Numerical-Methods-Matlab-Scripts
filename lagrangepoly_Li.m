function Li = lagrangepoly_Li (xi,i,input_x)
  i = i + 1;
  Li = 1;
  for j = 1:length(xi)
    if (j != i)
      Li = Li * (input_x-xi(j))/(xi(i)-xi(j));
    endif
  endfor
endfunction
