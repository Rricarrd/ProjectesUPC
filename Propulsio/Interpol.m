function y = Interpol(x1,y1,x2,y2,x)
%%Funció per interpolar
y = ((y2-y1)/(x2-x1))*(x-x1)+y1;
end

