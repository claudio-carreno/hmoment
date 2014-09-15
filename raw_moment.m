function valor = raw_moment(i,j,fila,columna,matriz)
acum = 0;
for x=1:fila
    for y=1:columna
        acum = acum + (x^i)*(y^j)*matriz(x,y);
    end
end
valor = acum;