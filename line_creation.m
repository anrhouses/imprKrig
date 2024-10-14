function li = line_creation( xa,ya,xb,yb,n_samp )

m = (yb-ya)/(xb-xa);
p = ya - m*xa;

li = zeros(n_samp,2);
for i=1:1:n_samp
    li(i,1) = xa + (i-1)*(xb - xa)/(n_samp-1);
    li(i,2) = m*li(i,1)+p;
end

end
