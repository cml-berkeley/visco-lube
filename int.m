function int = int(f)

int = 0;
interior = f(2:end-1,2:end-1);

int = int + sum(interior(:));
int = int + 0.5*(sum(f(2:end-1,1)) + sum(f(2:end-1,end)) + sum(f(1,2:end-1)) + sum(f(end,2:end-1)));
int = int + 0.25*(f(1,1)+f(1,end)+f(end,1)+f(end,end));
 
end