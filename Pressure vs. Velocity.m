dl=0; %initialize cumulative length change
for clock=1:clockmax
    t=clock*dt;
    if(use_v(t))
        dx=-v(t)*dt;
    else
        r=(P(t)-Px)/(Px+p1*sum(a));
        dx=(1/mu)*log(1+r);
    end
    x(find(a))=x(find(a))+dx;
    dl=dl+dx;
    pc=(beta*dt)*a+(alpha*dt)*(1-a);
    c=(rand(1,Nb)<pc)|(x>x1);
    a=xor(c,a);
    x(find(a&c))=x0;
    x(find(~a))=0;
    U=sum(a)/Nb;
    Px=sum(p1*(exp(mu*x(find(a)))-1));
end

