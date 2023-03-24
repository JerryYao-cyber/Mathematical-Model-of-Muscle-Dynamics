
Nb=10000; %number of crossbridges(arbitrary)
alpha=14; %(/s)rate of attachment 
beta=126; %(/s)rate of detachment
dt=0.01/(alpha+beta); %(s)time step
tmax=30; %(s)duration of simulation
clockmax=ceil(tmax/dt); %number of time steps
x0=5; %(nm)length of a new crossbridge
x1=10; %(nm)length at which crossbridge must break
p1=4; %(pN)crossbridge force constrant
mu=0.322; %(/nm)multiplyer of x in force
a=zeros(1,Nb); 
x=zeros(1,Nb);
for clock=1:clockmax
    x(find(a))=x(find(a))-v(clock*dt)*dt;
    pc=(beta*dt)*a+(alpha*dt)*(1-a);
    c=(rand(1,Nb)<pc)|(x>x1);
    a=xor(c,a);
    x(find(a&c))=x0;
    x(find(~a))=0;
    U=sum(a)/Nb;
    P=sum(p1*(exp(mu*x)-1));
end

