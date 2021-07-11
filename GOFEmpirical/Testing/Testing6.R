# Test functions: gof.logistic; gof.weibull; gof.logistic.bootstrap; gof.weibull.bootstrap
# Ask Richard: gof.uniform.boostrap

x=rlogis(100,5,13)
system.time(gof.logistic(x))
x=rlogis(1000,5,13)
system.time(gof.logistic(x))

x=rlogis(100,5,13)
system.time(gof.logistic.bootstrap(x))
x=rlogis(1000,5,13)
system.time(gof.logistic.bootstrap(x))

x=rweibull(100,5,13)
system.time(gof.weibull(x))
x=rweibull(1000,5,13)
system.time(gof.weibull(x))

x=rweibull(100,5,13)
system.time(gof.weibull.bootstrap(x))
x=rweibull(1000,5,13)
system.time(gof.weibull.bootstrap(x))


