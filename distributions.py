import math
from numpy import random as r
import scipy.stats as st
import scipy.special as sp



'''
This package intends to include random number generation for as many statistical distributions as possible. Where possible, the pdf, cdf, median, mode, mean, skewness, kurtosis, variance, standard deviation, and MLE estimation of parameters will be provided.
Any and all corrections and/or additions (as comments) are appreciated.
'''

'''
An asterisk after the name means that I could not find the source I used to obtain the RNG function
'''

class Distribution:
    @staticmethod
    def pdf():
        raise NotImplementedError("Should have implemented pdf function.")
    @staticmethod
    def cdf():
        raise NotImplementedError("Should have implemented cdf function.")
    @staticmethod
    def random():
        raise NotImplementedError("Should have implemented random number generator function.")
    @staticmethod
    def mean():
        raise NotImplementedError("Should have implemented mean function.")
    @staticmethod
    def median():
        raise NotImplementedError("Should have implemented median function.")
    @staticmethod
    def mode():
        raise NotImplementedError("Should have implemented median function.")
    @staticmethod
    def variance():
        raise NotImplementedError("Should have implemented variance function.")
    @staticmethod
    def stddev():
        raise NotImplementedError("Should have implemented standard deviation function.")
    @staticmethod
    def kurtosis():
        raise NotImplementedError("Should have implemented kurtosis function.")
    @staticmethod
    def entropy():
        raise NotImplementedError("Should have implemented entropy function.")
    @staticmethod
    def skewness():
        raise NotImplementedError("Should have implemented skewness function.")
    @staticmethod
    def ppf():
        raise NotImplementedError("Should have implemented ppf function.")


def qexp(x,q): #q-exponential
    if(q!=0):
        return math.pow(1+(1-q)*x,1/(1-q))
    return math.exp(x)

def qlog(x,q): #q-logarithm
    if(x>=0 and q==1):
        return math.log(x)
    if(x>=0 and q!=1):
        return (math.pow(x,1-q)-1)/(1-q)

def rg0():
    n=r.random()
    while(n==0):
        n=r.random()
    return n


'''
Dice roller
'''
def dice(numdice=1,sides=6,expllow=0,explhigh=0,ahigh=0,alow=0):
    sum_=0
    c=0
    d=numdice
    while(c<d):
        a=r.randint(1,sides+1)
        while(expllow!=0 and a<=expllow):
            b=r.randint(1,sides+1)
            a+=b
            if(b>expllow):
                break
        while(explhigh!=0 and a>=explhigh):
            b=r.randint(1,sides+1)
            a+=b
            if(b<explhigh):
                break
        if(ahigh!=0 and a>=ahigh):
            d+=1
        if(alow!=0 and a<=alow):
            d+=1
        sum_+=a
        c+=1
    return sum_


'''
Alpha distribution
'''
class alpha(Distribution):
    @staticmethod
    def random(aa):
        n=rg0()
        return 1/(aa-st.norm.ppf(n*st.norm.cdf(aa)))
    @staticmethod
    def pdf(aa,x):
        return 1/(x**2*st.norm.cdf(aa)*math.sqrt(2*math.pi))*math.exp(-1/2*(aa-1/x)**2)
    @staticmethod
    def cdf(aa,x):
        return st.norm.cdf(aa-1/x)/st.norm.cdf(aa)
    @staticmethod
    def median(aa):
        n=rg0()
        return 1/(aa-st.norm.ppf(st.norm.cdf(aa)/2))
    @staticmethod
    def ppf(aa):
        n=rg0()
        return 1/(aa-st.norm.ppf(q*st.norm.cdf(aa)))

'''
Amoroso distribution
'''
class amoroso(Distribution):
    @staticmethod
    def random(a,theta,aa,bb):
        if(theta>0 and x<a):
            raise InvalidInputError("x must be greater than a")
        if(theta<0 and x>a):
            raise InvalidInputError("x must be less than a")
        return a+theta*math.pow(gamma.random(aa,1),1/bb)
    @staticmethod
    def pdf(a,theta,aa,bb):
        if(theta>0 and x<a):
            raise InvalidInputError("x must be greater than a")
        if(theta<0 and x>a):
            raise InvalidInputError("x must be less than a")
        return 1/math.gamma(aa)*math.abs(bb/theta)*math.pow((x-a)/theta,aa*bb-1)*math.exp(-math.pow((x-a)/theta,bb))
    @staticmethod
    def mean(a,theta,aa,bb):
        if(aa+1/bb>=0):
            return a+theta*math.gamma(aa+1/bb)/math.gamma(aa)
        return None
    @staticmethod
    def mode(a,theta,aa,bb):
        if(aa*bb>=1):
            return a+theta*math.pow(aa-1/bb,1/bb)
        if(aa*bb<=1):
            return a
    @staticmethod
    def variance(a,theta,aa,bb):
        if(aa+2/bb>=0):
            return theta**2*(math.gamma(aa+2/bb)/math.gamma(aa)-(math.gamma(aa+1/bb)/math.gamma(aa))**2)
        return None
    @staticmethod
    def stddev(a,theta,aa,bb):
        if(aa+2/bb>=0):
            return theta*math.sqrt(math.gamma(aa+2/bb)/math.gamma(aa)-(math.gamma(aa+1/bb)/math.gamma(aa))**2)
        return None

'''
Anglit distribution
'''
class anglit(Distribution):
    @staticmethod
    def pdf(x):
        return math.cos(2*x)
    @staticmethod
    def cdf(x):
        return math.sin(x+math.pi/4)**2
    @staticmethod
    def random():
        n=r.random()
        return math.asin(math.sqrt(n))-math.pi/4
    @staticmethod
    def median():
        return math.asin(math.sqrt(1/2))-math.pi/4
    @staticmethod
    def ppf(q):
        return math.asin(math.sqrt(q))-math.pi/4

'''
Arcsine distribution
'''
class arcsin(Distribution):
    @staticmethod
    def pdf(x):
        if(x<0 or x>1):
            raise InvalidInputError("Must be between 0 and 1")
        return 1/(math.pi*math.sqrt(x*(1-x)))
    @staticmethod
    def cdf(x):
        if(x<0 or x>1):
            raise InvalidInputError("Must be between 0 and 1")
        return math.asin(2*x-1)/math.pi+1/2
    @staticmethod
    def random():
        return math.pow(math.sin(math.pi/2*r.random()),2)
    @staticmethod
    def mean():
        return 1/2
    @staticmethod
    def median():
        return 1/2
    @staticmethod
    def variance():
        return 1/8
    @staticmethod
    def stddev():
        return math.sqrt(1/8)
    @staticmethod
    def ppf(q):
        math.pow(math.sin(math.pi/2*q),2)

'''
Arctan distribution
'''
class arctan(Distribution):
    @staticmethod
    def pdf(lmbda,phi,x):
        if(lmbda<0):
            raise InvalidInputError("Must be between 0 and 1")
        return lmbda/((math.atan(lmbda*phi)+1/(2*math.pi))*(1+lmbda**2*((x-phi)**2)))
    @staticmethod
    def cdf(lmbda,phi,x):
        if(lmbda<0):
            raise InvalidInputError("Must be between 0 and 1")
        return 2*((math.atan(lmbda*phi)-math.atan(-x*lmbda+lmbda*phi))/(2*math.atan(lmbda*phi)+math.pi))
    @staticmethod
    def random(lmbda,phi):
        if(lmbda<0):
            raise InvalidInputError("Must be between 0 and 1")
        n=r.random()
        return (lmbda*phi+math.tan(-math.atan(lmbda*phi)+n*math.atan(lmbda*phi)+n*math.pi/2))/lmbda
    @staticmethod
    def kurtosis(lmbda,phi):
        if(lmbda<0):
            raise InvalidInputError("Must be between 0 and 1")
        return None
    @staticmethod
    def mean(lmbda,phi):
        if(lmbda<0):
            raise InvalidInputError("Must be between 0 and 1")
        return None
    @staticmethod
    def variance(lmbda,phi):
        if(lmbda<0):
            raise InvalidInputError("Must be between 0 and 1")
        return float("infinity")
    @staticmethod
    def stddev(lmbda,phi):
        if(lmbda<0):
            raise InvalidInputError("Must be between 0 and 1")
        return float("infinity")
    @staticmethod
    def ppf(lmbda,phi,q):
        return (lmbda*phi+math.tan(-math.atan(lmbda*phi)+q*math.atan(lmbda*phi)+q*math.pi/2))/lmbda

'''
Asymmetric Laplace distribution
'''
class asymlaplace(Distribution): #asymmetric laplace
    @staticmethod
    def pdf(m,lmbda,kppa,x):
        if(kppa<=0 or lmbda<=0):
            raise InvalidInputError("kppa and lambda must be greater than 0")
        if(x>=m):
            return lmbda/(kppa+1/kppa)*math.exp(-lmbda*kppa*(x-m))
        return lmbda/(kppa+1/kppa)*math.exp((lmbda/kppa)*(x-m))
    @staticmethod
    def cdf(m,lmbda,kppa,x):
        if(kppa<=0 or lmbda<=0):
            raise InvalidInputError("kppa and lambda must be greater than 0")
        if(x>m):
            return 1-1/(1+kppa**2)*math.pow(-lmbda*kppa*(x-m))
        return kppa**2/(1+kppa**2)*math.pow(lmbda/kppa*(x-m))
    @staticmethod
    def random(m,kppa,lmbda):
        if(kppa<=0 or lmbda<=0):
            raise InvalidInputError("kppa and lambda must be greater than 0")
        n=rg0()
        while(n==1):
            n=rg0()
        if(n<=(math.pow(k,2))/(1+math.pow(k,2))):
            return k/math.sqrt(2)*math.log((1+math.pow(k,2))/math.pow(k,2) * n)
        return 1/(k*math.sqrt(2))*math.log((1+math.pow(k,2))*(1-n))
    @staticmethod
    def mean(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise InvalidInputError("kppa and lambda must be greater than 0")
        return m+(1-kppa**2)/(lmbda*kppa)
    @staticmethod
    def median(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise InvalidInputError("kppa and lambda must be greater than 0")
        return m+(kppa/lmbda)*math.log((1+kppa**2)/(2*kppa**2))
    @staticmethod
    def variance(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise InvalidInputError("kppa and lambda must be greater than 0")
        return (1+math.pow(kppa,4))/(lmbda**2*kppa**2)
    @staticmethod
    def stddev(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise InvalidInputError("kppa and lambda must be greater than 0")
        return math.sqrt((1+math.pow(kppa,4))/(lmbda**2*kppa**2))
    @staticmethod
    def skewness(m,lmbda,kppa):
        if(kppa<=0 or lmbda<=0):
            raise InvalidInputError("kppa and lambda must be greater than 0")
        return 2*(1-math.pow(kppa,6))/(math.pow(math.pow(kppa,4)+1,3/2))
    @staticmethod
    def entropy(m,kppa,lmbda):
        if(kppa<=0 or lmbda<=0):
            raise InvalidInputError("kppa and lambda must be greater than 0")
        return math.log(math.exp(1)*(1+kppa**2)/(kppa*lmbda))
    @staticmethod
    def ppf(m,kppa,lmbda,q):
        if(kppa<=0 or lmbda<=0):
            raise InvalidInputError("kppa and lambda must be greater than 0")
        if(q<=(math.pow(k,2))/(1+math.pow(k,2))):
            return k/math.sqrt(2)*math.log((1+math.pow(k,2))/math.pow(k,2)*q)
        return 1/(k*math.sqrt(2))*math.log((1+math.pow(k,2))*(1-q))

'''
Balding-Nichols distribution
'''
class baldingnichols(Distribution):
    @staticmethod
    def random(f,p):
        if(F<=0 or F>=1 or p<=0 or p>=1):
            raise InvalidInputError("all inputs must be between 0 and 1 exclusive")
        return beta.random(p*(1-f)/f,(1-p)*(1-f)/f)
    @staticmethod
    def mean(f,p):
        if(F<=0 or F>=1 or p<=0 or p>=1):
            raise InvalidInputError("all inputs must be between 0 and 1 exclusive")
        return p
    @staticmethod
    def mode(f,p):
        if(F<=0 or F>=1 or p<=0 or p>=1):
            raise InvalidInputError("all inputs must be between 0 and 1 exclusive")
        return (F-(1-F)*p)/(3*F-1)
    @staticmethod
    def variance(f,p):
        if(F<=0 or F>=1 or p<=0 or p>=1):
            raise InvalidInputError("all inputs must be between 0 and 1 exclusive")
        return F*p*(1-p)
    @staticmethod
    def stddev(f,p):
        if(F<=0 or F>=1 or p<=0 or p>=1):
            raise InvalidInputError("all inputs must be between 0 and 1 exclusive")
        return math.sqrt(F*p*(1-p))
    @staticmethod
    def skewness(f,p):
        if(F<=0 or F>=1 or p<=0 or p>=1):
            raise InvalidInputError("all inputs must be between 0 and 1 exclusive")
        return 2*F*(1-2*p)/((1+F)*math.sqrt(F*(1-p)*p))

'''
Banach distribution*
'''
class banach(Distribution):
    @staticmethod
    def random(x,n):
        n=rg0()
        sum_=0
        counter=0
        while(sum_<n):
            sum_+=math.pow(2,x-2*n)*(math.factorial(2*n-x)/(math.factorial(n)*math.factorial(n-x)))
            counter+=1
        return counter

'''
Bartlett distribution*
'''
class bartlett(Distribution):
    @staticmethod
    def random(a,q):
        return poisson.random(a)*geometric.random(q)

'''
Bates distribution
'''
class bates(Distribution):
    @staticmethod
    def random():
        if(n%1!=0 or n<1):
            raise InvalidInputError("n must be a positive integer greater than 0")
        avg=0
        for _ in range(n):
            avg+=r.random()
        return avg/n
    @staticmethod
    def kurtosis(n):
        if(n%1!=0 or n<1):
            raise InvalidInputError("n must be a positive integer greater than 0")
        return -6/(5*n)
    @staticmethod
    def mean(n):
        if(n%1!=0 or n<1):
            raise InvalidInputError("n must be a positive integer greater than 0")
        return 1/2
    @staticmethod
    def variance(n):
        if(n%1!=0 or n<1):
            raise InvalidInputError("n must be a positive integer greater than 0")
        return 1/(12*n)
    @staticmethod
    def stddev(n):
        if(n%1!=0 or n<1):
            raise InvalidInputError("n must be a positive integer greater than 0")
        return math.sqrt(1/(12*n))
    @staticmethod
    def skewness(n):
        if(n%1!=0 or n<1):
            raise InvalidInputError("n must be a positive integer greater than 0")
        return 0

'''
Behrens-Fisher distribution
'''
class behrensfisher(Distribution):
    @staticmethod
    def random(n1,n2,theta):
        return st.t.rvs(n2)*math.cos(theta)-st.t.rvs(n1)*math.sin(theta)

'''
Benford distribution*
'''
class benford(Distribution):
    @staticmethod
    def random(b):
        return 1/(math.pow(b,1-r.random()))

'''
Benini distribution
'''
class benini(Distribution):
    @staticmethod
    def pdf(aa,bb,sigma):
        if(aa<=0 or bb<=0 or sigma<=0 or x<=sigma):
            raise InvalidInputError("aa, bb, and gamma must be greater than 0")
        return math.exp(-aa*math.log(x/sigma)-bb*(math.log(x/sigma))**2)*(aa/x+2*bb*math.log(x/sigma)/x)
    @staticmethod
    def cdf(aa,bb,sigma):
        if(aa<=0 or bb<=0 or sigma<=0 or x<=sigma):
            raise InvalidInputError("aa, bb, and gamma must be greater than 0")
        return 1-math.exp(-aa*math.log(x/sigma)-bb*(math.log(x/sigma))**2)
    @staticmethod
    def median(aa,bb,sigma):
        if(aa<=0 or bb<=0 or sigma<=0 or x<=sigma):
            raise InvalidInputError("aa, bb, and gamma must be greater than 0")
        return sigma*math.exp((-aa+math.sqrt(aa**2+bb*math.log(16)))/(2*bb))

'''
Benktander Weibull distribution
'''
class benktanderweibull(Distribution):
    @staticmethod
    def pdf(a,b,x):
        if(b<=0 or b>1 or x<1):
            raise InvalidInputError("b must be greater than 0 and less than or equal to 1, and x must be greater than or equal to 1")
        return math.exp(a*(1-math.pow(x,b))/b)*math.pow(x,b-2)*(1-b+a*math.pow(x,b))
    @staticmethod
    def cdf(a,b,x):
        if(b<=0 or b>1 or x<1):
            raise InvalidInputError("b must be greater than 0 and less than or equal to 1, and x must be greater than or equal to 1")
        return 1-math.exp(a*(1-math.pow(x,b))/b)*math.pow(x,b-1)
    @staticmethod
    def random(a,b):
        if(b<=0 or b>1):
            raise InvalidInputError("b must be greater than 0 and less than or equal to 1")
        n=r.random()
        return math.pow((b-1)*sp.lambertw(-(a*math.exp(-a/(b+1))*math.pow(1-n,b/(b+1)))/(b-1))/a,1/b)
    @staticmethod
    def mean(a,b):
        if(b<=0 or b>1):
            raise InvalidInputError("b must be greater than 0 and less than or equal to 1")
        return 1+1/a
    @staticmethod
    def ppf(a,b,q):
        if(b<=0 or b>1):
            raise InvalidInputError("b must be greater than 0 and less than or equal to 1")
        return math.pow((b-1)*sp.lambertw(-(a*math.exp(-a/(b+1))*math.pow(1-q,b/(b+1)))/(b-1))/a,1/b)

'''
Bernoulli distribution
'''
class bernoulli(Distribution):
    @staticmethod
    def pdf(k,p):
        if((k!=0 and k!=1) or p<0 or p>1):
            raise InvalidInputError("k must be 0 or 1, and p must be between 0 and 1 inclusive")
        if(k==0):
            return 1-p
        return p
    @staticmethod
    def cdf(k,p):
        if((k!=0 and k!=1) or p<0 or p>1):
            raise InvalidInputError("k must be 0 or 1, and p must be between 0 and 1 inclusive")
        if(k==0):
            return 1-p
        return 1
    @staticmethod
    def random(p):
        if(p<0 or p>1):
            raise InvalidInputError("p must be between 0 and 1 inclusive")
        if(r.random()>=p):
            return 1
        return 0
    @staticmethod
    def kurtosis(p):
        if(p<0 or p>1):
            raise InvalidInputError("p must be between 0 and 1 inclusive")
        return (1-6*p*(1-p))/(p*(1-p))
    @staticmethod
    def mean(p):
        if((k!=0 and k!=1) or p<0 or p>1):
            raise InvalidInputError("k must be 0 or 1, and p must be between 0 and 1 inclusive")
        return p
    @staticmethod
    def median(p):
        if(p<0 or p>1):
            raise InvalidInputError("p must be between 0 and 1 inclusive")
        if(p<1/2):
            return 0
        if(p==1/2):
            return 0.5
        if(p==1):
            return 1
    @staticmethod
    def variance(p):
        if(p<0 or p>1):
            raise InvalidInputError("p must be between 0 and 1 inclusive")
        return p*(1-p)
    @staticmethod
    def stddev(p):
        if(p<0 or p>1):
            raise InvalidInputError("p must be between 0 and 1 inclusive")
        return math.sqrt(p*(1-p))
    @staticmethod
    def entropy(p):
        if(p<0 or p>1):
            raise InvalidInputError("p must be between 0 and 1 inclusive")
        return 1/(p*(1-p))
    @staticmethod
    def skewness(p):
        if(p<0 or p>1):
            raise InvalidInputError("p must be between 0 and 1 inclusive")
        return (1-2*p)/math.sqrt(p*(1-p))
    @staticmethod
    def ppf(p,q):
        if(p<0 or p>1):
            raise InvalidInputError("p must be between 0 and 1 inclusive")
        if(q>=p):
            return 1
        return 0
    @staticmethod
    def mle(x):
        return {'p':sum(x)/len(x)}

'''
Beta distribution
'''
class beta(Distribution):
    @staticmethod
    def random(aa,bb):
        if(aa<=0 or bb<=0):
            raise InvalidInputError("aa and bb must be bigger than 0")
        x=gamma.random(aa)
        return x/(x+gamma.random(x+y))
    @staticmethod
    def kurtosis(aa,bb):
        if(aa<=0 or bb<=0):
            raise InvalidInputError("aa and bb must be bigger than 0")
        return 6*((aa-bb)**2*(aa+bb+1)-aa*bb*(aa+bb+2))/(aa*bb*(aa+bb+2)*(aa+bb+3))
    @staticmethod
    def mean(aa,bb):
        if(aa<=0 or bb<=0):
            raise InvalidInputError("aa and bb must be bigger than 0")
        return aa/(aa+bb)
    @staticmethod
    def median(aa,bb):
        if(aa<=0 or bb<=0):
            raise InvalidInputError("aa and bb must be bigger than 0")
    @staticmethod
    def variance(aa,bb):
        if(aa<=0 or bb<=0):
            raise InvalidInputError("aa and bb must be bigger than 0")
        return aa*bb/((aa+bb)**2*(aa+bb+1))
    @staticmethod
    def stddev(aa,bb):
        if(aa<=0 or bb<=0):
            raise InvalidInputError("aa and bb must be bigger than 0")
        return math.sqrt(aa*bb/((aa+bb)**2*(aa+bb+1)))
    @staticmethod
    def skewness(aa,bb):
        if(aa<=0 or bb<=0):
            raise InvalidInputError("aa and bb must be bigger than 0")
        return 2*(bb-aa)*math.sqrt(aa+bb+1)/((aa+bb+2)*math.sqrt(aa*bb))

'''
Beta-binomial distribution
'''
class betabinomial(Distribution):
    @staticmethod
    def random(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise InvalidInputError("aa and bb must be bigger than 0 and n must be a positive integer")
        return binomial.random(beta.random(aa,bb),n)
    @staticmethod
    def kurtosis(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise InvalidInputError("aa and bb must be bigger than 0 and n must be a positive integer")
        u=(aa+bb)**2*(1+aa+bb)/(n*aa*bb*(aa+bb+2)*(aa+bb+3)*(aa+bb+n))
        v=(aa+bb)*(aa+bb-1+6*n)+3*aa*bb*(n-2)+6*n**2-3*aa*bb*n*(6-n)/(aa+bb)-18*aa*bb*n**2/((aa+bb)**2)
        return u*v
    @staticmethod
    def mean(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise InvalidInputError("aa and bb must be bigger than 0 and n must be a positive integer")
        return n*aa/(aa+bb)
    @staticmethod
    def variance(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise InvalidInputError("aa and bb must be bigger than 0 and n must be a positive integer")
        return n*aa*bb*(aa+bb+n)/((aa+bb)**2*(aa+bb+1))
    @staticmethod
    def stddev(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise InvalidInputError("aa and bb must be bigger than 0 and n must be a positive integer")
        return math.sqrt(n*aa*bb*(aa+bb+n)/((aa+bb)**2*(aa+bb+1)))
    @staticmethod
    def skewness(aa,bb,n):
        if(aa<=0 or bb<=0 or n%1!=0 or n<1):
            raise InvalidInputError("aa and bb must be bigger than 0 and n must be a positive integer")
        return (aa+bb+2*n)*(bb-aa)/(aa+bb+2)*math.sqrt((1+aa+bb)/(n*aa*bb)*(n+aa+bb))

'''
Beta-exponential Burr XII distribution*
'''
class betaexpburr12(Distribution): #beta-exponential burr XII
    @staticmethod
    def random(a,b,c,d,k):
        u=beta.random(a,b)
        return math.pow(math.pow(1-math.pow(u,1/b),-1/k)-1,1/c)

'''
Beta-exponential distribution
'''
class betaexp(Distribution): #beta-exponential
    @staticmethod
    def random(a,g):
        return math.log(beta.random(a,g))

'''
Beta-Gompertz distribution
'''
class betagompertz(Distribution): 
    @staticmethod
    def random(g,t,a,b):
        return (1/g)*math.log(1-g/t*math.log(1-beta.random(a,b)))
    @staticmethod
    def pdf(g,t,a,b,x):
            return t*math.exp(g*x)*math.exp(-b*t/g*(math.exp(g*x)-1))/sp.beta(a,b)*math.pow(1-math.exp(-t/g*(math.exp(g*x)-1)),a-1)

'''
Beta-Pascal distribution*
'''
class betapascal(Distribution):
    @staticmethod
    def random(k,m,n):
        return negbin.random(k,beta.random(m,n))

'''
Beta prime distribution
'''
class betaprime(Distribution):
    @staticmethod
    def random(aa,bb):
        if(aa<=0 or bb<=0):
            raise InvalidInputError("aa and bb must be greater than 0")
        return gamma.random(aa,1)/gamma.random(bb,1)
    @staticmethod
    def mean(aa,bb):
        if(aa<=0 or bb<=0):
            raise InvalidInputError("aa and bb must be greater than 0")
        if(beta>1):
            return aa/(bb-1)
        return None
    @staticmethod
    def mode(aa,bb):
        if(aa<=0 or bb<=0):
            raise InvalidInputError("aa and bb must be greater than 0")
        if(aa>=1):
            return (aa-1)/(bb+1)
        return 0

'''
Bhattacharya negative binomial distribution*
'''
class bhattacharyanegbin(Distribution): #bhattacharya negative binomial
    @staticmethod
    def random(a,b,k):
        return negbin.random(k+negbin.random(a,b/(b+1)),(b+1)/(b+2))

'''
Binomial distribution
'''
class binomial(Distribution):
    @staticmethod
    def random(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise InvalidInputError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        count=0
        for _ in range(n):
            if r.random()<p:
                count+=1
        return count
    @staticmethod
    def kurtosis(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise InvalidInputError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return (1-6*p*(1-p))/(n*p*(1-p))
    @staticmethod
    def mean(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise InvalidInputError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return n*p
    @staticmethod
    def median(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise InvalidInputError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return int(round(n*p))
    @staticmethod
    def variance(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise InvalidInputError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return n*p*(1-p)
    @staticmethod
    def stddev(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise InvalidInputError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return math.sqrt(n*p*(1-p))
    @staticmethod
    def skewness(n,p):
        if(n%1!=0 or n<0 or p<0 or p>1):
            raise InvalidInputError("n must be a non-negative integer and p must be between 0 and 1 inclusive")
        return (1-2*p)/math.sqrt(n*p*(1-p))

'''
Birnbaum-Saunders distribution
'''
class birnbaumsaunders(Distribution): #standard fatigue-life
    '''Fatigue-life'''
    @staticmethod
    def random(gmma):
        n=r.random()
        return 1/4.0*math.pow(gmma*st.norm.ppf(n)+math.sqrt(4+math.pow(gmma*st.norm.ppf(n),2)),2)
    @staticmethod
    def pdf(gmma,x):
        return (math.sqrt(x)+math.sqrt(1/x))/(2*gmma*x)*st.norm.pdf((math.sqrt(x)-math.sqrt(1/x))/gmma)
    @staticmethod
    def cdf(gmma,x):
        return st.norm.cdf((math.sqrt(x)-math.sqrt(1/x))/gmma)
    @staticmethod
    def mean(gmma):
        return 0
    @staticmethod
    def median(gmma):
        return 1/4.0*math.pow(gmma*st.norm.ppf(1/2)+math.sqrt(4+math.pow(gmma*st.norm.ppf(1/2),2)),2)
    @staticmethod
    def ppf(gmma,q):
        return 1/4.0*math.pow(gmma*st.norm.ppf(q)+math.sqrt(4+math.pow(gmma*st.norm.ppf(q),2)),2)

'''
Boltzmann distribution
'''
class boltzmann(Distribution):
    @staticmethod
    def random(N,l):
        return math.ceil(-math.log(1-rg0()*(1-math.exp(-l*N)))/l-1)
    @staticmethod
    def pdf(N,l,x):
        return (1-math.exp(-l))/(1-math.exp(-l*N))*math.exp(-l*x)
    @staticmethod
    def cdf(N,l,x):
        if(x>=N-1):
            return 1
        return (1-math.exp(-l*(math.floor(x)+1)))/(1-math.exp(-l*N))
    @staticmethod
    def median(N,l):
        return math.ceil(-math.log(1-(1-math.exp(-l*N))/2)/l-1)
    @staticmethod
    def ppf(N,l):
        return math.ceil(-math.log(1-q*(1-math.exp(-l*N)))/l-1)

'''
Bounded pareto distribution
'''
class boundedpareto(Distribution):
    @staticmethod
    def pdf(aa,H,L,x):
        if(aa<=0):
            raise InvalidInputError("aa must be bigger than 0")
        if(L<=0):
            raise InvalidInputError("L must be bigger than 0")
        if(H<=L):
            raise InvalidInputError("H must be bigger than L")
        if(x<H or x>L):
            raise InvalidInputError("x must be between H and L inclusive")
        return (aa*math.pow(L,aa)*math.pow(x,-aa-1))/(1-math.pow(L/H,aa))
    @staticmethod
    def cdf(aa,H,L,x):
        if(aa<=0):
            raise InvalidInputError("aa must be bigger than 0")
        if(L<=0):
            raise InvalidInputError("L must be bigger than 0")
        if(H<=L):
            raise InvalidInputError("H must be bigger than L")
        if(x<H or x>L):
            raise InvalidInputError("x must be between H and L inclusive")
        return (1-math.pow(L,aa)*math.pow(x,-aa))/(1-math.pow(L/H,aa))
    @staticmethod
    def random(aa,H,L):
        if(aa<=0):
            raise InvalidInputError("aa must be bigger than 0")
        if(L<=0):
            raise InvalidInputError("L must be bigger than 0")
        if(H<=L):
            raise InvalidInputError("H must be bigger than L")
        n=r.random()
        return math.pow(-(n*math.pow(H,aa)-n*math.pow(L,aa)-math.pow(H,aa))/(math.pow(H,aa)*math.pow(L,aa)),-1/aa)
    @staticmethod
    def mean(aa,H,L):
        if(aa<=0):
            raise InvalidInputError("aa must be bigger than 0")
        if(L<=0):
            raise InvalidInputError("L must be bigger than 0")
        if(H<=L):
            raise InvalidInputError("H must be bigger than L")
        if(aa!=1):
            return math.pow(L,aa)/(1-math.pow(L/H,aa))*(aa/(aa-1))*(1/math.pow(L,aa-1)-1/math.pow(H,aa-1))
        return None
    @staticmethod
    def median(aa,H,L):
        if(aa<=0):
            raise InvalidInputError("aa must be bigger than 0")
        if(L<=0):
            raise InvalidInputError("L must be bigger than 0")
        if(H<=L):
            raise InvalidInputError("H must be bigger than L")
        return L*math.pow(1-1/2.0*(1-math.pow(L/H,aa)),-1/aa)
    @staticmethod
    def ppf(aa,H,L,q):
        if(aa<=0):
            raise InvalidInputError("aa must be bigger than 0")
        if(L<=0):
            raise InvalidInputError("L must be bigger than 0")
        if(H<=L):
            raise InvalidInputError("H must be bigger than L")
        return math.pow(-(q*math.pow(H,aa)-q*math.pow(L,aa)-math.pow(H,aa))/(math.pow(H,aa)*math.pow(L,aa)),-1/aa)

'''
Bradford distribution*
'''
class bradford(Distribution):
    @staticmethod
    def random(min_,max_,theta):
        n=r.random()
        return (-min_*math.pow(theta+1,n)+theta*min_+max_*math.pow(theta+1,n)+min_-max_)/theta

'''
Bramwell-Holdsworth-Pinton distribution*
'''
class bramwellholdsworthpinton(Distribution):
    @staticmethod
    def random(nu,lmbda):
        return gammaexp.random(nu,lmbda,math.pi/2)

'''
Burr II distribution
'''
class burr2(Distribution):
    @staticmethod
    def pdf(r,mu,sigma,x):
        if(r<=0 or sigma<=0):
            raise InvalidInputError("r and sigma must be greater than 0")
        return 1/sigma*r*math.pow(1+math.exp((-x+mu)/sigma),-r-1)/math.exp((x-mu)/sigma)
    @staticmethod
    def cdf(r,mu,sigma,x):
        if(r<=0 or sigma<=0):
            raise InvalidInputError("r and sigma must be greater than 0")
        return 1/math.pow(1+math.exp((x-mu)/sigma),r)
    @staticmethod
    def random(r,mu,sigma): #standard burr ii
        if(r<=0 or sigma<=0):
            raise InvalidInputError("r and sigma must be bigger than 0")
        return mu+sigma*math.log(math.pow(1/rg0(),1/r)-1)
    @staticmethod
    def median(r,mu,sigma):
        return mu+sigma*math.log(math.pow(2,1/r)-1)
    @staticmethod
    def ppf(r,mu,sigma,q):
        if(r<=0 or sigma<=0):
            raise InvalidInputError("r and sigma must be bigger than 0")
        return mu+sigma*math.log(math.pow(1/q,1/r)-1)

'''
Burr III distribution
'''
class burr3(Distribution):
    @staticmethod
    def random(r,k,l,s):
        return s*math.pow(math.pow(rg0(),-1/r)-1,-1/k)+l
    @staticmethod
    def pdf(r,k,l,s,x):
        return 1/s*r*k*math.pow((x-l)/s,r*k-1)/math.pow(1+math.pow(x,k),r+1)
    @staticmethod
    def cdf(r,k,l,s,x):
        return math.pow(1+math.pow((x-l)/s,-k),-r)
    @staticmethod
    def median(r,k,l,s):
        return math.pow(math.pow(1/2,-1/r)-1,-1/k)+l
    @staticmethod
    def ppf(r,k,l,s,q):
        return s*math.pow(math.pow(q,-1/r)-1,-1/k)+l

'''
Burr V distribution
'''
class burr5(Distribution):
    @staticmethod
    def random(r,k,l,s):
        return s*math.atan(-math.log((math.pow(rg0(),-1/r)-1)/k))+l
    @staticmethod
    def pdf(r,k,l,s,x):
        return r*k*math.pow(1+k/math.exp(math.tan((x-l)/s)),-r-1)*math.sec((x-l)/s)**2/math.exp(math.tan(x))
    @staticmethod
    def cdf(r,k,l,s,x):
        return math.pow(1+k*math.exp(-math.tan((x-l)/s)),-r)
    @staticmethod
    def median(r,k,l,s):
        return s*math.atan(-math.log((math.pow(1/2,-1/r)-1)/k))+l
    @staticmethod
    def ppf(r,k,l,s,q):
        return s*math.atan(-math.log((math.pow(q,-1/r)-1)/k))+l

'''
Burr VI distribution
'''
class burr6(Distribution):
    @staticmethod
    def random(r,k,l,s):
        return s*math.asinh(-math.log((math.pow(rg0(),-1/r)-1)/k))+l
    @staticmethod
    def pdf(r,k,l,s,x):
        return r*k*math.pow(1+k/math.exp(math.sinh(x)),-r-1)*math.cosh(x)/math.exp(math.sinh(x))
    @staticmethod
    def cdf(r,k,l,s,x):
        return math.pow(1+k*math.exp(-math.sinh((x-l)/s)),-r)
    @staticmethod
    def median(r,k,l,s):
        return s*math.asinh(-math.log((math.pow(1/2,-1/r)-1)/k))+l
    @staticmethod
    def ppf(r,k,l,s,q):
        return s*math.asinh(-math.log((math.pow(q,-1/r)-1)/k))+l

'''
Burr VII distribution
'''
class burr7(Distribution):
    @staticmethod
    def random(r):
        return s*math.atanh(math.pow(rg0()*math.pow(2,r),1/r)-1)+l
    @staticmethod
    def pdf(r,l,s,x):
            return 1/s*r*(x**2)*math.pow(1+math.tanh((x-l)/s),r-1)/math.pow(2,r)
    @staticmethod
    def cdf(r,l,s,x):
            return math.pow(2,-r)*math.pow(1+math.tanh((x-l)/s),r)
    @staticmethod
    def median(r,l,s):
        return s*math.atanh(math.pow(math.pow(2,r)/2,1/r)-1)+l
    @staticmethod
    def ppf(r,l,s,q):
        return s*math.atanh(math.pow(q*math.pow(2,r),1/r)-1)+l

'''
Burr VIII distribution
'''
class burr8(Distribution):
    @staticmethod
    def random(r,l,s):
        return s*math.log(math.tan(math.pi/2*math.pow(rg0(),1/r)))+l
    @staticmethod
    def pdf(r,l,s,x):
        return 1/s*r*math.exp((x-l)/s)*math.pow(2/math.pi,r)*math.pow(math.atan(math.exp((x-l)/s)),r-1)/(1+math.exp(2*(x-l)/s))
    @staticmethod
    def cdf(r,l,s,x):
        return math.pow(2/math.pi*math.atan(math.exp((x-l)/s)),r)
    @staticmethod
    def median(r,l,s):
        return s*math.log(math.tan(math.pi/2*math.pow(1/2,1/r)))+l
    @staticmethod
    def ppf(r,l,s,q):
        return s*math.log(math.tan(math.pi/2*math.pow(q,1/r)))+l

'''
Burr X distribution
'''
class burr10(Distribution):
    @staticmethod
    def random(r,k,l,s):
        return math.log(math.pow((2/(1-rg0())-2)/k+1,1/r)-1)*s+l
    @staticmethod
    def pdf(r,k,l,s,x):
        xx=(x-l)/s
        return 1/s*2*math.exp(xx)*math.pow(1+math.exp(xx),r-1)*k*r/math.pow(2+(-1+math.pow(1+math.exp(xx),r))*k,2)
    @staticmethod
    def cdf(r,k,l,s,x):
        return 1-2/(2+k*math.pow(1+math.exp((x-l)/s),r)-1)
    @staticmethod
    def median(r,k,l,s):
        return math.log(math.pow(2/k+1,1/r)-1)*s+l
    @staticmethod
    def ppf(r,k,l,s,q):
        return math.log(math.pow((2/(1-q)-2)/k+1,1/r)-1)*s+l

'''
Burr XII distribution
'''
class burr12(Distribution):
    @staticmethod
    def pdf(c,k,x):
        if(c<=0 or k<=0 or x<=0):
            raise InvalidInputError("c, k, and x must be bigger than 0")
        return c*k*math.pow(x,c-1)/math.pow(1+math.pow(x,c),k+1)
    @staticmethod
    def cdf(k,c,x):
        if(c<=0 or k<=0 or x<=0):
            raise InvalidInputError("c, k, and x must be bigger than 0")
        return 1-math.pow(1+math.pow(x,c),-k)
    @staticmethod
    def random(c,k):
        if(c<=0 or k<=0):
            raise InvalidInputError("c and k must be bigger than 0")
        return math.pow(math.pow(1-u,-1/k)-1,1/c)
    @staticmethod
    def median(c,k):
        if(c<=0 or k<=0):
            raise InvalidInputError("c and k must be bigger than 0")
        math.pow(math.pow(2,1/k)-1,1/c)
    @staticmethod
    def ppf(c,k,q):
        if(c<=0 or k<=0):
            raise InvalidInputError("c and k must be bigger than 0")
        math.pow(math.pow(1-q,-1/k)-1,1/c)

'''
Cauchy distribution
'''
class cauchy(Distribution):
    @staticmethod
    def pdf(gmma,x0,x):
        if(gmma<=0):
            raise InvalidInputError("gamma must be bigger than 0")
        return 1/(math.pi*gmma*(1+((x-x0)/gmma)**2))
    @staticmethod
    def cdf(gmma,x0,x):
        if(gmma<=0):
            raise InvalidInputError("gamma must be bigger than 0")
        return 1/math.pi*math.atan((x-x0)/gmma)+1/2
    @staticmethod
    def random(gmma,x0):
        if(gmma<=0):
            raise InvalidInputError("gamma must be bigger than 0")
        return x0+gmma*math.tan(math.pi*(r.random()-1/2))
    @staticmethod
    def kurtosis(gmma,x0):
        return None
    @staticmethod
    def mean(gmma,x0):
        return None
    @staticmethod
    def median(gmma,x0):
        if(gmma<=0):
            raise InvalidInputError("gamma must be bigger than 0")
        return x0+gmma*math.tan(0)
    @staticmethod
    def mode(gmma,x0):
        if(gmma<=0):
            raise InvalidInputError("gamma must be bigger than 0")
    @staticmethod
    def variance(gmma,x0):
        return None
    @staticmethod
    def stddev(gmma,x0):
        return None
    @staticmethod
    def entropy(gmma,x0):
        if(gmma<=0):
            raise InvalidInputError("gamma must be bigger than 0")
        return math.log(gmma)+math.log(4*math.pi)
    @staticmethod
    def skewness(gmma,x0):
        return None
    @staticmethod
    def ppf(gmma,x0):
        if(gmma<=0):
            raise InvalidInputError("gamma must be bigger than 0")
        return x0+gmma*math.tan(math.pi*(q-1/2))

'''
Chi distribution
'''
class chi(Distribution):
    @staticmethod
    def random(k):
        return math.sqrt(chi2.random(k))
    @staticmethod
    def pdf(k,x):
        return math.pow(2,1-k/2)*math.pow(x,k-1)*math.exp(-x**2/2)/math.gamma(k/2)
    @staticmethod
    def mean(k):
        return math.sqrt(2)*math.gamma((k+1)/2)/math.gamma(k/2)
    @staticmethod
    def mode(k):
        if(k>=1):
            return math.sqrt(k-1)
    @staticmethod
    def variance(k):
        return k-(math.sqrt(2)*math.gamma((k+1)/2)/math.gamma(k/2))**2
    @staticmethod
    def stddev(k):
        return math.sqrt(k-(math.sqrt(2)*math.gamma((k+1)/2)/math.gamma(k/2))**2)
    @staticmethod
    def skewness(k):
        return (k-(math.sqrt(2)*math.gamma((k+1)/2)/math.gamma(k/2))**2)/chi.stddev(k)**3*(1-2*chi.variance(k))

'''
Chi-square distribution
'''
class chi2(Distribution): #Chi-squared
    @staticmethod
    def random(k):
        avg_=0
        for _ in range(k):
            avg_+=math.pow(normal.random(0,1),2)
        return avg_
    @staticmethod
    def pdf(k,x):
        return 1/(math.pow(2,k/2)*math.gamma(k/2))*math.pow(x,k/2-1)*math.exp(-x/2)
    @staticmethod
    def kurtosis(k):
        return 12/k
    @staticmethod
    def mean(k):
        return k
    @staticmethod
    def median(k):
        return k*math.pow(1-2/(9*k),3)
    @staticmethod
    def mode(k):
        return max(k-2,0)
    @staticmethod
    def variance(k):
        return 2*k
    @staticmethod
    def stddev(k):
        return math.sqrt(2*k)
    @staticmethod
    def skewness(k):
        return math.sqrt(8/k)

'''
Chi-square exponential distribution*
'''
class chi2exp(Distribution): #chi-squared exponential
    @staticmethod
    def random(k):
        return gammaexp.random(math.log(2),1,k/2)

'''
Complementary exponentiated inverse Weibull logarithmic distribution
'''
class compexpinvweibulllog(Distribution):
    @staticmethod
    def random(b,l,t):
        return math.pow(-math.log((1-(1-l)*rg0())/l)/t,-1/b)
    @staticmethod
    def pdf(b,l,t,x):
        return t*b*l/(-math.log(1-l)*(1-l*math.exp(-t*math.pow(x,-b))))*math.exp(y,-b-1)*math.exp(-t*math.pow(y,-b))
    @staticmethod
    def cdf(b,l,t,x):
        return (math.exp(l*math.exp(-t*math.pow(x,-b)))-1)/(math.exp(l)-1)
    @staticmethod
    def median(b,l,t):
        return math.pow(-math.log((1-(1-l)/2)/l)/t,-1/b)
    @staticmethod
    def ppf(b,l,t,q):
        return math.pow(-math.log((1-(1-l)*q)/l)/t,-1/b)

'''
Complementary exponentiated exponential-geometric lifetime distribution
'''
class complementaryexponentiatedexpgeolifetime(Distribution):
    @staticmethod
    def random(a,l,t):
        u=rg0()
        return -(math.log(1-math.pow(u/(t*(1-u)+u),1/a)))/l
    @staticmethod
    def pdf(a,l,t,x):
        return a*l*t*math.exp(-l*x)*math.pow(1-math.exp(-l*x),a-1)/math.pow(1-(1-t)*math.pow(1-math.exp(-l*x),a),2)
    @staticmethod
    def cdf(a,l,t,x):
        z=1-math.exp(-l*x)
        return 1-(1-math.pow(z,a))/(1-(1-t)*math.pow(z,a))
    @staticmethod
    def median(a,l,t):
        return -(math.log(1-math.pow(1/(t+1),1/a)))/l
    @staticmethod
    def ppf(a,l,t,q):
        return -(math.log(1-math.pow(q/(t*(1-q)+q),1/a)))/l

'''
Complementary exponentiated inverse Weibull poisson distribution
'''
class complementaryexponentiatedinvweibullpoisson(Distribution):
    @staticmethod
    def random(b,l,t):
        return math.pow(-math.log(math.log(rg0()*(math.exp(l)-1)+1)/l)/t,-1/b)
    @staticmethod
    def pdf(b,l,t,x):
        return t*b*l/(math.exp(l)-1)*math.pow(y,-b-1)*math.exp(-math.pow(x,-t*b))*math.exp(l*math.exp(-t*math.pow(x,-b)))
    @staticmethod
    def cdf(b,l,t,x):
        return (math.exp(l*math.exp(-t*math.pow(x,-b)))-1)/(math.exp(l)-1)
    @staticmethod
    def median(b,l,t):
        return math.pow(-math.log(math.log((math.exp(l)-1)/2+1)/l)/t,-1/b)
    @staticmethod
    def ppf(b,l,t,q):
        return math.pow(-math.log(math.log(q*(math.exp(l)-1)+1)/l)/t,-1/b)

'''
Compound Poisson distribution*
'''
class compoundpoisson(Distribution):
    @staticmethod
    def random(lmbda,mu):
        return poisson.random(mu*poisson.random(lmbda))

'''
Cosine distribution
'''
class cosine(Distribution):
    @staticmethod
    def random():
        return st.cosine.rvs()
    @staticmethod
    def pdf(x):
        return (1+math.cos(x))/(2*math.pi)
    @staticmethod
    def cdf(x):
        return (math.pi+x+math.sin(x))/(2*math.pi)

'''
Dagum distribution
'''
class dagum(Distribution):
    @staticmethod
    def pdf(p,a,b,x):
        if(p<=0 or a<=0 or b<=0 or x<=0):
            raise InvalidInputError("p, a, b, and x must be bigger than 0")
        return a*p/x*(math.pow(x/b,a*p)/math.pow(math.pow(x/b,a)+1,p+1))
    @staticmethod
    def cdf(p,a,b,x):
        if(p<=0 or a<=0 or b<=0 or x<=0):
            raise InvalidInputError("p, a, b, and x must be bigger than 0")
        return math.pow(1+math.pow(x/b,-a),-p)
    @staticmethod
    def random(p,a,b):
        if(p<=0 or a<=0 or b<=0):
            raise InvalidInputError("p, a, and b must be bigger than 0")
        return b*math.pow(math.pow(r.random(),-1/p)-1,-1/a)
    @staticmethod
    def median(p,a,b):
        if(p<=0 or a<=0 or b<=0):
            raise InvalidInputError("p, a, and b must be bigger than 0")
        return b*math.pow(-1+math.pow(2,1/p),-1/a)
    @staticmethod
    def ppf(p,a,b,q):
        if(p<=0 or a<=0 or b<=0):
            raise InvalidInputError("p, a, and b must be bigger than 0")
        return b*math.pow(math.pow(q,-1/p)-1,-1/a)

'''
Discrete generalized Rayleigh distribution
'''
class discretegenrayleigh(Distribution):
    @staticmethod
    def random(a,l):
        p=math.exp(-(l**2))
        return 1-math.sqrt(math.log(1-math.pow(rg0(),1/a))/math.log(p))
    @staticmethod
    def pdf(a,l,x):
        p=math.exp(-(l**2))
        return math.pow(1-math.pow(p,(x+1)**2),a)-math.pow(1-math.pow(p,x**2),a)
    @staticmethod
    def cdf(a,l,x):
        p=math.exp(-(l**2))
        return math.pow(1-math.pow(p,(x+1)**2),a)
    @staticmethod
    def median(a,l):
        p=math.exp(-(l**2))
        return 1-math.sqrt(math.log(1-math.pow(1/2,1/a))/math.log(p))
    @staticmethod
    def ppf(a,l,q):
        p=math.exp(-(l**2))
        return 1-math.sqrt(math.log(1-math.pow(q,1/a))/math.log(p))

'''
Discrete Weibull distribution
'''
class discreteweibull(Distribution):
    @staticmethod
    def random(aa,bb):
        n=r.random()
        while n==1:
            n=r.random()
        return math.pow(-math.log(1-n),1/bb)*aa-1
    @staticmethod
    def pdf(aa,bb,x):
        return math.exp(-math.pow(x/aa,bb))-math.exp(-math.pow((x+1)/aa,bb))
    @staticmethod
    def cdf(aa,bb,x):
        return 1-math.exp(-math.pow((x+1)/aa,bb))
    @staticmethod
    def median(aa,bb):
        return math.pow(-math.log(1/2),1/bb)*aa-1
    @staticmethod
    def ppf(aa,bb,q):
        return math.pow(-math.log(1-q),1/bb)*aa-1

'''
Double Weibull distribution
'''
class doubleweibull(Distribution):
    @staticmethod
    def random(a,b,c):
        y=r.random()
        if y<=0.5:
            return -b*math.pow(math.log(1/(2*y)),1/c)+a
        return b*math.pow(math.log(1/(2*y-1)),1/c)+a
    @staticmethod
    def pdf(a,b,c,x):
        z=(x-a)/b
        return c/2*math.pow(math.abs(z),c-1)*math.exp(-math.pow(math.abs(z),c))
    @staticmethod
    def cdf(a,b,c,x):
        z=(x-a)/b
        if(z<=0):
            return 1/2*math.exp(-math.pow(math.abs(z),c))
        return 1-1/2*math.exp(-math.pow(math.abs(z),c))
    @staticmethod
    def median(a,b,c):
        return a
    @staticmethod
    def mean(a,b,c):
        return a
    @staticmethod
    def variance(a,b,c):
        return math.gamma((c+2)/c*b*b)
    @staticmethod
    def stddev(a,b,c):
        return math.sqrt(math.gamma((c+2)/c*b*b))
    @staticmethod
    def skewness(a,b,c):
        return 0
    @staticmethod
    def ppf(a,b,c,q):
        if q<=0.5:
            return -b*math.pow(math.log(1/(2*q)),1/c)+a
        return b*math.pow(math.log(1/(2*q-1)),1/c)+a

'''
Erlang distribution
'''
class erlang(Distribution):
    @staticmethod
    def pdf(k,lmbda,x):
        if(k%1!=0 or k<=0):
            raise InvalidInputError("k must be a positive integer")
        if(lmbda<=0 or x<=0):
            raise InvalidInputError("lmbda and x must be bigger than 0")
        return math.pow(lmbda,k)*math.pow(x,k-1)*math.exp(-lmbda*x)/math.factorial(k-1)
    @staticmethod
    def random(k,lmbda):
        if(k%1!=0 or k<=0):
            raise InvalidInputError("k must be a positive integer")
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bigger than 0")
        total=1
        for _ in range(k):
            total+=exp.random(lmbda/(float(k)))
        return -math.log(total)/lmbda
    @staticmethod
    def kurtosis(k,lmbda):
        if(k%1!=0 or k<=0):
            raise InvalidInputError("k must be a positive integer")
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bigger than 0")
        return 6/k
    @staticmethod
    def mean(k,lmbda):
        if(k%1!=0 or k<=0):
            raise InvalidInputError("k must be a positive integer")
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bigger than 0")
        return k/lmbda
    @staticmethod
    def mode(k,lmbda):
        if(k%1!=0 or k<=0):
            raise InvalidInputError("k must be a positive integer")
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bigger than 0")
        return 1/lmbda*(k-1)
    @staticmethod
    def variance(k,lmbda):
        if(k%1!=0 or k<=0):
            raise InvalidInputError("k must be a positive integer")
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bigger than 0")
        return k/lmbda**2
    @staticmethod
    def stddev(k,lmbda):
        if(k%1!=0 or k<=0):
            raise InvalidInputError("k must be a positive integer")
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bigger than 0")
        return math.sqrt(k)/lmbda
    @staticmethod
    def skewness(k,lmbda):
        if(k%1!=0 or k<=0):
            raise InvalidInputError("k must be a positive integer")
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bigger than 0")
        return 2/math.sqrt(k)

'''
Exponentiated discrete Weibull distribution
'''
class expdiscweibull(Distribution): #exponentiated discrete weibull
    @staticmethod
    def random(p,a,g):
        return math.pow(math.pow(math.log(1-rg0()),(1/g))/math.log(p),1/a)-1
    @staticmethod
    def pdf(p,a,g,x):
        return math.pow(1-math.pow(p,math.pow(x+1,a)),g) - math.pow(1-math.pow(p,math.pow(x,a)))
    @staticmethod
    def cdf(p,a,g,x):
        return math.pow(1-math.pow(p,math.pow(x+1,a)),g)
    @staticmethod
    def median(p,a,g):
        return math.pow(math.pow(math.log(1/2),(1/g))/math.log(p),1/a)-1
    @staticmethod
    def ppf(p,a,g,q):
        return math.pow(math.pow(math.log(1-q),(1/g))/math.log(p),1/a)-1

'''
Exponentiated generalized extreme value distribution
'''
class expgenextremevalue(Distribution): #exponential generalized extreme value
    @staticmethod
    def random(b,k):
        if(k==0):
            return -math.log(-math.log(rg0())/b)
        return (1-math.pow(-math.log(rg0())/b,k))/k
    @staticmethod
    def pdf(b,k,x):
        if(k==0):
            return b*math.exp(-b*math.exp(-x))*math.exp(-x)
        return b*math.pow(1-k*x,1/k-1)*math.exp(-b*math.sqrt(1-k*x))
    @staticmethod
    def cdf(b,k,x):
        if(k==0):
            return math.exp(-b*math.exp(-x))
        return math.exp(-b*math.pow(1-k*x,1/k))
    @staticmethod
    def mean(b,k):
        if(k!=0):
            return 1/k*(1-math.pow(b,-k)*math.gamma(k+1))
    @staticmethod
    def median(b,k):
        if(k==0):
            return -math.log(-math.log(1/2)/b)
        return (1-math.pow(-math.log(1/2)/b,k))/k
    @staticmethod
    def ppf(b,k,q):
        if(k==0):
            return -math.log(-math.log(q)/b)
        return (1-math.pow(-math.log(q)/b,k))/k

'''
Exponentiated generalized Frechet distribution
'''
class expgenfrechet(Distribution): #exponentiated generalized frechet
    @staticmethod
    def random(a,b,l,s):
        return s*l/math.log(1-math.pow(1-math.pow(rg0(),1/b),1/a))
    @staticmethod
    def pdf(a,b,l,s,x):
        return a*b*l*math.pow(s,l)*math.pow(x,-l-1)*math.exp(-math.pow(s/x,l))*math.pow(1-math.exp(-math.pow(s/x,l)),a-1)*math.pow(1-math.pow(1-math.exp(-math.pow(s/x,l)),a),b-1)
    @staticmethod
    def cdf(a,b,l,s,x):
        return math.pow(1-math.pow(1-math.exp(-math.pow(s/x,l)),a),b)
    @staticmethod
    def median(a,b,l,s):
        return s*l/math.log(1-math.pow(1-math.pow(1/2,1/b),1/a))
    @staticmethod
    def ppf(a,b,l,s,q):
        return s*l/math.log(1-math.pow(1-math.pow(q,1/b),1/a))

'''
Exponentiated generalized Gumbel distribution
'''
class expgengumbel(Distribution): #exponentiated generalized gumbel
    @staticmethod
    def random(a,b,m,s):
        return m-s*math.log(-math.log(1-math.pow(1-math.pow(rg0(),1/b),1/a)))
    @staticmethod
    def pdf(a,b,m,s,x):
        return math.pow(1-math.pow(1-math.exp(-math.exp(-(x-m)/s)),a),b)
    @staticmethod
    def cdf(a,b,m,s,x):
        return a*b/s*math.exp(-((x-m)/s+math.exp(-(x-m)/s)))*math.pow(1-math.exp(-math.exp(-(x-m)/s)),a-1)*math.pow(1-math.pow(1-math.exp(-math.exp(-(x-m)/s)),a),b-1)
    @staticmethod
    def median(a,b,m,s):
        return m-s*math.log(-math.log(1-math.pow(1-math.pow(1/2,1/b),1/a)))
    @staticmethod
    def ppf(a,b,m,s,q):
        return m-s*math.log(-math.log(1-math.pow(1-math.pow(q,1/b),1/a)))
'''
Exponentiated generalized inverse exponential distribution
'''
class expgeninvexp(Distribution): #exponentiated generalized inverse exponential
    @staticmethod
    def random(a,b,l):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(rg0(),1/b),1/a)),-1)
    @staticmethod
    def median(a,b,l):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(1/2,1/b),1/a)),-1)
    @staticmethod
    def ppf(a,b,l,q):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(q,1/b),1/a)),-1)

'''
Exponentiated generalized inverse Weibull distribution
'''
class expgeninvweibull(Distribution): #exponentiated generalized inverse weibull
    @staticmethod
    def random(a,b,l,t):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(rg0(),1/b),1/a)),-1/t)
    @staticmethod
    def pdf(a,b,l,t,x):
        return a*b*t*math.pow(l,t)*math.pow(x,-t-1)*math.exp(-math.pow(l/x,t))*math.pow(1-math.exp(-math.pow(l/x,t)),a-1)*math.pow(1-math.pow(1-math.exp(math.pow(l/x,t)),a),b-1)
    @staticmethod
    def cdf(a,b,l,t,x):
        return math.pow(1-math.pow(1-math.exp(math.pow(l/x,t)),a),b)
    @staticmethod
    def median(a,b,l,t):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(1/2,1/b),1/a)),-1/t)
    @staticmethod
    def ppf(a,b,l,t,q):
        return l*math.pow(-math.log(1-math.pow(1-math.pow(q,1/b),1/a)),-1/t)

'''
Exponentiated generalized linear exponential distribution
'''
class expgenlinearexp(Distribution): #exponentiated generalized linear exponential
    @staticmethod
    def random(a,b,c,d):
        u=rg0()
        if(b!=0):
            return -a/b+1/b*math.sqrt(a**2+2*b*math.pow(-math.log(1-math.pow(u,1/d)),1/c))
        return 1/a*math.pow(-math.log(1-math.pow(u,1/d)),1/c)
    @staticmethod
    def pdf(a,b,c,d,x):
        return c*d*(a+b*x)*math.pow(a*x+b/2*x**2,c-1)*math.pow(1-math.exp(-math.pow(a*x+b/2*x**2,c)),d-1)*math.exp(-math.pow(a*x+b/2*x**2,c))
    @staticmethod
    def cdf(a,b,c,d,x):
        return math.pow(1-math.exp(-math.pow(a*x+b/2*x**2,c)),d)
    @staticmethod
    def median(a,b,c,d):
        if(b!=0):
            return -a/b+1/b*math.sqrt(a**2+2*b*math.pow(-math.log(1-math.pow(1/2,1/d)),1/c))
        return 1/a*math.pow(-math.log(1-math.pow(1/2,1/d)),1/c)
    @staticmethod
    def ppf(a,b,c,d,q):
        if(b!=0):
            return -a/b+1/b*math.sqrt(a**2+2*b*math.pow(-math.log(1-math.pow(q,1/d)),1/c))
        return 1/a*math.pow(-math.log(1-math.pow(q,1/d)),1/c)

'''
Exponentiated generalized normal distribution
'''
class expgennormal(Distribution): #exponentiated generalized normal
    @staticmethod
    def random(a,b,m,s):
        u=rg0()
        return st.norm.ppf(s*(1-math.pow(1-math.pow(u,1/b),1/a))+m,0,1)
    @staticmethod
    def median(a,b,m,s):
        return st.norm.ppf(s*(1-math.pow(1-math.pow(1/2,1/b),1/a))+m,0,1)
    @staticmethod
    def ppf(a,b,m,s,q):
        return st.norm.ppf(s*(1-math.pow(1-math.pow(q,1/b),1/a))+m,0,1)

'''
Exponential-geometric distribution
'''
class expgeo(Distribution):
    @staticmethod
    def random(l,t):
        u=rg0()
        return -math.log(t*(1-u)/(t*(1-u)+u))/l
    @staticmethod
    def pdf(l,t,x):
        return l*t*math.exp(-l*x)/math.pow(math.exp(-l*x)*(1-t)+t,2)
    @staticmethod
    def cdf(l,t,x):
        return 1-math.exp(-l*x)/(math.exp(-l*x)*(1-t)+t)
    @staticmethod
    def median(l,t):
        return -math.log(t/(t+1))/l
    @staticmethod
    def mean(l,t):
        return -math.log(t)/(l*(1-t))
    @staticmethod
    def mode(l,t):
        return 1/l*math.log((1-t)/t)
    @staticmethod
    def ppf(l,t,q):
        return -math.log(t*(1-q)/(t*(1-q)+q))/l

'''
Exponentiated Kumaraswamy distribution
'''
class expkumaraswamy(Distribution): #exponentiated kumaraswamy
    @staticmethod
    def random(a,b,g):
        return math.pow(1-math.pow(1-math.pow(rg0(),1/g),1/b),1/a)
    @staticmethod
    def median(a,b,g):
        return math.pow(1-math.pow(1-math.pow(1/2,1/g),1/b),1/a)
    @staticmethod
    def ppf(a,b,g,q):
        return math.pow(1-math.pow(1-math.pow(q,1/g),1/b),1/a)

'''
Exponentiated Kumaraswamy-Dagum distribution
'''
class expkumaraswamydagum(Distribution): #exponentiated kumaraswamy-dagum
    @staticmethod
    def random(a,d,l,p,t):
        return math.pow(l,1/d)*math.pow(math.pow(1-math.pow((1-math.pow(rg0(),1/t)),1/p),-1/a)-1,-1/d)
    @staticmethod
    def pdf(a,d,l,p,t,x):
        return math.pow(1-math.pow(1-math.pow(1+l*math.pow(x,-d),-a),p),t)
    @staticmethod
    def cdf(a,d,l,p,t,x):
        return a*l*d*p*t*math.pow(x,-d-1)*math.pow(1+l*math.pow(x,-d),-a-1)*math.pow(1-math.pow(1+l*math.pow(x,-d),-a),p-1)*math.pow(1-math.pow(1-math.pow(1+l*math.pow(x,-d),-a),p),t-1)
    @staticmethod
    def median(a,d,l,p,t):
        return math.pow(l,1/d)*math.pow(math.pow(1-math.pow((1-math.pow(1/2,1/t)),1/p),-1/a)-1,-1/d)
    @staticmethod
    def ppf(a,d,l,p,t,q):
        return math.pow(l,1/d)*math.pow(math.pow(1-math.pow((1-math.pow(q,1/t)),1/p),-1/a)-1,-1/d)

'''
Exponentiated Kumaraswamy-inverse Weibull distribution
'''
class expkumaraswamyinvweibull(Distribution):
    @staticmethod
    def random(a,b,e,l,t):
        return a/math.pow(-math.log(math.pow(1-math.pow(1-math.pow(rg0(),1/t),1/e),1/l)),1/b)
    @staticmethod
    def pdf(a,b,e,l,t,x):
        b*l*t*e*math.pow(a,b)*math.pow(x,-b-1)*math.exp(-l*math.pow(a/x,b))*math.pow(1-math.exp(-l*math.pow(a/x,b)),e-1)
    @staticmethod
    def cdf(a,b,e,l,t,x):
        return math.pow(1-math.pow(1-math.exp(-l*math.pow(a/x,b)),e),t)
    @staticmethod
    def median(a,b,e,l,t):
        return a/math.pow(-math.log(math.pow(1-math.pow(1-math.pow(1/2,1/t),1/e),1/l)),1/b)
    @staticmethod
    def ppf(a,b,e,l,t,q):
        return a/math.pow(-math.log(math.pow(1-math.pow(1-math.pow(q,1/t),1/e),1/l)),1/b)

'''
Exponential-logistic distribution
'''
class explog(Distribution): #exponentiated logistic
    @staticmethod
    def pdf(p,bb,x):
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        if(x<0):
            raise InvalidInputError("x must be non-negative")
        return 1/(-math.log(p))*(bb*(1-p)*math.exp(-bb*x))/(1-(1-p)*math.exp(-bb*x))
    @staticmethod
    def cdf(p,bb,x):
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        if(x<0):
            raise InvalidInputError("x must be non-negative")
        return 1-math.log(1-(1-p)*math.exp(-bb*x))/math.log(p)
    @staticmethod
    def random(p,bb):
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return 1/bb * math.log((1-p)/(1-math.pow(p,r.random())))
    @staticmethod
    def median(p,bb):
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return math.log(1+math.sqrt(p))/bb
    @staticmethod
    def mode(p,bb):
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return 0
    @staticmethod
    def ppf(p,bb,q):
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return 1/bb * math.log((1-p)/(1-math.pow(p,q)))

'''
Exponential-Lomax distribution
'''
class explomax(Distribution): #exponential lomax
    @staticmethod
    def random(a,b,g):
        return b*(math.pow(-math.log(1-rg0())/l,1/a)-1)
    @staticmethod
    def pdf(a,b,g,x):
        return a*l/b*math.pow(b/(x+b),-a+1)*math.exp(-l*math.pow(b/(x+b),-a))
    @staticmethod
    def cdf(a,b,g,x):
        return 1-math.exp(-l*math.pow(b/(x+b),-a))
    @staticmethod
    def median(a,b,g):
        return b*(math.pow(-math.log(1/2)/l,1/a)-1)
    @staticmethod
    def ppf(a,b,g,q):
        return b*(math.pow(-math.log(1-q)/l,1/a)-1)

'''
Exponentially-modified normal distribution
'''
class expmodnorm(Distribution): #exponentially modified normal
    @staticmethod
    def random(mu,sigma2,lmbda):
        return normal.random(mu,sigma2)+exp.random(lmbda)
    @staticmethod
    def pdf(mu,lmbda,sigma2,x):
        return lmbda/2*math.exp(lmbda/2*(2*mu+lmbda*sigma2-2*x))*math.erfc((mu+lmbda*sigma2-x)/(math.sqrt(2*sigma2)))
    @staticmethod
    def cdf(mu,lmbda,sigma2,x):
        u=lmbda*(x-mu)
        v=lmbda*math.sqrt(sigma2)
        return st.norm.cdf(u,0,v)-math.exp(-u+v**2/2+math.log(st.norm.cdf(u,v**2,v)))
    @staticmethod
    def kurtosis(mu,lmbda,sigma2):
        return 3*(1+2/(sigma2*lmbda**2)+3/(lmbda**4*sigma2**2))/math.pow(1+1/(lmbda**2*sigma2),2)-3
    @staticmethod
    def mean(mu,lmbda,sigma2):
        return mu+1/lmbda
    @staticmethod
    def variance(mu,lmbda,sigma2):
        return sigma2+1/lmbda**2
    @staticmethod
    def stddev(mu,lmbda,sigma2):
        return math.sqrt(sigma2+1/lmbda**2)
    @staticmethod
    def skewness(mu,lmbda,sigma2):
        return 2/math.pow(sigma2*lmbda**2,3/2)*math.pow(1+1/(sigma2*lmbda**2),-3/2)

'''
Exponential-normal distribution
'''
class expnorm(Distribution): #exponential normal
    @staticmethod
    def random(k):
        return st.exponnorm.rvs(K)

'''
Exponential distribution
'''
class exponential(Distribution):
    @staticmethod
    def pdf(lmbda,x):
        if(x<0 or lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0 and x must be non-negative")
        return lmbda*math.exp(-lmbda*x)
    @staticmethod
    def cdf(lmbda,x):
        if(x<0 or lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0 and x must be non-negative")
        return 1-math.exp(-lmbda*x)
    @staticmethod
    def random(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0")
        return -math.log(rg0())/lmbda
    @staticmethod
    def kurtosis(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0")
        return 6
    @staticmethod
    def mean(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0")
        return 1/lmbda
    @staticmethod
    def median(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0")
        return math.log(2)/lmbda
    @staticmethod
    def mode(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0")
        return 0
    @staticmethod
    def variance(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0")
        return 1/(lmbda**2)
    @staticmethod
    def stddev(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0")
        return 1/lmbda
    @staticmethod
    def entropy(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0")
        return 1-math.log(lmbda)
    @staticmethod
    def skewness(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0")
        return 2
    @staticmethod
    def ppf(lmbda,q):
        if(lmbda<=0):
            raise InvalidInputError("lmbda must be bgger than 0")
        return -math.log(q)/lmbda
    @staticmethod
    def mle(x):
        return {'lambda':len(x)/sum(x)}

'''
Exponentiated exponential-binomial distribution
'''
class exponentiatedexpbin(Distribution): #exponentiated exponential binomial
    @staticmethod
    def random(a,l,n,t):
        return -1/l*math.log(1-math.pow(1/t*(1-math.pow(1-(1-math.pow(1-t,n))*rg0(),1/n)),1/a))
    @staticmethod
    def pdf(a,l,n,t,x):
        return l*a*n*t*math.exp(-l*x)*math.pow(1-math.exp(-l*x),a-1)*math.pow(1-t*math.pow(1-math.exp(-l*x),a),n-1)/(1-math.pow(1-t,n))
    @staticmethod
    def cdf(a,l,n,t,x):
        return (1-math.pow(1-t*math.pow(1-math.exp(-l*x),a),n))/(1-math.pow(1-t,n))
    @staticmethod
    def median(a,l,n,t):
        return -1/l*math.log(1-math.pow(1/t*(1-math.pow(1-(1-math.pow(1-t,n))/2,1/n)),1/a))
    @staticmethod
    def ppf(a,l,n,t,q):
        return -1/l*math.log(1-math.pow(1/t*(1-math.pow(1-(1-math.pow(1-t,n))*q,1/n)),1/a))

'''
Exponentiated exponential-geometric distribution
'''
class exponentiatedexpgeo(Distribution): #exponentiated exponential-geometric
    @staticmethod
    def random(a,l,t):
        u=rg0()
        return -math.log(1-math.pow((u*t)/(1-u+u*t),1/a))/l
    @staticmethod
    def pdf(a,l,t,x):
        return a*l*t*math.exp(-l*x)*math.pow(1-math.exp(-l*x),a-1)/math.pow(1-(1-t)*math.pow(1-math.exp(-l*x),a),2)
    @staticmethod
    def cdf(a,l,t,x):
        return 1-(1-math.pow(1-math.exp(-l*x),a))/(1-(1-t)*math.pow(1-math.exp(-l*x),a))
    @staticmethod
    def median(a,l,t):
        return -math.log(1-math.pow(t/(1+t),1/a))/l
    @staticmethod
    def ppf(a,l,t,q):
        return -math.log(1-math.pow((q*t)/(1-q+q*t),1/a))/l

'''
Exponentiated Frechet distribution
'''
class exponentiatedfrechet(Distribution):
    @staticmethod
    def random(a,l,s):
        u=rg0()
        return s/math.pow(-math.log(1-math.pow(1-u,1/a)),1/l)
    @staticmethod
    def pdf(a,l,s,x):
        return a*l*math.pow(s/l)/math.pow(x,-l-1)*math.pow(1-math.exp(-math.pow(s/x,l)),a-1)*math.exp(-math.pow(s/x,l))
    @staticmethod
    def cdf(a,l,s,x):
        return 1-math.pow(1-math.exp(-math.pow(s/x,l)),a)
    @staticmethod
    def median(a,l,s):
        return s/math.pow(-math.log(1-math.pow(1/2,1/a)),1/l)
    @staticmethod
    def ppf(a,l,s,q):
        return s/math.pow(-math.log(1-math.pow(1-q,1/a)),1/l)

'''
Exponentiated Gompertz distribution
'''
class exponentiatedgompertz(Distribution):
    @staticmethod
    def random(a,l,t):
        return math.log(1-math.log(1-math.pow(rg0(),1/t))/l)/a
    @staticmethod
    def median(a,l,t):
        return math.log(1-math.log(1-math.pow(1/2,1/t))/l)/a
    @staticmethod
    def ppf(a,l,t,q):
        return math.log(1-math.log(1-math.pow(q,1/t))/l)/a

'''
Exponentiated Gumbel distribution
'''
class exponentiatedgumbel(Distribution):
    @staticmethod
    def random(a,m,s):
        return -s*math.log(-math.log(1-math.pow(1-rg0(),1/a)))+m
    @staticmethod
    def cdf(a,m,s,x):
        return 1-math.pow(1-math.exp(-math.exp(-(x-m)/s)),a)
    @staticmethod
    def median(a,m,s):
        return -s*math.log(-math.log(1-math.pow(1/2,1/a)))+m
    @staticmethod
    def ppf(a,m,s,q):
        return -s*math.log(-math.log(1-math.pow(1-q,1/a)))+m

'''
Exponentiated Lomax-Poisson distribution
'''
class exponentiatedlomaxpoisson(Distribution):
    @staticmethod
    def random(a,b,g,l):
        u=rg0()
        return (math.pow(-math.pow(-math.log(-u*(math.exp(l)-1)+math.exp(l))/l,1/a)+1,-1/g)-1)/b
    @staticmethod
    def pdf(a,b,g,l,x):
        return l*a*g*b*math.exp(l)*math.pow(1-math.pow(1+b*x,-g),a-1)*math.exp(-l*math.pow(1-math.pow(1+b*x,-y),a))/((math.exp(l)-1)*math.pow(1+b*x,g+1))
    @staticmethod
    def cdf(a,b,g,l,x):
        return math.exp(l)*(1-math.exp(-l*math.pow(1-math.pow(1+b*x,-g),a)))/(math.exp(l)-1)
    @staticmethod
    def median(a,b,g,l):
        return (math.pow(-math.pow(-math.log((math.exp(l)/2-1)+math.exp(l))/l,1/a)+1,-1/g)-1)/b
    @staticmethod
    def ppf(a,b,g,l,q):
        return (math.pow(-math.pow(-math.log(-q*(math.exp(l)-1)+math.exp(l))/l,1/a)+1,-1/g)-1)/b

'''
Exponentiated modified Weibull extension distribution
'''
class exponentiatedmodweibullext(Distribution): #exponentiated modified weibull extension
    @staticmethod
    def random(a,b,g,l):
        return a*math.pow(-math.log(1-1/(a*l)*math.log(1-math.pow(rg0(),1/g))),a/b)
    @staticmethod
    def pdf(a,b,g,l,x):
        return l*b*g*math.pow(x/a,b-1)*math.exp(math.pow(x/a,b)+l*a*(1-math.exp(math.pow(x/a,b))))*math.pow(1-math.exp(l*a*(1-math.exp(math.pow(x/a,b)))),g-1)
    @staticmethod
    def cdf(a,b,g,l,x):
        return math.pow(1-math.exp(l*a*(1-math.exp(math.pow(x/a,b)))),g)
    @staticmethod
    def median(a,b,g,l):
        return a*math.pow(-math.log(1-1/(a*l)*math.log(1-math.pow(1/2,1/g))),a/b)
    @staticmethod
    def ppf(a,b,g,l,q):
        return a*math.pow(-math.log(1-1/(a*l)*math.log(1-math.pow(q,1/g))),a/b)

'''
Exponentiated Weibull-exponential distribution
'''
class exponentiatedweibullexp(Distribution): #exponentiated weibull exponential
    @staticmethod
    def random(a,c,g):
        return -math.log(1-math.pow(1-math.exp(-g*math.pow(-math.log(1-rg0()),1/a)),1/c))
    @staticmethod
    def pdf(a,c,g,x):
        return c*a/g*math.exp(-x)*math.pow(1-math.exp(-x),c-1)/(1-math.pow(1-math.exp(-x),c))*math.pow(-math.log(1-math.pow(1-math.exp(-x),c))/g,a-1)*math.exp(-math.pow(-math.log(1-math.pow(1-math.exp(-x),c))/g),a)
    @staticmethod
    def cdf(a,c,g,x):
        return 1-math.exp(-math.pow(-math.log(1-math.pow(1-math.exp(-x),c))/g,a))
    @staticmethod
    def median(a,c,g):
        return -math.log(1-math.pow(1-math.exp(-g*math.pow(-math.log(1/2),1/a)),1/c))
    @staticmethod
    def ppf(a,c,g,q):
        return -math.log(1-math.pow(1-math.exp(-g*math.pow(-math.log(1-q),1/a)),1/c))

'''
Exponentiated Weibull-logistic distribution
'''
class exponentiatedweibulllog(Distribution): #exponentiated weibull logarithmic
    @staticmethod
    def random(a,b,g,t):
        u=rg0()
        return math.pow(-math.log(1-math.pow((1-math.pow(1-t,u))/t,1/a)),1/g)/b
    @staticmethod
    def pdf(a,b,g,t,x):
        return a*t*g*math.pow(b,g)*math.pow(x,g-1)*math.exp(-math.pow(b*x,g))*math.pow(1,math.exp(-math.pow(b*x,g)),a-1)/(math.log(1-t)*(t*math.pow(1-math.exp(-math.pow(b*x,g)),a)-1))
    @staticmethod
    def cdf(a,b,g,t,x):
        return math.log(1-t*math.pow(1-math.exp(-math.pow(b*x,g)),a))/math.log(1-t)
    @staticmethod
    def median(a,b,g,t):
        return math.pow(-math.log(1-math.pow((1-math.pow(1-t,1/2))/t,1/a)),1/g)/b
    @staticmethod
    def ppf(a,b,g,t,q):
        return math.pow(-math.log(1-math.pow((1-math.pow(1-t,q))/t,1/a)),1/g)/b

'''
Exponential transmuted generalized Rayleigh distribution
'''
class exptransmutedgenrayleigh(Distribution): #exponential transmuted generalized rayleigh
    @staticmethod
    def random(a,b,d,l):
        i=1+l-math.sqrt((1+l)**2-4*l*math.pow(rg0(),1/d))/(2*l)
        return math.sqrt(-math.log(1-math.pow(i,1/a)))/b
    @staticmethod
    def pdf(a,b,d,l,x):
        2*a*d*b**2*x*math.exp(-(b*x)**2)*math.pow(1-math.exp(-(b*x)**2),a*d-1)*(1+l-2*l*math.pow(1-math.exp(-(b*x)**2),a))*math.pow(1+l-l*math.pow(1-math.exp(-(b*x)**2),a),d-1)
    @staticmethod
    def cdf(a,b,d,l,x):
        return math.pow(1-math.exp(-(b*x)**2),a*d)*math.pow(1+l-l*math.pow(1-math.exp(-(b*x)**2),a),d)
    @staticmethod
    def median(a,b,d,l):
        i=1+l-math.sqrt((1+l)**2-4*l*math.pow(1/2,1/d))/(2*l)
        return math.sqrt(-math.log(1-math.pow(i,1/a)))/b
    @staticmethod
    def ppf(a,b,d,l,q):
        i=1+l-math.sqrt((1+l)**2-4*l*math.pow(q,1/d))/(2*l)
        return math.sqrt(-math.log(1-math.pow(i,1/a)))/b

'''
Exponential transmuted Weibull distribution
'''
class exptransmutedweibull(Distribution): #exponential transmuted weibull
    @staticmethod
    def random(a,b,l,n):
        i=math.sqrt(1/l*(1-math.pow(rg0(),1/n))+((1-l)/l)**2/4)
        return a*math.pow(-math.log(i-((1-l)/l)/2),1/b)
    @staticmethod
    def pdf(a,b,l,n,x):
        return n*b/a*math.pow(x/a,b-1)*math.exp(-math.pow(x/a,b))*(1-l+2*l*math.exp(-math.pow(x/a,b)))*math.pow(1+(l-1)*math.exp(-math.pow(x/a,b))-l(math.exp(-2(math.pow(x/a,b)))),n-1)
    @staticmethod
    def cdf(a,b,l,n,x):
        return math.pow(1+(l-1)*math.exp(-math.pow(x/a,b))-l(math.exp(-2(math.pow(x/a,b)))),n)
    @staticmethod
    def median(a,b,l,n):
        i=math.sqrt(1/l*(1-math.pow(1/2,1/n))+((1-l)/l)**2/4)
        return a*math.pow(-math.log(i-((1-l)/l)/2),1/b)
    @staticmethod
    def ppf(a,b,l,n,q):
        i=math.sqrt(1/l*(1-math.pow(q,1/n))+((1-l)/l)**2/4)
        return a*math.pow(-math.log(i-((1-l)/l)/2),1/b)

'''
Exponentiated Weibull distribution
'''
class expweibull(Distribution): #exponentiated weibull
    @staticmethod
    def cdf(aa,k,lmbda,x):
        if(aa<=0 or k<=0 or lmbda<=0 or x<=0):
            raise InvalidInputError("All inputs must be positive")
        return math.pow(1-math.exp(-math.pow(x/lmbda,k)),aa)
    @staticmethod
    def pdf(aa,k,lmbda,x):
        if(aa<=0 or k<=0 or lmbda<=0 or x<=0):
            raise InvalidInputError("All inputs must be positive")
        return aa*(k/lmbda)*math.pow(x/lmbda,k-1)*math.pow(1-math.exp(-math.pow(x/lmbda,k)),aa-1)*math.exp(-math.pow(x/lmbda,k))
    @staticmethod
    def random(aa,k,lmbda):
        if(aa<=0 or k<=0 or lmbda<=0):
            raise InvalidInputError("All inputs must be positive")
        return aa*math.pow(-math.log(1-math.pow(r.random(),1/k)),1/lmbda)
    @staticmethod
    def median(aa,k,lmbda):
        if(aa<=0 or k<=0 or lmbda<=0):
            raise InvalidInputError("All inputs must be positive")
        return aa*math.pow(-math.log(1-math.pow(1/2,1/k)),1/lmbda)
    @staticmethod
    def ppf(aa,k,lmbda,q):
        if(aa<=0 or k<=0 or lmbda<=0):
            raise InvalidInputError("All inputs must be positive")
        return aa*math.pow(-math.log(1-math.pow(q,1/k)),1/lmbda)

'''
Exponentiated Weibull-geometric distribution
'''
class expweibullgeo(Distribution): #exponentiated weibull-geometric
    @staticmethod
    def random(a,b,g,t):
        u=rg0()
        return math.pow(-math.log(1-math.pow(u/(1-t*(1-u)),1/a)),1/g)/b
    @staticmethod
    def pdf(a,b,g,t,x):
        (1-t)*a*g*math.pow(b,g)*math.pow(x,g-1)*math.exp(-math.pow(b*x,g))*math.pow(1-math.exp(-math.pow(b*x,g)),a-1)/(1-t*math.pow(1-math.exp(-math.pow(b*x,g)),a))**2
    @staticmethod
    def cdf(a,b,g,t,x):
        return (1-t)*math.pow(1-math.exp(-math.pow(b*x,g)),a)/(1-t*math.pow(1-math.exp(-math.pow(b*x,g)),a))
    @staticmethod
    def median(a,b,g,t):
        return math.pow(-math.log(1-math.pow(1/(2-t),1/a)),1/g)/b
    @staticmethod
    def ppf(a,b,g,t,q):
        return math.pow(-math.log(1-math.pow(q/(1-t*(1-q)),1/a)),1/g)/b

'''
Exponentiated Weibull-Poisson distribution
'''
class expweibullpoisson(Distribution): #exponentiated weibull-poisson
    @staticmethod
    def random(a,b,g,t):
        u=rg0()
        return 1/b*math.pow(-math.log(1-math.pow(1/t*math.log(u*(math.exp(t)-1)+1),1/a)),1/g)
    @staticmethod
    def pdf(a,b,g,t,x):
        return a*g*t*math.pow(b,g)*math.pow(x,g-1)/(math.exp(t)-1)*math.exp(-math.pow(b*x,g))*math.pow(1-math.exp(-math.pow(b*x,g)),a-1)*math.exp(t*math.pow(1-math.exp(-math.pow(b*x,g)),a))
    @staticmethod
    def cdf(a,b,g,t,x):
        return (math.exp(t*math.pow(1-math.exp(-math.pow(b*x,g)),a))-1)/(math.exp(t)-1)
    @staticmethod
    def median(a,b,g,t):
        return 1/b*math.pow(-math.log(1-math.pow(1/t*math.log((math.exp(t)-1)/2+1),1/a)),1/g)
    @staticmethod
    def ppf(a,b,g,t,q):
        return 1/b*math.pow(-math.log(1-math.pow(1/t*math.log(q*(math.exp(t)-1)+1),1/a)),1/g)

'''
Extended generalized exponential distribution
'''
class extendedgenexp(Distribution): #extended generalized exponential
    @staticmethod
    def random(a,b,l):
        if(b==0):
            return -(1/l)*math.ln(1-math.pow(rg0(),1/a))
        return (1/(b*l))*(1-math.pow(1-math.pow(rg0(),1/a),b))
    @staticmethod
    def pdf(a,b,l,x):
        if(b==0):
            return a*l*math.pow(1-math.exp(-l*x),a-1)*math.exp(-l*x)
        return a*l*math.pow(1-math.pow(1-b*l*x,1/b),a-1)*math.pow(1-b*l*x,1/b-1)
    @staticmethod
    def cdf(a,b,l,x):
        if(b==0):
            return math.pow(1-math.exp(-l*x),a)
        return math.pow(1-math.pow(1-b*l*x,1/b),a)
    @staticmethod
    def median(a,b,l):
        if(b==0):
            return -(1/l)*math.ln(1-math.pow(1/2,1/a))
        return (1/(b*l))*(1-math.pow(1-math.pow(1/2,1/a),b))
    @staticmethod
    def ppf(a,b,l,q):
        if(b==0):
            return -(1/l)*math.ln(1-math.pow(q,1/a))
        return (1/(b*l))*(1-math.pow(1-math.pow(q,1/a),b))

'''
Extreme value distribution
'''
class extremevalue(Distribution):
    @staticmethod
    def pdf(mu,sigma,x):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return 1/sigma*math.exp(-(x-mu)/sigma)*math.exp(-math.exp(-(x-mu)/sigma))
    @staticmethod
    def cdf(mu,sigma,x):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return math.exp(-(x-mu)/sigma)
    @staticmethod
    def random(mu,sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return mu-sigma*math.log(-math.log(rg0()))
    @staticmethod
    def kurtosis(mu,sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return 2.4
    @staticmethod
    def mean(mu,sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return mu+sigma*np.euler_gamma
    @staticmethod
    def median(mu,sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return mu-sigma*math.log(math.log(2))
    @staticmethod
    def mode(mu,sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return mu
    @staticmethod
    def variance(mu,sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return sigma**2*math.pi**2/6
    @staticmethod
    def stddev(mu,sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return sigma*math.pi/math.sqrt(6)
    @staticmethod
    def entropy(mu,sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return math.log(sigma)+np.euler_gamma+1
    @staticmethod
    def skewness(mu,sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return 12*math.sqrt(6)*np.special.zeta(3,1)/(math.pi**3)
    @staticmethod
    def ppf(mu,sigma,q):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return mu-sigma*math.log(-math.log(q))

'''
Extreme value maximum distribution*
'''
class extremevaluemax(Distribution):
    @staticmethod
    def random(a,b):
        return a-b*math.log(math.log(1/rg0()))

'''
Extreme value minimum distribution*
'''
class extremevaluemin(Distribution):
    @staticmethod
    def random(a,b):
        return -a+b*math.log(math.log(1/rg0()))

'''
F distribution
'''
class f(Distribution):
    @staticmethod
    def random(d1,d2):
        if(d1<=0 or d2<=0):
            raise InvalidInputError("All inputs must be positive")
        return (chi2.random(n1)/n1)/(chi2.random(n2)/n2)
    @staticmethod
    def kurtosis(d1,d2):
        if(d1<=0 or d2<=0):
            raise InvalidInputError("All inputs must be positive")
        if(d2>8):
            return 12*(d1*(5*d2-22)*(d1+d2-2)+(d**2-4)*(d**2-2)**2)/(d1*(d2-6)*(d2-8)*(d1+d2-2))
    @staticmethod
    def mean(d1,d2):
        if(d1<=0 or d2<=0):
            raise InvalidInputError("All inputs must be positive")
        if(d2>2):
            return d2/(d2-2)
    @staticmethod
    def mode(d1,d2):
        if(d1<=0 or d2<=0):
            raise InvalidInputError("All inputs must be positive")
        if(d1>2):
            return (d1-2)/d1*d2/(d2+2)
    @staticmethod
    def variance(d1,d2):
        if(d1<=0 or d2<=0):
            raise InvalidInputError("All inputs must be positive")
        if(d2>4):
            return (2*d2**2*(d1+d2-2))/(d1*(d2-2)**2*(d2-4))
    @staticmethod
    def stddev(d1,d2):
        if(d1<=0 or d2<=0):
            raise InvalidInputError("All inputs must be positive")
        if(d2>4):
            return math.sqrt((2*d2**2*(d1+d2-2))/(d1*(d2-2)**2*(d2-4)))
    @staticmethod
    def skewness(d1,d2):
        if(d1<=0 or d2<=0):
            raise InvalidInputError("All inputs must be positive")
        if(d2>6):
            return ((2*d1+d2-2)*math.sqrt(8*(d2-4)))/((d2-6)*math.sqrt(d1*(d1+d2-2)))

'''
Feller-Pareto distribution*
'''
class fellerpareto(Distribution):
    @staticmethod
    def random(mu,sigma,gmma,d1,d2):
        return mu+sigma*math.pow(gamma.random(d1,1)/gamma.random(d2,1),gmma)

'''
Fisher-z distribution
'''
class fisherz(Distribution):
    @staticmethod
    def random(n,m):
        return math.log(f.random(n,m))/2
    @staticmethod
    def pdf(n,m,z):
        return 2*math.pow(n,n/2)*math.pow(m,m/2)/sp.beta(n/2,m/2)*math.exp(n*z)/math.pow(n*math.exp(2*z)+m,(n+m)/2)
    @staticmethod
    def mean(n,m):
        return m/(m-2)
    @staticmethod
    def mode(n,m):
        return m/(m+2)*(n-2)/n

'''
Fisk distribution
'''
class fisk(Distribution):
    @staticmethod
    def pdf(a,b,c,x):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        if(x<=a):
            raise InvalidInputError("x must be bigger than or equal to a")
        return c/b*math.pow((x-a)/b,c-1)*math.pow(1+math.pow((x-a)/b,c),-2)
    @staticmethod
    def cdf(a,b,c,x):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        if(x<=a):
            raise InvalidInputError("x must be bigger than or equal to a")
        return 1/(1+math.pow((y-a)/b,-c))
    @staticmethod
    def random(a,b,c):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        u=r.random()
        return a+b*math.pow(u/(1-u),1/c)
    @staticmethod
    def kurtosis(a,b,c):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        if(c>4):
            p=math.pi/c
            return -3*(p*b*math.csc(p))**4+12*b*((p*b)**3)*math.csc(p)**2*math.csc(2*p)-12*(p**2)*(b**4)*math.csc(p)*math.csc(3*p)+4**p*(b**4)*math.csc(4*p)
        return None
    @staticmethod
    def mean(a,b,c):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        if(c>1):
            return a+math.pi*b/c*math.csc(math.pi/c)
        return None
    @staticmethod
    def median(a,b,c):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        return a+b
    @staticmethod
    def mode(a,b,c):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        return a+b*math.pow((c-1)/(c+1),1/c)
    @staticmethod
    def variance(a,b,c):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        if(c>2):
            return 2*math.pi*(b**2)/c*math.csc(2*math.pi/c)-(math.pi**2)*(b**2)/(c**2)*(math.csc(math.pi/c)**2)
        return None
    @staticmethod
    def stddev(a,b,c):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        if(c>2):
            return math.sqrt(2*math.pi*(b**2)/c*math.csc(2*math.pi/c)-(math.pi**2)*(b**2)/(c**2)*(math.csc(math.pi/c)**2))
        return None
    @staticmethod
    def skewness(a,b,c):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        if(c>3):
            return (2*math.pi/c*math.csc(math.pi/c))**3-6*b*(math.pi*b/c)**2*math.csc(math.pi/c)*math.csc(2*math.pi/c)+3*math.pi*(b**3)/c*math.csc(3*math.pi/c)
        return None
    @staticmethod
    def ppf(a,b,c,q):
        if(b<=0 or c<=0):
            raise InvalidInputError("b and c must be bigger than 0")
        return a+b*math.pow(q/(1-q),1/c)

'''
Friemer, Mudholkar, Kollia, and Lin generalized Tukey-lambda distribution
'''
class fmlkl(Distribution): #Friemer, Mudholkar, Kollia, and Lin generalized tukey-lambda
    @staticmethod
    def random(l1,l2,l3,l4):
        n=r.random()
        return l1+(1/l2)*((math.pow(n,l3)-1)/l3-((math.pow(1-n,l4)-1)/l4))
    @staticmethod
    def median(l1,l2,l3,l4):
        return l1+(1/l2)*((math.pow(1/2,l3)-1)/l3-((math.pow(1-n,l4)-1)/l4))
    @staticmethod
    def quantile(l1,l2,l3,l4,q):
        return l1+(1/l2)*((math.pow(q,l3)-1)/l3-((math.pow(1-q,l4)-1)/l4))

'''
Folded normal distribution
'''
class foldednormal(Distribution):
    @staticmethod
    def random(m,s):
        return abs(normal.random(m,s))
    @staticmethod
    def pdf(m,s,x):
        c=math.abs(m)/s
        return math.sqrt(2/math.pi)*math.cosh(c*x)*math.exp(-(x**2+c**2)/2)
    @staticmethod
    def cdf(m,s,x):
        c=math.abs(m)/s
        return st.norm.cdf(x-c)-st.norm.cdf(-x-c)
    @staticmethod
    def mean(m,s):
        c=math.abs(m)/s
        return math.sqrt(2/math.pi)*math.exp(-c**2/2)+c*math.erf(c/math.sqrt(2))

'''
Frechet distribution
'''
class frechet(Distribution):
    @staticmethod
    def pdf(aa,m,s,x):
        if(aa<=0 or s<=0):
            raise InvalidInputError("aa and s must be positive")
        if(x<=m):
            raise InvalidInputError("x must be bigger than m")
        return aa/s*math.pow((x-m)/s,-1-aa)*math.exp(-math.pow((x-m)/s,-aa))
    @staticmethod
    def cdf(aa,m,s,x):
        if(aa<=0 or s<=0):
            raise InvalidInputError("aa and s must be positive")
        if(x<=m):
            raise InvalidInputError("x must be bigger than m")
        return math.exp(-math.pow((x-m)/s,-aa))
    @staticmethod
    def random(aa,m,s):
        if(aa<=0 or s<=0):
            raise InvalidInputError("aa and s must be positive")
        n=rg0()
        return -math.pow(math.log(n),-1/aa)*(m*math.pow(-math.log(n),1/aa)+s)
    @staticmethod
    def median(aa,m,s):
        if(aa<=0 or s<=0):
            raise InvalidInputError("aa and s must be positive")
        return m+s/math.pow(math.log(2),1/aa)
    @staticmethod
    def mode(aa,m,s):
        if(aa<=0 or s<=0):
            raise InvalidInputError("aa and s must be positive")
        return m+s*math.pow(aa/(1+aa),1/aa)
    @staticmethod
    def entropy(aa,m,s):
        if(aa<=0 or s<=0):
            raise InvalidInputError("aa and s must be positive")
        return 1+np.euler_gamma/aa+np.euler_gamma+math.log(s/aa)
    @staticmethod
    def ppf(aa,m,s,q):
        if(aa<=0 or s<=0):
            raise InvalidInputError("aa and s must be positive")
        return -math.pow(math.log(q),-1/aa)*(m*math.pow(-math.log(q),1/aa)+s)

'''
Gamma distribution
'''
class gamma(Distribution):
    @staticmethod
    def random(k,t):
        u=r.random()
        x=-2*math.log(1-math.pow(u,1/k))
        v=r.random()
        while(v>(math.pow(x,k/2.0)*math.exp(-x/2.0))/(math.pow(2,k-1)*math.pow((1-math.exp(-x/2.0)),k-1))):
            u=r.random()
            x=-2*math.log(1-math.pow(u,1/k))
            v=r.random()
        return t*v
    @staticmethod
    def pdf(k,t,x):
        return 1/(math.gamma(k)*math.pow(t,k)*math.pow(x,k-1)*math.exp(-x/t))
    @staticmethod
    def kurtosis(k,t):
        return 6/k
    @staticmethod
    def mean(k,t):
        return k*t
    @staticmethod
    def mode(k,t):
        if(k>=1):
            return (k-1)*t
        return None
    @staticmethod
    def variance(k,t):
        return k*t**2
    @staticmethod
    def stddev(k,t):
        return math.sqrt(k)*t
    @staticmethod
    def skewness(k,t):
        return 2/math.sqrt(k)

'''
Gamma-exponential distribution*
'''
class gammaexp(Distribution): #gamma-exponential
    @staticmethod
    def random(nu,lmbda,a):
        return -math.log(amoroso.random(0,math.exp(nu),a,1/lmbda))

'''
Gamma-shifted Gompertz distribution
'''
class gammashiftedgompertz(Distribution):
    @staticmethod
    def random(b,aa,bb):
        return max(exp.random(b),b-gamma.random(aa,bb)*math.log(-math.log(r.random())))

'''
G and H distribution
'''
class gandh(Distribution):
    @staticmethod
    def random(g,h):
        n=rg0()
        return math.exp(g*st.norm.ppf(n)-1)*(math.exp((h*math.pow(st.norm.ppf(n),2))/2)/g)
    @staticmethod
    def median(g,h):
        return math.exp(g*st.norm.ppf(1/2)-1)*(math.exp((h*math.pow(st.norm.ppf(1/2),2))/2)/g)
    @staticmethod
    def ppf(g,h,q):
        return math.exp(g*st.norm.ppf(q)-1)*(math.exp((h*math.pow(st.norm.ppf(q),2))/2)/g)

'''
Gauss hypergeometric distribution
'''
class gausshyper(Distribution):
    @staticmethod
    def random(a,b,c,z):
        return st.gausshyper.rvs(a,b,c,z)

'''
Gauss-Kuzmin distribution*
'''
class gausskuzmin(Distribution):
    @staticmethod
    def random():
        a=r.random()
        i=1
        while True:
            if(a<1-math.log(1-1/((1+i)**2),2)):
                i+=1
            else:
                return i

'''
Generalized beta prime distribution*
'''
class genbetaprime(Distribution): #generalized beta-prime
    @staticmethod
    def random(a,s,aa,gmma,bb):
        return a+s*math.pow((gamma.random(aa,1)/gamma.random(gmma,1)),1/bb)

'''
Generalized exponential distribution*
'''
class genexp(Distribution): #generalized exponential
    @staticmethod
    def random(l,a):
        return -math.log(1-math.pow(rg0(),1/a))/l
    @staticmethod
    def median(l,a):
        return -math.log(1-math.pow(1/2,1/a))/l
    @staticmethod
    def ppf(l,a,q):
        return -math.log(1-math.pow(q,1/a))/l

'''
Generalized exponentiated Poisson distribution
'''
class genexppoisson(Distribution): #generalized exponentiated poisson
    @staticmethod
    def random(a,b,l):
        return -math.log(1+math.log(1-math.pow(rg0(),1/a)*(1-math.exp(-l)))/l)/b
    @staticmethod
    def pdf(a,b,l,x):
        return a*l*b/math.pow(1-math.exp(-l),a)*math.pow(1-math.exp(-l+l*math.exp(-b*x)),a-1)*math.exp(-l-b*x+l*math.exp(-b*x))
    @staticmethod
    def cdf(a,b,l,x):
        return math.pow((1-math.exp(-l+l*math.exp(-b*x)))/(1-math.exp(-l)),a)
    @staticmethod
    def median(a,b,l):
        return -math.log(1+math.log(1-math.pow(1/2,1/a)*(1-math.exp(-l)))/l)/b
    @staticmethod
    def ppf(a,b,l,q):
        return -math.log(1+math.log(1-math.pow(q,1/a)*(1-math.exp(-l)))/l)/b

'''
Generalized extreme value distribution
'''
class genextremevalue(Distribution):
    @staticmethod
    def pdf(mu,sigma,xi,x):
        if(xi==0):
            return extremevalue.pdf(mu,sigma,x)
        if(xi>0 and x<mu-sigma/xi):
            raise InvalidInputError("x must be bigger than or equal to mu-sigma/xi")
        if(xi<0 and x>mu-sigma/xi):
            raise InvalidInputError("x must be less than or equal to mu-sigma/xi")
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        t=math.pow(1+xi*(x-mu)/sigma,-1/xi)
        return 1/sigma*math.pow(t,xi+1)*math.exp(t)
    @staticmethod
    def cdf(mu,sigma,xi,x):
        if(xi==0):
            return extremevalue.cdf(mu,sigma,x)
        if(xi>0 and x<mu-sigma/xi):
            raise InvalidInputError("x must be bigger than or equal to mu-sigma/xi")
        if(xi<0 and x>mu-sigma/xi):
            raise InvalidInputError("x must be less than or equal to mu-sigma/xi")
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        t=math.pow(1+xi*(x-mu)/sigma,-1/xi)
        return math.exp(-t)
    @staticmethod
    def random(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return mu+sigma*(math.pow(math.log(1/rg0()),-xi)-1)/xi
    @staticmethod
    def kurtosis(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        if(xi>=1/4):
            return float("infinity")
        return (math.gamma(1-4*xi)-4*math.gamma(1-xi)*math.gamma(1-3*xi)+6*math.gamma(1-2*xi)*(math.gamma(1-xi)**2)-3*(math.gamma(1-xi)**4))/math.pow(math.gamma(1-2*xi)-math.gamma(1-xi)**2,2)-3
    @staticmethod
    def mean(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        if(xi>=1):
            return float("infinity")
        return mu+sigma*(math.gamma(1-xi)-1)/xi
    @staticmethod
    def median(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return mu+sigma*(math.pow(math.log(2),-xi)-1)/xi
    @staticmethod
    def mode(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return mu+sigma*(math.pow(1+xi,-xi)-1)/xi
    @staticmethod
    def variance(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        if(xi>=1/2):
            return float("infinity")
        return sigma**2*(math.gamma(1-2*xi)-math.gamma(1-xi)**2)/xi**2
    @staticmethod
    def stddev(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        if(xi>=1/2):
            return float("infinity")
        return sigma/xi*math.sqrt(math.gamma(1-2*xi)-math.gamma(1-xi)**2)
    @staticmethod
    def entropy(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return math.log(sigma)+np.euler_gamma*(xi+1)+1
    @staticmethod
    def skewness(mu,sigma,xi):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        if(xi>=1/3):
            return float("infinity")
        return (math.abs(xi)/xi)*(math.gamma(1-3*xi)-3*math.gamma(1-xi)*math.gamma(1-2*xi)+2*math.gamma(1-xi)**3)/math.pow(math.gamma(1-2*xi)-math.gamma(1-xi)**2,3/2)
    @staticmethod
    def ppf(mu,sigma,xi,q):
        if(xi==0):
            return extremevalue.random(mu,sigma)
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return mu+sigma*(math.pow(math.log(1/q),-xi)-1)/xi

'''
Generalized gamma distribution
'''
class gengamma(Distribution): #generalized gamma
    @staticmethod
    def random(a,c):
        return st.gengamma.rvs(a,c)

'''
Generalized Gompertz distribution
'''
class gengompertz(Distribution): #generalized gompertz
    @staticmethod
    def random(c,l,t):
        return math.log(1-c/l*math.log(1-math.pow(rg0(),t)))/c
    @staticmethod
    def pdf(c,l,t,x):
        return t*l*math.exp(c*x)*math.exp(-l/c*(math.exp(c*x)-1))*math.pow(1-math.exp(-l/c*(math.exp(c*x)-1)),t-1)
    @staticmethod
    def cdf(c,l,t,x):
        return math.pow(1-math.exp(-l/c*(math.exp(c*x)-1)),t)
    @staticmethod
    def median(c,l,t):
        return math.log(1-c/l*math.log(1-math.pow(1/2,t)))/c
    @staticmethod
    def ppf(c,l,t,q):
        return math.log(1-c/l*math.log(1-math.pow(q,t)))/c

'''
Generalized Gumbel distribution*
'''
class gengumbel(Distribution): #generalized gumbel
    @staticmethod
    def random(u,lmbda,n):
        return gammaexp.random(u-lmbda*math.log(n),lmbda,n)

'''
Generalized inverse generalized exponential distribution
'''
class geninvgenexp(Distribution): #generalized inverse generalized exponential
    @staticmethod
    def random(a,g,l):
        return l*g/(math.log(1-math.pow(1-rg0(),1/a)))
    @staticmethod
    def pdf(a,g,l,x):
        return a*l*g*math.pow(x,-2)*math.exp(-g*l/x)*math.pow(1-math.exp(-y*l/x),a-1)
    @staticmethod
    def cdf(a,g,l,x):
        return 1-math.pow(1-math.exp(-g*l/x),a)
    @staticmethod
    def median(a,g,l):
        return l*g/(math.log(1-math.pow(1/2,1/a)))
    @staticmethod
    def ppf(a,g,l,q):
        return l*g/(math.log(1-math.pow(1-q,1/a)))

'''
Generalized inverse Weibull distribution
'''
class geninvweibull(Distribution):
    @staticmethod
    def random(a,b,g):
        return a*math.pow(-math.log(rg0())/g,-1/b)
    @staticmethod
    def pdf(a,b,g,x):
        return g*b*math.pow(a,b)*math.pow(x,-b-1)*math.exp(-g*math.pow(a/x,b))
    @staticmethod
    def cdf(a,b,g,x):
        return math.exp(-g*math.pow(a/x,b))
    @staticmethod
    def median(a,b,g):
        return a*math.pow(-math.log(1/2)/g,-1/b)
    @staticmethod
    def ppf(a,b,g,q):
        return a*math.pow(-math.log(q)/g,-1/b)

'''
Generalized inverted exponential distribution
'''
class geninvexp(Distribution):
    @staticmethod
    def random(a,l):
        return -l/math.log(1-math.pow(1-rg0(),1/a))
    @staticmethod
    def pdf(a,l,x):
        return (a*l/x**2)*math.exp(-l/x)*math.pow(1-math.exp(-l/x),a-1)
    @staticmethod
    def cdf(a,l,x):
        return 1-math.pow(1-math.exp(-l/x),a)
    @staticmethod
    def median(a,l):
        return -l/math.log(1-math.pow(1/2,1/a))
    @staticmethod
    def ppf(a,l,q):
        return -l/math.log(1-math.pow(1-q,1/a))

'''
Generalized lambda distribution
'''
class genlambda(Distribution):
    @staticmethod
    def random(l1,l2,l3,l4):
        n=r.random()
        return l1+(math.pow(n,l3)-math.pow((1-n),l4))/l2
    @staticmethod
    def median(l1,l2,l3,l4):
        return l1+(math.pow(1/2,l3)-math.pow(1/2,l4))/l2
    @staticmethod
    def quantile(l1,l2,l3,l4,q):
        return l1+(math.pow(q,l3)-math.pow((1-q),l4))/l2

'''
Generalized Lindley distribution
'''
class genlindley(Distribution): #generalized lindley
    @staticmethod
    def random(a,t,b):
        u=rg0()
        v=gamma.random(a,t)
        vv=gamma.random(a+1,t)
        if(u<=(t/(b+t))):
            return v
        return vv

'''
Generalized linear failure rate-geometric distribution
'''
class genlinearfailurerategeo(Distribution): #generalized linear failure rate-geometric
    @staticmethod
    def random(a,b,aa,p):
        u=rg0()
        return -a/b+math.sqrt(a**2-2*b*math.log(1-math.pow((1-p)*u)/(1-p*u),1/aa))/b
    @staticmethod
    def median(a,b,aa,p):
        u=rg0()
        return -a/b+math.sqrt(a**2-2*b*math.log(1-math.pow((1-p)/2)/(1-p/2),1/aa))/b
    @staticmethod
    def ppf(a,b,aa,p,q):
        u=rg0()
        return -a/b+math.sqrt(a**2-2*b*math.log(1-math.pow((1-p)*q)/(1-p*q),1/aa))/b

'''
Generalized logistic distribution
'''
class genlogistic(Distribution): #generalized logistic
    @staticmethod
    def random(aa,mu,sigma):
        n=rg0()
        return mu-sigma*math.log((1-math.pow(n,aa))/math.pow(n,aa))
    @staticmethod
    def pdf(aa,mu,sigma,x):
        return aa/math.pow(math.exp((x-mu)/sigma)*(1+math.exp(-(x-mu)/sigma)),aa+1)
    @staticmethod
    def cdf(aa,mu,sigma,x):
        return 1/math.pow(1+math.exp(-(x-mu)/sigma),alpha)
    @staticmethod
    def median(aa,mu,sigma):
        return mu-sigma*math.log((1-math.pow(1/2,aa))/math.pow(1/2,aa))
    @staticmethod
    def ppf(aa,mu,sigma,q):
        return mu-sigma*math.log((1-math.pow(q,aa))/math.pow(q,aa))

'''
Generalized logistic V distribution*
'''
class genlogistic5(Distribution): #generalized logistic V
    @staticmethod
    def random(aa,mu,sigma):
        return mu+sigma*(1/aa)*(1-math.pow(1/rg0()-1,aa))
    @staticmethod
    def pdf(aa,mu,sigma,x):
        z=(x-mu)/sigma
        return math.pow(1-aa*z,(1/aa)-1)/(1+math.pow(1-aa*z,1/aa))**2
    @staticmethod
    def cdf(aa,x):
        z=(x-mu)/sigma
        return 1/(1+math.pow(1-aa*z,1/aa))
    @staticmethod
    def median(aa,mu,sigma):
        return mu+sigma*(1/aa)*(1-math.pow(1,aa))
    @staticmethod
    def ppf(aa,mu,sigma,q):
        return mu+sigma*(1/aa)*(1-math.pow(1/q-1,aa))

'''
Generalized odd log-logistic distribution
'''
class genoddloglogistic(Distribution): #generalized odd log-logistic
    @staticmethod
    def random(a,t):
        u=rg0()
        return math.pow(math.pow(u/(1-u),1/a)/(1+math.pow(u/(1-u),1/a)),1/t)
    @staticmethod
    def median(a,t):
        return math.pow(1/2,1/t)
    @staticmethod
    def ppf(a,t,q):
        return math.pow(math.pow(q/(1-q),1/a)/(1+math.pow(q/(1-q),1/a)),1/t)

'''
Generalized Pareto distribution
'''
class genpareto(Distribution):
    @staticmethod
    def pdf(mu,sigma,xi,x):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        if(x<mu and xi>=0):
            raise InvalidInputError("x must be bigger than or equal to mu")
        if(xi<0 and (x<mu or x>(mu-sigma/xi))):
            raise InvalidInputError("x must be between mu and mu-sigma/xi")
        return 1/sigma*math.pow(1+xi*(x-mu)/sigma,-1/xi+1)
    @staticmethod
    def cdf(mu,sigma,xi,x):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        if(x<mu and xi>=0):
            raise InvalidInputError("x must be bigger than or equal to mu")
        if(xi<0 and (x<mu or x>(mu-sigma/xi))):
            raise InvalidInputError("x must be between mu and mu-sigma/xi")
        return 1-math.pow(1+xi*(x-mu)/sigma,-1/xi)
    @staticmethod
    def random(mu,sigma,xi):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        if(xi==0):
            return mu-sigma*math.log(rg0())
        return mu+(sigma*(math.pow(rg0(),-xi)-1)/xi)
    @staticmethod
    def kurtosis(mu,sigma,xi):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        if(xi<1/4):
            return 3*(1-2*xi)*(2*(xi**2)+xi+3)/((1-3*xi)*(1-4*xi))-3
        return None
    @staticmethod
    def mean(mu,sigma,xi):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        if(xi<1):
            return mu+sigma/(1+xi)
        return None
    @staticmethod
    def median(mu,sigma,xi):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        return mu+sigma*(math.pow(2,xi)-1)/xi
    @staticmethod
    def variance(mu,sigma,xi):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        if(xi<1/2):
            return sigma**2/((1-xi)**2*(1-2*xi))
        return None
    @staticmethod
    def stddev(mu,sigma,xi):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        if(xi<1/2):
            return sigma/math.sqrt((1-xi)**2*(1-2*xi))
        return None
    @staticmethod
    def entropy(mu,sigma,xi):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        return math.log(sigma)+xi+1
    @staticmethod
    def skewness(mu,sigma,xi):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        if(xi<1/3):
            return 2*(1+xi)*math.sqrt(1-2*xi)/(1-3*xi)
        return None
    @staticmethod
    def ppf(mu,sigma,xi,q):
        if(sigma<=0):
            raise InvalidInputError("sigma must be greater than 0")
        if(xi==0):
            return mu-sigma*math.log(rg0())
        return mu+(sigma*(math.pow(q,-xi)-1)/xi)

'''
Generalized Pearson VII distribution*
'''
class genpearson7(Distribution): #generalized pearson VII
    @staticmethod
    def random(a,s,m,b):
        n=r.random()
        if(n<0.5):
            return -halfgenpearson7.random(a,s,m,b)
        return halfgenpearson7.random(a,s,m,b)

'''
Generalized Rayleigh distribution
'''
class genrayleigh(Distribution):
    @staticmethod
    def random(a,l):
        return math.sqrt(-math.log(1-math.pow(rg0(),1/a)))/l
    @staticmethod
    def pdf(a,l):
        return 2*a*l**2*x*math.exp(-(l*x)**2)*math.pow(1-math.exp(-(l*x)**2),a-1)
    @staticmethod
    def cdf(a,l):
        return math.pow(1-math.exp(-(l*x)**2),a)
    @staticmethod
    def median(a,l):
        return math.sqrt(-math.log(1-math.pow(1/2,1/a)))/l
    @staticmethod
    def ppf(a,l,q):
        return math.sqrt(-math.log(1-math.pow(q,1/a)))/l

'''
Generalized Topp-Leone distribution
'''
class gentoppleone(Distribution): #generalized topp-leone
    @staticmethod
    def random(aa,bb):
        if(aa==1):
            return math.pow(r.random(),1/bb)
        return (-aa+math.sqrt(math.pow(aa,2)-4*(1-aa)*(-math.pow(r.random(),1/bb))))/(2-2*aa)
    @staticmethod
    def pdf(aa,bb,x):
        return bb*math.pow(aa*x-(aa-1)*x**2,bb-1)*(aa-2*(aa-1)*x)
    @staticmethod
    def cdf(aa,bb,x):
        return math.pow(x,bb)*math.pow(aa-(aa-1)*x,bb)
    @staticmethod
    def median(aa,bb):
        if(aa==1):
            return math.pow(1/2,1/bb)
        return (-aa+math.sqrt(math.pow(aa,2)-4*(1-aa)*(-math.pow(1/2,1/bb))))/(2-2*aa)
    @staticmethod
    def ppf(aa,bb,q):
        if(aa==1):
            return math.pow(q,1/bb)
        return (-aa+math.sqrt(math.pow(aa,2)-4*(1-aa)*(-math.pow(q,1/bb))))/(2-2*aa)

'''
Generalized Tukey-lambda distribution*
'''
class gentukeylambda(Distribution): #generalized tukey-lambda
    @staticmethod
    def random(a,d):
        u=r.random()
        return (math.pow(u,a-d)-1)/(a-d)-(math.pow(1-u,a+d)-1)/(a+d)
    @staticmethod
    def median(a,d):
        return (math.pow(1/2,a-d)-1)/(a-d)-(math.pow(1/2,a+d)-1)/(a+d)
    @staticmethod
    def ppf(a,d,q):
        return (math.pow(q,a-d)-1)/(a-d)-(math.pow(1-q,a+d)-1)/(a+d)

'''
Generalized type I extreme value distribution*
'''
class gentypeIextremevalue(Distribution): #generalized type I extreme value
    @staticmethod
    def random(a,b,m,s):
        return -s*math.log(b/s*(math.pow(rg0(),-1/a)-1))+m
    @staticmethod
    def median(a,b,m,s):
        return -s*math.log(b/s*(math.pow(1/2,-1/a)-1))+m
    @staticmethod
    def ppf(a,b,m,s,q):
        return -s*math.log(b/s*(math.pow(q,-1/a)-1))+m

'''
Generalized Weibull-exponential distribution
'''
class genweibullexp(Distribution): #generalized weibull-exponential
    @staticmethod
    def random(a,c,t,g):
        return -math.log(1-math.pow(1-math.exp(math.pow(-math.log(1-t*rg0())/g,1/a)),1/c))
    @staticmethod
    def pdf(a,c,g,t,x):
        z=1-math.exp(-t*x)
        return c*a/g*t*math.exp(-t*x)*math.pow(z,c-1)/(1-math.pow(z,c))*math.pow(-math.log(1-math.pow(z,c))/g,a-1)*math.pow(math.exp-(-math.log(1-math.pow(1-math.exp(-t*x),c))/g),a)
    @staticmethod
    def cdf(a,c,g,t,x):
        return 1-math.pow(math.exp-(-math.log(1-math.pow(1-math.exp(-t*x),c))/g),a)
    @staticmethod
    def median(a,c,t,g):
        return -math.log(1-math.pow(1-math.exp(math.pow(-math.log(1-t/2)/g,1/a)),1/c))
    @staticmethod
    def ppf(a,c,t,g,q):
        return -math.log(1-math.pow(1-math.exp(math.pow(-math.log(1-t*q)/g,1/a)),1/c))

'''
Geometric distribution
'''
class geometric(Distribution):
    @staticmethod
    def pdf(p,k):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        if(k%1!=0 or k<1):
            raise InvalidInputError("k must be a positive integer")
        return math.pow(1-p,k-1)*p
    @staticmethod
    def cdf(p,k):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        if(k%1!=0 or k<1):
            raise InvalidInputError("k must be a positive integer")
        return 1-math.pow(1-p,k)
    @staticmethod
    def random(p):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        return math.floor(math.log(rg0())/math.log(1-p))
    @staticmethod
    def kurtosis(p):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        return 6+(p**2)/(1-p)
    @staticmethod
    def mean(p):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        return 1/p
    @staticmethod
    def median(p):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        return math.ceil(-1/math.log(1-p,2))
    @staticmethod
    def mode(p):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        return 1
    @staticmethod
    def variance(p):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        return (1-p)/(p**2)
    @staticmethod
    def stddev(p):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        return math.sqrt(1-p)/p
    @staticmethod
    def entropy(p):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        return (-(1-p)*math.log(1-p,2)-p*math.log(p,2))/p
    @staticmethod
    def skewness(p):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        return (2-p)/math.sqrt(1-p)
    @staticmethod
    def ppf(p,q):
        if(p<=0 or p>1):
            raise InvalidInputError("p must be greater than zero and less than or equal to 1")
        return math.floor(math.log(q)/math.log(1-p))

'''
Geometric extreme exponential distribution*
'''
class geometricextremeexp(Distribution): #geometric extreme-exponential
    @staticmethod
    def random(gmma):
        n=r.random()
        while n==1:
            n=r.random()
        return math.log((gmma)/(1-n)+1-gmma)

'''
Gompertz distribution
'''
class gompertz(Distribution):
    @staticmethod
    def pdf(eta,b,x):
        if(eta<=0 or b<=0):
            raise InvalidInputError("eta and b must be positive")
        if(x<0):
            raise InvalidInputError("x must be non-negative")
        return b*eta*math.exp(b*x)*math.exp(eta)*math.exp(-eta*math.exp(b*x))
    @staticmethod
    def cdf(eta,b,x):
        if(eta<=0 or b<=0):
            raise InvalidInputError("eta and b must be positive")
        if(x<0):
            raise InvalidInputError("x must be non-negative")
        return 1-math.exp(-eta*(math.exp(b*x)-1))
    @staticmethod
    def random(eta,b):
        if(eta<=0 or b<=0):
            raise InvalidInputError("eta and b must be positive")
        return math.log((-math.log(1-r.random())/eta)+1)/b
    @staticmethod
    def median(eta,b):
        if(eta<=0 or b<=0):
            raise InvalidInputError("eta and b must be positive")
        return (1/b)*math.log(math.log(1/2)/eta+1)
    @staticmethod
    def ppf(eta,b,q):
        if(eta<=0 or b<=0):
            raise InvalidInputError("eta and b must be positive")
        return math.log((-math.log(1-q)/eta)+1)/b

'''
Gompertz-Makeham distribution*
'''
class gompertzmakeham(Distribution):
    @staticmethod
    def random(lmbda,xi):
        q=rg0()
        return math.log(1-math.log(1-q)/xi)/lmbda
    @staticmethod
    def median(lmbda,xi):
        return math.log(1-math.log(1/2)/xi)/lmbda
    @staticmethod
    def random(lmbda,xi,q):
        return math.log(1-math.log(1-q)/xi)/lmbda

'''
Gumbel distribution
'''
class gumbel(Distribution):
    @staticmethod
    def pdf(mu,bb,x):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        z=(x-mu)/bb
        return 1/bb*math.exp(-(z+math.exp(-z)))
    @staticmethod
    def cdf(mu,bb,x):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return math.exp(-math.exp(-(x-mu)/bb))
    @staticmethod
    def random(mu,bb):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return mu-bb*math.log(-math.log(rg0()))
    @staticmethod
    def kurtosis(mu,bb):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return 12/5
    @staticmethod
    def mean(mu,bb):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return mu+bb*np.euler_gamma
    @staticmethod
    def median(mu,bb):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return mu-bb*math.log(math.log(2))
    @staticmethod
    def mode(mu,bb):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return mu
    @staticmethod
    def variance(mu,bb):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return math.pi**2/6*(bb**2)
    @staticmethod
    def stddev(mu,bb):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return math.pi*bb/math.sqrt(6)
    @staticmethod
    def entropy(mu,bb):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return math.log(bb)+np.euler_gamma+1
    @staticmethod
    def skewness(mu,bb):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return 12*math.sqrt(6)*np.special.zeta(3,1)/(math.pi**3)
    @staticmethod
    def ppf(mu,bb,q):
        if(bb<=0):
            raise InvalidInputError("bb must be positive")
        return mu-bb*math.log(-math.log(q))

'''
Gumbel 2 distribution*
'''
class gumbel2(Distribution):
    @staticmethod
    def random(a,b):
        return math.log(rg0())/(a*b)
    @staticmethod
    def median(a,b):
        return math.log(1/2)/(a*b)
    @staticmethod
    def ppf(a,b,q):
        return math.log(q)/(a*b)

'''
Half-Cauchy distribution
'''
class halfcauchy(Distribution):
    @staticmethod
    def random(x,gmma):
        return x+gmma*math.tan(math.pi*(r.random()/2.0))
    @staticmethod
    def median(x,gmma):
        return x+gmma*math.tan(math.pi*(1/4))
    @staticmethod
    def ppf(x,gmma,q):
        return x+gmma*math.tan(math.pi*(q/2.0))

'''
Half generalized Pearson V distribution*
'''
class halfgenpearson7(Distribution): #half generalized pearson VII
    @staticmethod
    def random(a,s,bb,m):
        return genbetaprime.random(a,s,1/bb,m-1/bb,bb)

'''
Half-Laha distribution*
'''
class halflaha(Distribution):
    @staticmethod
    def random(a,s):
        return halfgenpearson7.random(a,s,1,4)

'''
Half-logistic distribution
'''
class halflogistic(Distribution):
    @staticmethod
    def random(mu,s):
        n=logistic.random(mu,s)
        while(n<0):
            n=logistic.random(mu,s)
        return n

'''
Half-normal distribution
'''
class halfnormal(Distribution):
    @staticmethod
    def pdf(sigma2,x):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        if(x<0):
            raise InvalidInputError("x must be non-negative")
        return math.sqrt(2/(math.pi*sigma2))*math.exp(-x**2/(2*sigma2))
    @staticmethod
    def cdf(sigma2,x):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        if(x<0):
            raise InvalidInputError("x must be non-negative")
        return math.erf(math.sqrt(x**2/(2*sigma2)))
    @staticmethod
    def random(sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        n=normal.random(0,sigma2)
        while(n<0):
            n=normal.random(0,sigma2)
        return n
    @staticmethod
    def mean(sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return math.sqrt(2*sigma2/math.pi)
    @staticmethod
    def variance(sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return sigma2*(1-2/math.pi)
    @staticmethod
    def stddev(sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return math.sqrt(sigma2*(1-2/math.pi))
    @staticmethod
    def entropy(sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return 1/2*math.log(math.pi*sigma2/2)+1/2
    @staticmethod
    def skewness(sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")

'''
Harris extended exponential distribution
'''
class harrisextexp(Distribution): #harris extended exponential
    @staticmethod
    def random(k,l,t):
        return math.pow(l*k,-1)*math.log((1-t)+t*math.pow(1-rg0(),-k))
    @staticmethod
    def pdf(l,k,t,x):
        return l*math.pow(t,1/k)*math.exp(-l*x)/math.pow(1-(1-t)*math.exp(-l*k*x),1+1/k)
    @staticmethod
    def cdf(l,k,t,x):
        return 1-math.pow(t*math.exp(-l*k*x)/(1-(1-t)*math.exp(-l*k*x)),1/k)
    @staticmethod
    def random(k,l,t):
        return math.pow(l*k,-1)*math.log((1-t)+t*math.pow(1/2,-k))
    @staticmethod
    def ppf(l,k,t,q):
        return math.pow(l*k,-1)*math.log((1-t)+t*math.pow(1-q,-k))

'''
Hermite distribution
'''
class hermite(Distribution):
    @staticmethod
    def random(a1,a2):
        return poisson.random(a1)+2*poisson.random(a2)
    @staticmethod
    def kurtosis(a1,a2):
        return (a1+16*a2)/(a1+4*a2)**2
    @staticmethod
    def mean(a1,a2):
        return a1+2*a2
    @staticmethod
    def variance(a1,a2):
        return a1+4*a2
    @staticmethod
    def stddev(a1,a2):
        return math.sqrt(a1+4*a2)
    @staticmethod
    def skewness(a1,a2):
        return (a1+8*a2)/math.pow(a1+4*a2,3/2)

'''
Hyperbolic secant distribution
'''
class hyperbolicsecant(Distribution):
    @staticmethod
    def pdf(x):
        return 1/2*math.sech(math.pi*x/2)
    @staticmethod
    def cdf(x):
        return 2/math.pi*math.arctan(math.exp(math.pi*x/2))
    @staticmethod
    def random():
        return 2/math.pi*math.log(math.tan(math.pi*r.random()/2))
    @staticmethod
    def kurtosis():
        return 2
    @staticmethod
    def mean():
        return 0
    @staticmethod
    def median():
        return 0
    @staticmethod
    def mode():
        return 0
    @staticmethod
    def variance():
        return 1
    @staticmethod
    def stddev():
        return 1
    @staticmethod
    def entropy():
        return 4/math.pi*0.915965594177219015054603
    @staticmethod
    def skewness():
        return 0
    @staticmethod
    def ppf(q):
        return 1/2*math.sech(math.pi*q/2)

'''
Hyperexponential distribution
'''
class hyperexp(Distribution): #hyperexponential
    @staticmethod
    def random(lmbdas,probs):
        a=r.random()
        probs2=probs/sum(probs)
        counter=0
        while(a>sum(probs2[:counter+1])):
            counter+=1
        return exp.random(lmbdas[counter])

'''
Hypergeometric distribution
'''
class hypergeo(Distribution): #hypergeometric
    @staticmethod
    def random(n,N,K):
        num_success=0
        for i in range(n):
            d=r.random()
            if(d<(K-num_success)/(N-i)):
                num_success+=1
        return num_success
    @staticmethod
    def kurtosis(n,N,K):
        return 1/(n*K*(N-K)*(N-n)*(N-2)*(N-3))*((N-1)*N**2*(N*(N+1)-6*K*(N-K)-6*n*(N-n))+6*n*K*(N-K)*(N-n)*(5*N-6))
    @staticmethod
    def mean(n,N,K):
        return n*K/N
    @staticmethod
    def mode(n,N,K):
        return math.floor((n+1)*(K+1)/(N+2))
    @staticmethod
    def variance(n,N,K):
        return n*K/N*(N-K)/N*(N-n)/(N-1)
    @staticmethod
    def stddev(n,N,K):
        return math.sqrt(n*K/N*(N-K)/N*(N-n)/(N-1))
    @staticmethod
    def skewness(n,N,K):
        return (N-2*K)*math.sqrt(N-1)*(N-2*n)/(math.sqrt(n*K*(N-K)*(N-n))*(N-2))

'''
Hypoexponential distribution*
'''
class hypoexp(Distribution): #hypoexponential
    @staticmethod
    def random(lmbdas):
        sum_=0
        for i in range(lmbdas):
            sum_+=exp.random(lmbdas[i])
        return sum_

'''
Inverse chi-square distribution
'''
class invchi2(Distribution): #inverse chi-squared
    @staticmethod
    def random(nu):
        k=chi2.random(nu)
        while(k==0):
            k=chi2.random(nu)
        return 1/k
    @staticmethod
    def pdf(nu,x):
        return math.pow(2,-nu/2)/math.gamma(nu/2)*math.pow(x,-nu/2-1)*math.exp(-1/(2*x))
    @staticmethod
    def kurtosis(nu):
        return 12*(5*nu-22)/((nu-6)*(nu-8))
    @staticmethod
    def mean(nu):
        if(nu>2):
            return 1/(nu-2)
    @staticmethod
    def mode(nu):
        return 1/(nu+2)
    @staticmethod
    def variance(nu):
        return 2/((nu-2)**2*(nu-4))
    @staticmethod
    def stddev(nu):
        return math.sqrt(2/((nu-2)**2*(nu-4)))
    @staticmethod
    def skewness(nu):
        return 4/(nu-6)*math.sqrt(2*(nu-4))

'''
Inverse exponential distribution*
'''
class invexp(Distribution): #inverse exponential
    @staticmethod
    def random(theta):
        return invgamma.random(theta,1)

'''
Inverse gamma distribution
'''
class invgamma(Distribution): #inverse gamma
    @staticmethod
    def random(a,b):
        return 1/gamma.random(a,b)
    @staticmethod
    def pdf(a,b,x):
        return math.pow(b,a)/math.gamma(a)*math.pow(x,-a-1)*math.exp(-b/x)
    @staticmethod
    def kurtosis(a,b):
        if(a>4):
            return (30*a-66)/((a-3)*(a-4))
    @staticmethod
    def mean(a,b):
        if(a>1):
            return b/(a-1)
    @staticmethod
    def mode(a,b):
        return b/(a+1)
    @staticmethod
    def variance(a,b):
        if(a>2):
            return b**2/((a-1)**2*(a-2))
    @staticmethod
    def stddev(a,b):
        if(a>2):
            return math.sqrt(b**2/((a-1)**2*(a-2)))
    @staticmethod
    def skewness(a,b):
        if(a>3):
            return 4*math.sqrt(a-2)/(a-3)

'''
Inverse normal distribution*
'''
class invnormal(Distribution): #inverse normal
    @staticmethod
    def random(mu,lmbda):
        nu=normal.random(0,1)
        y=nu*nu
        x=mu+(mu*mu*y/(2*lmbda))-(mu/(2*lmbda))*math.sqrt(4*mu*lmbda*y+mu*mu*y*y)
        z=r.random()
        if(z<=(mu/(mu+x))):
            return x
        return mu*mu/x

'''
Inverse paralogistic distribution*
'''
class invparalogistic(Distribution): #inverse paralogistic
    @staticmethod
    def random(bb):
        return genbetaprime.random(0,1,bb,1,bb)

'''
Inverse Weibull distribution*
'''
class invweibull(Distribution): #inverse weibull
    @staticmethod
    def random(c):
        return st.invweibull.rvs(c)

'''
Irwin-Hall distribution
'''
class irwinhall(Distribution):
    @staticmethod
    def random(n):
        if(n<=0 or n%1!=0):
            raise InvalidInputError("n must be a positive integer")
        sum_=0
        for _ in range(n):
            sum_+=r.random()
        return sum_
    @staticmethod
    def kurtosis(n):
        if(n<=0 or n%1!=0):
            raise InvalidInputError("n must be a positive integer")
        return -6/(5*n)
    @staticmethod
    def mean(n):
        if(n<=0 or n%1!=0):
            raise InvalidInputError("n must be a positive integer")
        return n/2
    @staticmethod
    def median(n):
        if(n<=0 or n%1!=0):
            raise InvalidInputError("n must be a positive integer")
        return n/2
    @staticmethod
    def mode(n):
        if(n<=0 or n%1!=0):
            raise InvalidInputError("n must be a positive integer")
        if(n==1):
            return 1
        return n/2
    @staticmethod
    def variance(n):
        if(n<=0 or n%1!=0):
            raise InvalidInputError("n must be a positive integer")
        return n/12
    @staticmethod
    def stddev(n):
        if(n<=0 or n%1!=0):
            raise InvalidInputError("n must be a positive integer")
        return math.sqrt(n/12)
    @staticmethod
    def skewness(n):
        if(n<=0 or n%1!=0):
            raise InvalidInputError("n must be a positive integer")
        return 0

'''
Johnson SB distribution
'''
class johnsonsb(Distribution):
    @staticmethod
    def random(delta,gmma,xi,lmbda):
        v=(normal.random(0,1)-gmma)/delta
        return xi+lmbda*(1/(1+math.exp(-v)))
    @staticmethod
    def pdf(delta,gmma,xi,lmbda,x):
        return delta*lmbda*math.exp(-1/2*(gmma+delta*math.log((x-xi)/(-x+xi+lmbda)))**2)/(math.sqrt(2*math.pi)*(x-xi)*(-x+xi+lmbda))
    @staticmethod
    def cdf(delta,gmma,xi,lmbda,x):
        if(xi<x and x<xi+lmbda/2):
            return 1/2*math.erfc(-(gmma+delta*math.log((x-xi)/(-x+xi+lmbda)))/math.sqrt(2))
        if(xi+lmbda/2<=x and x<xi+lmbda):
            return 1/2*(1+math.erf((gmma+delta*math.log((x-xi)/(-x+xi+lmbda)))/math.sqrt(2)))
    @staticmethod
    def median(delta,gmma,xi,lmbda):
        return xi+lmbda/(1+math.exp(gmma/delta))

'''
Johnson SL distribution
'''
class johnsonsl(Distribution):
    @staticmethod
    def random(delta,gmma,xi,lmbda):
        v=(normal.random(0,1)-gmma)/delta
        return xi+lmbda*math.exp(v)
    @staticmethod
    def pdf(delta,gmma,xi,lmbda,x):
        return delta*math.exp(-1/2*(gmma+delta*math.log((x-xi)/lmbda)))/(math.sqrt(2*math.pi)*(x-xi))
    @staticmethod
    def cdf(delta,gmma,xi,lmbda,x):
        if(xi<x and x<=xi+lmbda):
            return 1/2*math.erfc(-(gmma+delta*math.log((x-xi)/lmbda))/math.sqrt(2))
        if(x>xi+lmbda):
            return 1/2*(1+math.erf((gmma+delta*math.log((x-xi)/lmbda))/math.sqrt(2)))
    @staticmethod
    def kurtosis(delta,gmma,xi,lmbda):
        return -3+math.exp(2/delta**2)*(3+math.exp(1/delta**2)*(2+math.exp(1/delta**2)))
    @staticmethod
    def mean(delta,gmma,xi,lmbda):
        return xi+lmbda*math.exp((1-2*gmma*delta)/(2*delta**2))
    @staticmethod
    def median():
        return xi+lmbda*math.exp(-gmma/delta)
    @staticmethod
    def variance(delta,gmma,xi,lmbda):
        return math.exp((1-2*gmma*delta)/(delta**2))*lmbda**2*(-1+math.exp(1/delta**2))
    @staticmethod
    def stddev():
        return math.sqrt(math.exp((1-2*gmma*delta)/(delta**2))*lmbda**2*(-1+math.exp(1/delta**2)))
    @staticmethod
    def skewness(delta,gmma,xi,lmbda):
        return (2+math.exp(1/delta**2))*math.sqrt(-1+math.exp(1/delta**2))

'''
Johnson SU distribution
'''
class johnsonsu(Distribution):
    @staticmethod
    def random(delta,gmma,xi,lmbda):
        return lmbda*math.sinh((normal.random(0,1)-gmma)/delta)+xi
    @staticmethod
    def pdf(delta,gmma,xi,lmbda,x):
        return delta/(gmma*math.sqrt(2*math.pi))*1/math.sqrt(1+((x-xi)/lmbda)**2)*math.exp(-1/2*math.pow(gmma+delta*math.asinh((x-xi)/lmbda),2))
    @staticmethod
    def cdf(delta,gmma,xi,lmbda,x):
        return st.norm.cdf(gmma+delta*math.asinh((x-xi)/lmbda))
    @staticmethod
    def kurtosis(delta,gmma,xi,lmbda):
        z=1/delta**2
        return (4*math.exp((2+2*gmma*delta)*z)*(2+math.exp(z))+4*math.exp((2+6*gmma*delta)*z)*(2+math.exp(z))+6*math.exp(4*gmma*z)*(1+2*math.exp(z))+math.exp(2*z)*(-3+math.exp(2*z)*(3+math.exp(z)*(2+math.exp(z))))+math.exp((2+8*gmma*delta)*z)*(-3+math.exp(2*z)*(3+math.exp(z)*(2+math.exp(z)))))/math.pow(math.exp(z)+2*math.exp(2*z)+math.exp((1+4*gmma*delta)*z),2)
    @staticmethod
    def mean(delta,gmma,xi,lmbda):
        return xi-lmbda*math.exp(1/(2*delta**2))*math.sinh(gmma/delta)
    @staticmethod
    def variance(delta,gmma,xi,lmbda):
        return lmbda**2/2*(math.exp(1/delta**2)-1)*(math.exp(1/delta**2)*math.cosh(2*gmma/delta)+1)
    @staticmethod
    def stddev(delta,gmma,xi,lmbda):
        return math.sqrt(lmbda**2/2*(math.exp(1/delta**2)-1)*(math.exp(1/delta**2)*math.cosh(2*gmma/delta)+1))
    @staticmethod
    def skewness(delta,gmma,xi,lmbda):
        i=math.exp(1/(2*delta**2))*math.sqrt(-1+math.exp(1/delta**2))*(3*math.exp(2*gmma/delta)-3*math.exp(4*gmma/delta)+(2+math.exp(1/delta**2)*math.exp(1/delta**2))-(2+math.exp(1/delta**2))*math.exp((1+6*gmma*delta)/delta**2))/math.pow(math.exp(1/delta**2)+2*math.exp(2*gmma/delta)+math.exp((1+4*gmma*delta)/delta**2),3/2)

'''
K distribution*
'''
class k(Distribution):
    @staticmethod
    def random(a,b):
        return math.sqrt(exp.random(1)*gamma.random(a,b/a))

'''
K-prime distribution
'''
class kprime(Distribution):
    @staticmethod
    def random(a,b,nu1,nu2):
        return (b*normal.random(0,1)+a*math.sqrt(chi2.random(nu1)/nu1))/math.sqrt(chi2.random(nu2)/nu2)

'''
Kappa distribution*
'''
class kappa(Distribution):
    @staticmethod
    def random(h,k,xi,aa):
        n=rg0()
        return xi+(aa/k)*(1-math.pow(1-math.pow(n,h),k))
    @staticmethod
    def median(h,k,xi,aa):
        return xi+(aa/k)*(1-math.pow(1-math.pow(1/2,h),k))
    @staticmethod
    def random(h,k,xi,aa,q):
        return xi+(aa/k)*(1-math.pow(1-math.pow(q,h),k))

'''
Kolmogorov-Smirnov I distribution
'''
class kolmogorovsmirnovone(Distribution):
    @staticmethod
    def random(n):
        return st.ksone.rvs(n)

'''
Kolmogorov-Smirnov II large distribution
'''
class kolmogovsmirnovtwolarge(Distribution):
    @staticmethod
    def random(n):
        return st.kstwobign.rvs(n)

'''
Kumaraswamy distribution
'''
class kumaraswamy(Distribution):
    @staticmethod
    def pdf(a,b,x):
        if(a<=0 or b<=0):
            raise InvalidInputError("a and b must be positive")
        if(x<0 or x>1):
            raise InvalidInputError("x must be between 0 and 1 inclusive")
        return a*b*math.pow(x,a-1)*math.pow(1-math.pow(x,a),b-1)
    @staticmethod
    def cdf(a,b,x):
        if(a<=0 or b<=0):
            raise InvalidInputError("a and b must be positive")
        if(x<0 or x>1):
            raise InvalidInputError("x must be between 0 and 1 inclusive")
        return 1-math.pow(1-math.pow(x,a),b)
    @staticmethod
    def random(a,b):
        if(a<=0 or b<=0):
            raise InvalidInputError("a and b must be positive")
        return math.pow(1-math.pow(1-r.random(),1/b),1/a)
    @staticmethod
    def mean(a,b):
        if(a<=0 or b<=0):
            raise InvalidInputError("a and b must be positive")
        return b*math.gamma(1+1/a)*math.gamma(b)/math.gamma(1+1/a+b)
    @staticmethod
    def median(a,b):
        if(a<=0 or b<=0):
            raise InvalidInputError("a and b must be positive")
        return math.pow(1-math.pow(2,-1/b),1/a)
    @staticmethod
    def mode(a,b):
        if(a<=0 or b<=0):
            raise InvalidInputError("a and b must be positive")
        if(a>=1 and b>=1 and (a!=1 or b!=1)):
            return math.pow((a-1)/(a*b-1),1/a)
        return None
    @staticmethod
    def variance(a,b):
        if(a<=0 or b<=0):
            raise InvalidInputError("a and b must be positive")
        return (b*math.gamma(1+2/a)*math.gamma(b)/math.gamma(1+b+2/a))-(b*math.gamma(1+1/a)*math.gamma(b)/math.gamma(1+b+1/a))**2
    @staticmethod
    def stddev(a,b):
        if(a<=0 or b<=0):
            raise InvalidInputError("a and b must be positive")
    @staticmethod
    def ppf(a,b,q):
        if(a<=0 or b<=0):
            raise InvalidInputError("a and b must be positive")
        return math.pow(1-math.pow(1-q,1/b),1/a)

'''
Kumaraswamy 4 distribution*
'''
class kumaraswamy4(Distribution):
    @staticmethod
    def random(a,b,c,d):
        '''alpha, beta, min, max'''
        n=rg0()
        while n==1:
            n=rg0()
        return c+(d-c)*math.pow(1-math.pow(n,1/b),1/a)
    @staticmethod
    def median(a,b,c,d):
        '''alpha, beta, min, max'''
        return c+(d-c)*math.pow(1-math.pow(1/2,1/b),1/a)
    @staticmethod
    def ppf(a,b,c,d,q):
        '''alpha, beta, min, max'''
        return c+(d-c)*math.pow(1-math.pow(q,1/b),1/a)

'''
Kumaraswamy-Dagum distribution
'''
class kumaraswamydagum(Distribution):
    @staticmethod
    def random(a,b,bb,d,l):
        u=rg0()
        return math.pow((math.pow(1-math.pow(1-u,1/b),1/(-bb*a))-1)/l,-1/d)
    @staticmethod
    def pdf(a,b,bb,d,l,x):
        return a*b*bb*l*d*math.pow(x,-d-1)*math.pow(1+l*math.pow(x,-d),-bb*a-1)*math.pow(1-math.pow(1+l*math.pow(x,-d),-bb*a),-b-1)
    @staticmethod
    def cdf(a,b,bb,d,l,x):
        return 1-math.pow(1-math.pow(1+l*math.pow(x,-d),-bb*a),b)
    @staticmethod
    def median(a,b,bb,d,l):
        return math.pow((math.pow(1-math.pow(1/2,1/b),1/(-bb*a))-1)/l,-1/d)
    @staticmethod
    def ppf(a,b,bb,d,l,q):
        return math.pow((math.pow(1-math.pow(1-q,1/b),1/(-bb*a))-1)/l,-1/d)

'''
Kumaraswamy-exponentiated Pareto distribution
'''
class kumaraswamyexponentiatedpareto(Distribution):
    @staticmethod
    def random(a,b,l,m,t):
        u=rg0()
        return l/math.pow(1-math.pow(1-math.pow(1-u,1/b),1/(t*a)),1/m)
    @staticmethod
    def pdf(a,b,l,m,t,x):
        return a*b*t*m*math.pow(l,m)/math.pow(x,m+1)*math.pow(1-math.pow(l/x,m),t*a-1)*math.pow(1-math.pow(1-math.pow(l/x,m),t*a),b-1)
    @staticmethod
    def cdf(a,b,l,m,t,x):
        return 1-math.pow(1-math.pow(1-math.pow(l/x,m),t*a),b)
    @staticmethod
    def median(a,b,l,m,t):
        return l/math.pow(1-math.pow(1-math.pow(1/2,1/b),1/(t*a)),1/m)
    @staticmethod
    def ppf(a,b,l,m,t,q):
        return l/math.pow(1-math.pow(1-math.pow(1-q,1/b),1/(t*a)),1/m)

'''
Kumaraswamy flexible Weibull extension distribution
'''
class kumaraswamyflexibleweibullextension(Distribution):
    @staticmethod
    def random(a,b,aa,bb):
        u=rg0()
        i=math.log(-math.log(1-math.pow(1-math.pow(1-u,1/b),1/a)))
        return 1/(2*aa)*(i+math.sqrt(i**2+4*aa*bb))
    @staticmethod
    def pdf(a,b,aa,bb,x):
        return a*b*(aa+bb/x**2)*math.exp(aa*x-bb/x)*math.exp(-math.exp(aa*x-bb/x))*math.pow(1-math.exp(-math.exp(aa*x-bb/x)),a-1)*math.pow(1-math.pow(1-math.exp(-math.exp(aa*x-bb/x)),a),b-1)
    @staticmethod
    def cdf(a,b,aa,bb,x):
        return 1-math.pow(1-math.pow(1-math.exp(-math.exp(aa*x-bb/x)),a),b)
    @staticmethod
    def median():
        i=math.log(-math.log(1-math.pow(1-math.pow(1/2,1/b),1/a)))
        return 1/(2*aa)*(i+math.sqrt(i**2+4*aa*bb))
    @staticmethod
    def ppf(a,b,aa,bb,q):
        i=math.log(-math.log(1-math.pow(1-math.pow(1-q,1/b),1/a)))
        return 1/(2*aa)*(i+math.sqrt(i**2+4*aa*bb))

'''
Kumaraswamy generalized exponentiated exponential distribution
'''
class kumaraswamygenexpexp(Distribution): #kumaraswamy generalized exponentiated exponential
    @staticmethod
    def random(a,b,c,l):
        return 1/(l*math.log(1-math.pow(1-math.pow(1-rg0(),1/b),1/(a*c))))
    def median(a,b,c,l):
        return 1/(l*math.log(1-math.pow(1-math.pow(1/2,1/b),1/(a*c))))
    def ppf(a,b,c,l,q):
        return 1/(l*math.log(1-math.pow(1-math.pow(1-q,1/b),1/(a*c))))

'''
Kumaraswamy generalized exponentiated Pareto distribution
'''
class kumaraswamygenexppareto(Distribution): #kumaraswamy generalized exponentiated pareto
    @staticmethod
    def random(a,b,l,t):
        return math.pow(1-math.pow(1-math.pow(1-rg0(),1/b),1/(a*t)),-1/l)-1
    @staticmethod
    def pdf(a,b,l,t,x):
        return a*b*t*l*math.pow(1-math.pow(1+x,-l),t*a-1)*math.pow(1+x,-l-1)*math.pow(1-math.pow(1-math.pow(1+x,-l),t*a),b-1)
    @staticmethod
    def cdf(a,b,l,t,x):
        return 1-math.pow(1-math.pow(1-math.pow(1+x,-l),t*a),b)
    @staticmethod
    def median(a,b,l,t):
        return math.pow(1-math.pow(1-math.pow(1/2,1/b),1/(a*t)),-1/l)-1
    @staticmethod
    def ppf(a,b,l,t,q):
        return math.pow(1-math.pow(1-math.pow(1-q,1/b),1/(a*t)),-1/l)-1

'''
Kumaraswamy generalized exponentiated Weibull distribution
'''
class kumaraswamygenexpweibull(Distribution): #kumaraswamy generalized exponentiated weibull
    @staticmethod
    def random(a,b,aa,bb):
        return -math.log(1-math.pow(1-math.pow(1-rg0(),1/b),1/a))/(aa+bb)
    @staticmethod
    def median(a,b,aa,bb):
        return -math.log(1-math.pow(1-math.pow(1/2,1/b),1/a))/(aa+bb)
    @staticmethod
    def ppf(a,b,aa,bb,q):
        return -math.log(1-math.pow(1-math.pow(1-q,1/b),1/a))/(aa+bb)

'''
Kumaraswamy generalized half-normal distribution
'''
class kumaraswamygenhalfnormal(Distribution): #kumaraswamy generalized half-normal
    @staticmethod
    def random(a,b,aa,t):
        return t*math.pow(st.norm.ppf((1+math.pow(1-math.pow(1-rg0(),1/b),1/a))/2),1/aa)
    @staticmethod
    def pdf(a,b,aa,t,x):
        pp=math.erf(math.pow(x/t,aa)/math.sqrt(2))
        return a*b*math.sqrt(2/math.pi)*(aa/x)*math.pow(x/t,aa)*math.exp(-1/2*math.pow(x/t,2*aa))*math.pow(pp-1,a-1)*math.pow(1-math.pow(pp-1,a),b-1)
    @staticmethod
    def cdf(a,b,aa,t,x):
        return 1-math.pow(1-math.pow(math.erf(math.pow(x/t,aa)/math.sqrt(2)),a),b)
    @staticmethod
    def median(a,b,aa,t):
        return t*math.pow(st.norm.ppf((1+math.pow(1-math.pow(1/2,1/b),1/a))/2),1/aa)
    @staticmethod
    def ppf(a,b,aa,t,q):
        return t*math.pow(st.norm.ppf((1+math.pow(1-math.pow(1-q,1/b),1/a))/2),1/aa)

'''
Kumaraswamy generalized power-Weibull distribution
'''
class kumaraswamygenpowerweibull(Distribution):#kumaraswamy generalized power-weibull
    @staticmethod
    def random(a,b,t,aa,l):
        return l*math.pow(math.pow(1-math.ln(1-math.pow(1-math.pow(1-rg0(),1/b),1/a)),1/t)-1,1/aa)
    @staticmethod
    def pdf(a,b,t,aa,l,x):
        p5=math.pow(1-math.pow(1-math.exp(1-math.pow(1+math.pow(x/l,aa),t)),a),b-1)
        p4=math.pow(1-math.exp(1-math.pow(1+math.pow(x/l,aa),t)),a-1)
        p3=math.exp(1-math.pow(1+math.pow(x/l,aa),t))
        p2=math.pow(1+math.pow(x/l,aa),t-1)
        p1=a*b*aa*math.pow(x,aa-1)/math.pow(l,aa)
        return p1*p2*p3*p4*p5
    @staticmethod
    def cdf(a,b,t,aa,l,x):
        return 1-math.pow(1-math.pow(1-math.exp(1-math.pow(1+math.pow(x/l,aa),t)),a),b)
    @staticmethod
    def median(a,b,t,aa,l):
        return l*math.pow(math.pow(1-math.ln(1-math.pow(1-math.pow(1/2,1/b),1/a)),1/t)-1,1/aa)
    @staticmethod
    def ppf(a,b,t,aa,l,q):
        return l*math.pow(math.pow(1-math.ln(1-math.pow(1-math.pow(1-q,1/b),1/a)),1/t)-1,1/aa)

'''
Kumaraswamy-Gumbel distribution
'''
class kumaraswamygumbel(Distribution):
    @staticmethod
    def random(a,b,m,s):
        return -s*math.log(-math.log(math.pow(1-math.pow(1-rg0(),1/b),1/a)))+m
    @staticmethod
    def median(a,b,m,s):
        return -s*math.log(-math.log(math.pow(1-math.pow(1/2,1/b),1/a)))+m
    @staticmethod
    def ppf(a,b,m,s,q):
        return -s*math.log(-math.log(math.pow(1-math.pow(1-q,1/b),1/a)))+m

'''
Kumaraswamy half-Cauchy distribution
'''
class kumaraswamyhalfcauchy(Distribution):
    @staticmethod
    def random(a,b,d):
        return d*math.tan(math.pi/2*math.pow(1-math.pow(1-rg0(),1/b),1/a))
    @staticmethod
    def pdf(a,b,d,x):
        return a*b*math.pow(2/math.pi,a)/d*1/(1+(x/d)**2)*math.pow(math.atan(x/d),a-1)*math.pow(1-math.pow(2/math.pi*math.atan(x/d),a),b-1)
    @staticmethod
    def cdf(a,b,d,x):
        return 1-math.pow(1-math.pow(2/math.pi*math.atan(x/d),a),b)
    @staticmethod
    def median(a,b,d):
        return d*math.tan(math.pi/2*math.pow(1-math.pow(1/2,1/b),1/a))
    @staticmethod
    def ppf(a,b,d,q):
        return d*math.tan(math.pi/2*math.pow(1-math.pow(1-q,1/b),1/a))

'''
Kumaraswamy inverse exponential distribution
'''
class kumaraswamyinvexp(Distribution): #kumaraswamy inverse exponential
    @staticmethod
    def random(a,b,l):
        return -a*l/math.log(1-math.pow(1-rg0(),1/b))
    @staticmethod
    def pdf(a,b,l,x):
        return a*b*l/x**2*math.exp(-a*l/x)*math.pow(1-math.exp(-a*l/x),b-1)
    @staticmethod
    def cdf(a,b,l,x):
        return 1-math.pow(1-math.exp(-a*l/x),b)
    @staticmethod
    def median(a,b,l):
        return -a*l/math.log(1-math.pow(1/2,1/b))
    @staticmethod
    def ppf(a,b,l,q):
        return -a*l/math.log(1-math.pow(1-q,1/b))

'''
Kumaraswamy inverse Weibull distribution
'''
class kumaraswamyinvweibull(Distribution): #kumaraswamy inverse weibull
    @staticmethod
    def random(a,aa,b,bb):
        return math.pow(-a*aa/math.log(1-math.pow(1-rg0(),1/b)),1/bb)
    @staticmethod
    def pdf(a,aa,b,bb,x):
        return a*b*aa*bb/math.pow(x,bb+1)*math.exp(-aa/math.pow(x,bb))*math.pow(math.exp(-a/math.pow(x,bb)),a-1)*math.pow(1-math.pow(math.exp(-aa/math.pow(x,bb)),a),b-1)
    @staticmethod
    def cdf(a,aa,b,bb,x):
        return 1-math.pow(1-math.exp(-a*aa/math.pow(x,bb)),b)
    @staticmethod
    def median(a,aa,b,bb):
        return math.pow(-a*aa/math.log(1-math.pow(1/2,1/b)),1/bb)
    @staticmethod
    def ppf(a,aa,b,bb,q):
        return math.pow(-a*aa/math.log(1-math.pow(1-q,1/b)),1/bb)

'''
Kumaraswamy-Kumaraswamy distribution
'''
class kumaraswamykumaraswamy(Distribution): #kumaraswamy-kumaraswamy
    @staticmethod
    def random(a,b,c,d):
        return math.pow(1-math.pow(1-math.pow(1-math.pow(1-rg0(),1/b),1/a),1/d),1/c)
    @staticmethod
    def median(a,b,c,d):
        return math.pow(1-math.pow(1-math.pow(1-math.pow(1/2,1/b),1/a),1/d),1/c)
    @staticmethod
    def ppf(a,b,c,d,q):
        return math.pow(1-math.pow(1-math.pow(1-math.pow(1-q,1/b),1/a),1/d),1/c)

'''
Kumaraswamy linear exponential distribution
'''
class kumaraswamylinearexp(Distribution): #kumaraswamy linear exponential
    @staticmethod
    def random(a,b,l,t):
        return (-l+math.sqrt(l**2-(2*t/a)*math.log(1-(1-math.pow(1-rg0(),1/b)))))/t
    @staticmethod
    def random(a,b,l,t):
        return (-l+math.sqrt(l**2-(2*t/a)*math.log(1-(1-math.pow(1/2,1/b)))))/t
    @staticmethod
    def random(a,b,l,t,q):
        return (-l+math.sqrt(l**2-(2*t/a)*math.log(1-(1-math.pow(1-q,1/b)))))/t

'''
Kumaraswamy log-logistic distribution
'''
class kumaraswamyloglogistic(Distribution):
    @staticmethod
    def random(a,b,aa,g):
        return aa*(math.pow(1-math.pow(1-math.pow(1-rg0(),1/b),1/a),-1/g)-1)
    @staticmethod
    def pdf(a,b,aa,g,x):
        return a*b*g/math.pow(aa,a*g)*math.pow(x,a*g-1)*math.pow(1+math.pow(x/aa,g),-a-1)*math.pow(1-math.pow(1-1/(1+(x/aa)*g),a),b-1)
    @staticmethod
    def cdf(a,b,aa,g,x):
        return 1-math.pow(1-math.pow(1-1/(1+(x/aa)*g),a),b)
    @staticmethod
    def median(a,b,aa,g):
        return aa*(math.pow(1-math.pow(1-math.pow(1/2,1/b),1/a),-1/g)-1)
    @staticmethod
    def ppf(a,b,aa,g,q):
        return aa*(math.pow(1-math.pow(1-math.pow(1-q,1/b),1/a),-1/g)-1)

'''
Kumaraswamy-Pareto distribution
'''
class kumaraswamypareto(Distribution):
    @staticmethod
    def random(a,b,bb,k):
        return bb/math.pow(1-math.pow(1-math.pow(1-rg0(),1/b),1/a),1/k)
    @staticmethod
    def pdf(a,b,bb,k,x):
        return a*b*k*math.pow(bb,k)/math.pow(x,k+1)*math.pow(1-math.pow(bb/x,k),a-1)*math.pow(1-math.pow(1-math.pow(bb/x,k),a),b-1)
    @staticmethod
    def cdf(a,b,bb,k,x):
        return 1-math.pow(1-math.pow(1-math.pow(bb/x,k),a),b)
    @staticmethod
    def median(a,b,bb,k):
        return bb/math.pow(1-math.pow(1-math.pow(1/2,1/b),1/a),1/k)
    @staticmethod
    def ppf(a,b,bb,k,q):
        return bb/math.pow(1-math.pow(1-math.pow(1-q,1/b),1/a),1/k)

'''
Kumaraswamy Weibull-Poisson distribution
'''
class kumaraswamyweibullpoisson(Distribution):
    @staticmethod
    def random(a,b,bb,c,l):
        u=rg0()
        return math.pow(-math.log(1+math.log(1-u*(1-math.exp(-l)))/l),1/c)/bb
    @staticmethod
    def pdf(a,b,bb,c,l,x):
        return 1-math.exp(-l*(1-math.exp(-math.pow(bb*x,c))))/(1-math.exp(-l))
    @staticmethod
    def cdf(a,b,bb,c,l,x):
        z=1-math.exp(-math.pow(bb*x,c))
        return l*a*b*c*math.pow(bb,c)/(math.exp(l)-1)*math.pow(x,c-1)*math.pow(z,a-1)*math.pow(1-math.pow(z,a),b-1)*math.exp(l*math.pow(1-math.pow(z,a),b)-math.pow(bb*x,c))
    @staticmethod
    def median(a,b,bb,c,l):
        return math.pow(-math.log(1+math.log(1-(1-math.exp(-l))/2)/l),1/c)/bb
    @staticmethod
    def ppf(a,b,bb,c,l,q):
        return math.pow(-math.log(1+math.log(1-q*(1-math.exp(-l)))/l),1/c)/bb

'''
Laha distribution*
'''
class laha(Distribution):
    @staticmethod
    def random(a,s):
        n=r.random()
        if(n<0.5):
            return -halflaha.random(a,s)
        return halflaha.random(a,s)

'''
Lambda-prime distribution
'''
class lambdaprime(Distribution):
    @staticmethod
    def random(t,nu):
        return normal.random(0,1)+t*math.sqrt(chi2.random(nu)/nu)

'''
Laplace distribution
'''
class laplace(Distribution):
    @staticmethod
    def pdf(mu,b,x):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        return 1/(2*b)*math.exp(-math.abs(x-mu)/b)
    @staticmethod
    def cdf(mu,b,x):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        if(x<mu):
            return 1/2*math.exp((x-mu)/b)
        return 1-1/2*math.exp(-(x-mu)/b)
    @staticmethod
    def random(mu,b):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        u=r.random()
        return mu-b*(math.abs(u-0.5)/(u-0.5))*math.log(1-2*math.abs(u-0.5))
    @staticmethod
    def kurtosis(mu,b):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        return 3
    @staticmethod
    def mean(mu,b):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        return mu
    @staticmethod
    def median(mu,b):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        return mu
    @staticmethod
    def mode(mu,b):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        return mu
    @staticmethod
    def variance(mu,b):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        return 2*(b**2)
    @staticmethod
    def stddev(mu,b):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        return b*math.sqrt(2)
    @staticmethod
    def entropy(mu,b):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        return math.log(2*b*math.exp(1))
    @staticmethod
    def skewness(mu,b):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        return 0
    @staticmethod
    def ppf(mu,b,q):
        if(b<=0):
            raise InvalidInputError("b must be positive")
        return mu-b*(math.abs(q-0.5)/(q-0.5))*math.log(1-2*math.abs(q-0.5))
    @staticmethod
    def mle(x):
        mu=np.median(x)
        b=np.average(np.abs(np.array(x)-mu))/len(x)
        return {'mu':mu,'b':b}

'''
Levy distribution
'''
class levy(Distribution):
    @staticmethod
    def pdf(mu,c,x):
        if(c<=0):
            raise InvalidInputError("c must be positive")
        if(x<mu):
            raise InvalidInputError("x must be greater than or equal to mu")
        return math.sqrt(c/(2*math.pi))*math.exp((-c)/(2*(x-mu)))/math.pow(x-mu,3/2)
    @staticmethod
    def cdf(mu,c,x):
        if(c<=0):
            raise InvalidInputError("c must be positive")
        if(x<mu):
            raise InvalidInputError("x must be greater than or equal to mu")
        return math.erfc(math.sqrt(c/(2*(x-mu))))
    @staticmethod
    def random(mu,c):
        if(c<=0):
            raise InvalidInputError("c must be positive")
        return c/math.pow(math.abs(st.norm.ppf(1-rg0()/2)),2)+mu
    @staticmethod
    def kurtosis(mu,c):
        if(c<=0):
            raise InvalidInputError("c must be positive")
        return None
    @staticmethod
    def mean(mu,c):
        if(c<=0):
            raise InvalidInputError("c must be positive")
        return float("infinity")
    @staticmethod
    def variance(mu,c):
        if(c<=0):
            raise InvalidInputError("c must be positive")
        return float("infinity")
    @staticmethod
    def stddev(mu,c):
        if(c<=0):
            raise InvalidInputError("c must be positive")
        return float("infinity")
    @staticmethod
    def entropy(mu,c):
        if(c<=0):
            raise InvalidInputError("c must be positive")
        return (1+3*np.euler_gamma+math.log(16*math.pi*c*c))/2
    @staticmethod
    def skewness(mu,c):
        if(c<=0):
            raise InvalidInputError("c must be positive")
        return None
    @staticmethod
    def ppf(mu,c):
        if(c<=0):
            raise InvalidInputError("c must be positive")
        return c/math.pow(math.abs(st.norm.ppf(1-q/2)),2)+mu

'''
Levy-l distribution
'''
class levyl(Distribution):
    @staticmethod
    def random():
        return st.levy_l.rvs()

'''
Levy stable distribution
'''
class levystable(Distribution):
    @staticmethod
    def random(a,b):
        return st.levy_stable.rvs(a,b)

'''
Log-Cauchy distribution
'''
class logcauchy(Distribution):
    @staticmethod
    def random(mu,sigma):
        return math.exp(mu+sigma*math.tan(math.pi*(r.random()-1/2)))
    @staticmethod
    def pdf(mu,sigma,x):
        return 1/(x*math.pi)*sigma/((math.log(x)-mu)**2+sigma**2)
    @staticmethod
    def cdf(mu,sigma,x):
        return 1/math.pi*math.atan((math.log(x)-mu)/sigma)+1/2
    @staticmethod
    def median(mu,sigma):
        return math.exp(mu)
    @staticmethod
    def variance(mu,sigma):
        return float("infinity")
    @staticmethod
    def ppf(mu,sigma,q):
        return math.exp(mu+sigma*math.tan(math.pi*(q-1/2)))

'''
Log-gamma distribution
'''
class loggamma(Distribution):
    @staticmethod
    def random(c):
        return st.loggamma.rvs(c)

'''
Log generalized Lindley-Weibull distribution
'''
class loggenlindleyweibull(Distribution): #log generalized lindley-weibull
    @staticmethod
    def random(a,b,t,g,c):
        u=rg0()
        v1=gamma.random(a,t)
        v2=gamma.random(a+1,t)
        if(u<=t/(b+t)):
            return g*math.pow(v1,1/c)
        return g*math.pow(v2,1/c)
    @staticmethod
    def pdf(a,b,t,g,c,x):
        return c*math.pow(t,a+1)/(g*(b+t)*math.gamma(a+1))*math.pow(x/g,c*a-1)*(a+b*math.pow(x/g,c))*math.exp(-t*math.pow(x/g,c))
    @staticmethod
    def mean(a,b,t,g,c):
        return g*math.pow(t,1-1/c)/((b+t)*math.gamma(a+1))*math.gamma(a+1/c)*(a+b/t*(a+1/c))
    @staticmethod
    def mode(a,b,t,g,c):
        return g*math.pow((c*(b-t)*a+b*(c-1))/(2*c*t*b)+math.sqrt(math.pow(-c*(b-t)*a-b*(c-1),2)-4*b*c*t*(-c*a*a+a))/(2*c*t*b),1/c)

'''
Logistic distribution
'''
class logistic(Distribution):
    @staticmethod
    def pdf(mu,s,x):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        z=-(x-mu)/s
        return math.exp(z)/(s*(1+math.exp(z))**2)
    @staticmethod
    def cdf(mu,s,x):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        z=-(x-mu)/s
        return 1/(1+math.exp(z))
    @staticmethod
    def random(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        n=rg0()
        while(n==1):
            n=rg0()
        return mu+s*math.log(n/(1-n))
    @staticmethod
    def kurtosis(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return 1.2
    @staticmethod
    def mean(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return mu
    @staticmethod
    def median(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return mu
    @staticmethod
    def mode(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return mu
    @staticmethod
    def variance(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return (s*math.pi)**2/3
    @staticmethod
    def stddev(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return s*math.pi/math.sqrt(3)
    @staticmethod
    def entropy(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return math.log(s)+2
    @staticmethod
    def skewness(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return 0
    @staticmethod
    def ppf(mu,s,q):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return mu+s*math.log(q/(1-q))

'''
Logistic Burr XII distribution
'''
class logisticburr12(Distribution):
    @staticmethod
    def random(c,k,l,s):
        return s*math.pow(math.exp(math.pow(math.pow(1/rg0()-1,-1/l),1/k))-1,1/c)
    @staticmethod
    def pdf(c,k,l,s,x):
        return l*c*k*math.pow(x,c-1)/(math.pow(s,c)*(1+math.pow(x/s,c)))*math.pow(math.log(math.pow(1+math.pow(x/s,c),k)),-l-1)*1/math.pow(1+math.pow(math.log(math.pow(1+math.pow(x/s,c),k)),-l),2)
    @staticmethod
    def cdf(c,k,l,s,x):
        return 1/(1+math.pow(math.log(math.pow(1+math.pow(x/s,c),k)),-l))
    @staticmethod
    def median(c,k,l,s):
        return s*math.pow(math.exp(math.pow(math.pow(1,-1/l),1/k))-1,1/c)
    @staticmethod
    def ppf(c,k,l,s,q):
        return s*math.pow(math.exp(math.pow(math.pow(1/q-1,-1/l),1/k))-1,1/c)

'''
Logistic exponential distribution
'''
class logisticexp(Distribution): #logistic exponential
    @staticmethod
    def random(b,l):
        return logisticweibull.random(1,b,l)
    @staticmethod
    def pdf(b,l,x):
        return logisticweibull.pdf(1,b,l)
    @staticmethod
    def cdf(b,l,x):
        return logisticweibull.cdf(1,b,l)
    @staticmethod
    def median(b,l):
        return logisticweibull.median(1,b,l)

'''
Logistic Frechet distribution
'''
class logisticfrechet(Distribution):
    @staticmethod
    def random(a,b,l):
        return -b/math.pow(math.log(1-math.exp(-math.pow(1/rg0()-1,-1/l))),1/a)
    @staticmethod
    def pdf(a,b,l,x):
        return l*a*math.pow(b,a)*math.pow(x,-a-1)*math.exp(-math.pow(b/x,a))/(1-math.exp(math.pow(-b/x,a)))*math.pow(-math.log(1-math.exp(math.pow(-b/x,a))),-l-1)*math.pow(1+math.pow(-math.log(1-math.exp(math.pow(-b/x,a))),-l),-2)
    @staticmethod
    def cdf(a,b,l,x):
        return 1/(1+math.pow(-math.log(1-math.exp(math.pow(-b/x,a))),-l))
    @staticmethod
    def median(a,b,l):
        return -b/math.pow(math.log(1-math.exp(-math.pow(1,-1/l))),1/a)
    @staticmethod
    def ppf(a,b,l,q):
        return -b/math.pow(math.log(1-math.exp(-math.pow(1/q-1,-1/l))),1/a)

'''
Logistic log-logistic distribution
'''
class logisticloglogistic(Distribution):
    @staticmethod
    def random(c,l,s):
        return logisticburr12.random(c,1,l,s)
    @staticmethod
    def pdf(c,l,s,x):
        return logisticburr12.pdf(c,1,l,s,x)
    @staticmethod
    def cdf(c,l,s,x):
        return logisticburr12.cdf(c,1,l,s,x)
    @staticmethod
    def median(c,l,s):
        return logisticburr12.median(c,1,l,s)

'''
Logistic Lomax distribution
'''
class logisticlomax(Distribution):
    @staticmethod
    def random(k,l,s):
        return logisticburr12.random(1,k,l,s)
    @staticmethod
    def pdf(k,l,s,x):
        return logisticburr12.pdf(1,k,l,s,x)
    @staticmethod
    def cdf(k,l,s,x):
        return logisticburr12.cdf(1,k,l,s,x)
    @staticmethod
    def median(k,l,s):
        return logisticburr12.median(1,k,l,s)

'''
Logistic Pareto distribution
'''
class logisticpareto(Distribution):
    @staticmethod
    def random(k,l,t):
        return t*math.pow(math.exp(math.pow(1/rg0()-1,-1/l)),1/k)
    @staticmethod
    def pdf(k,l,t,x):
        return l*k/x*math.pow(math.log(math.pow(x/t,k)),-l-1)/math.pow(1+math.pow(math.log(math.pow(x/t,k)),-l),2)
    @staticmethod
    def cdf(k,l,t,x):
        return 1/(1+math.pow(math.log(math.pow(x/t,k)),-l))
    @staticmethod
    def median(k,l,t):
        return t*math.pow(math.exp(math.pow(1,-1/l)),1/k)
    @staticmethod
    def ppf(k,l,t,q):
        return t*math.pow(math.exp(math.pow(1/q-1,-1/l)),1/k)

'''
Logistic Rayleigh distribution
'''
class logisticrayleigh(Distribution):
    @staticmethod
    def random(b,l):
        return logisticweibull.random(2,b,l)
    @staticmethod
    def pdf(b,l,x):
        return logisticweibull.pdf(2,b,l)
    @staticmethod
    def cdf(b,l,x):
        return logisticweibull.cdf(2,b,l)
    @staticmethod
    def median(b,l):
        return logisticweibull.median(2,b,l)

'''
Logistic uniform distribution
'''
class logisticuniform(Distribution):
    @staticmethod
    def random(l,t):
        return t*(1-math.exp(1-math.pow(1/rg0(),-1/l)))
    @staticmethod
    def pdf(l,t,x):
        return l/(t-x)*math.pow(math.log(t/(t-x)),-l-1)*math.pow(1+math.pow(math.log(t/(t-x)),-l),-2)
    @staticmethod
    def cdf(l,t,x):
        return 1/(1+math.pow(math.log(t/(t-x)),-l))
    @staticmethod
    def median(l,t):
        return t*(1-math.exp(1-math.pow(2,-1/l)))
    @staticmethod
    def ppf(l,t,q):
        return t*(1-math.exp(1-math.pow(1/q,-1/l)))

'''
Logistic Weibull distribution
'''
class logisticweibull(Distribution):
    @staticmethod
    def random(a,b,l):
        return math.pow((math.pow(1/rg0(),-1/l)-1)/b,1/a)
    @staticmethod
    def pdf(a,b,l,x):
        return l*a/math.pow(b,l)*math.pow(x,-l*a-1)*1/math.pow(1+math.pow(b*math.pow(x,a),-l),2)
    @staticmethod
    def cdf(a,b,l,x):
        return 1/(1+math.pow(b*math.pow(x,a),-l))
    @staticmethod
    def median(a,b,l):
        return math.pow((math.pow(2,-1/l)-1)/b,1/a)
    @staticmethod
    def ppf(a,b,l,q):
        return math.pow((math.pow(1/q,-1/l)-1)/b,1/a)

'''
Log-Laplace distribution
'''
class loglaplace(Distribution):
    @staticmethod
    def random(mu,b):
        return math.exp(laplace.random(mu,b))
    @staticmethod
    def pdf(mu,b,x):
        return 1/(2*b*x)*math.exp(-math.abs(math.log(x)-mu)/b)
    @staticmethod
    def cdf(mu,b,x):
        if(x>0):
            return (1+(1-math.exp(-math.abs(math.log(x)-mu)/b))*(math.abs(math.log(x)-mu)/(math.log(x)-mu)))

'''
Log-logistic distribution
'''
class loglogistic(Distribution):
    @staticmethod
    def random(aa,bb):
        n=r.random()
        while(n==1):
            n=r.random()
        return aa*math.pow((n/(1-n)),1/bb)
    @staticmethod
    def pdf(aa,bb,x):
        return (bb/aa)*math.pow(x/aa,bb-1)/math.pow(1+math.pow(x/aa,bb),2)
    @staticmethod
    def cdf(aa,bb,x):
        return 1/(1+math.pow(x/aa,-bb))
    @staticmethod
    def mean(aa,bb):
        if(bb>1):
            return aa*math.pi/bb/(math.sin(math.pi/bb))
        return None
    @staticmethod
    def median(aa,bb):
        return aa
    @staticmethod
    def mode(aa,bb):
        if(b>1):
            return aa*math.pow((bb-1)/(bb+1),1/bb)
        return 0
    @staticmethod
    def variance(aa,bb):
        if(bb>2):
            b=math.pi/bb
            return aa**2*(2*b/math.sin(2*b)-b**2/(math.sin(b)**2))
    @staticmethod
    def stddev(aa,bb):
        return aa*math.sqrt(2*b/math.sin(2*b)-b**2/(math.sin(b)**2))
    @staticmethod
    def ppf(aa,bb,q):
        return aa*math.pow((q/(1-q)),1/bb)

'''
Log-normal distribution
'''
class lognormal(Distribution):
    @staticmethod
    def random(mu,sigma):
        return math.exp(normal.random(mu,sigma))
    @staticmethod
    def pdf(mu,sigma,x):
        return 1/(x*sigma*math.sqrt(2*math.pi))*math.exp(-(math.log(x)-mu)**2/(2*sigma**2))
    @staticmethod
    def cdf(mu,sigma,x):
        return 1/2+1/2*math.erf((math.log(x)-mu)/math.sqrt(2)*sigma)
    @staticmethod
    def kurtosis(mu,sigma):
        return math.exp(4*sigma**2)+2*math.exp(3*sigma**2)+3*math.exp(2*sigma**2)-6
    @staticmethod
    def mean(mu,sigma):
        return math.exp(mu+sigma**2/2)
    @staticmethod
    def median(mu,sigma):
        return math.exp(mu)
    @staticmethod
    def mode(mu,sigma):
        return math.exp(mu-sigma**2)
    @staticmethod
    def variance(mu,sigma):
        return (math.exp(sigma**2)-1)*math.exp(2*mu+sigma**2)
    @staticmethod
    def stddev(mu,sigma):
        return math.sqrt(math.exp(sigma**2)-1)*math.exp(2*mu+sigma**2)
    @staticmethod
    def entropy(mu,sigma):
        return math.log(sigma*math.exp(mu+1/2)*math.sqrt(2*math.pi))
    @staticmethod
    def skewness(mu,sigma):
        return (math.exp(sigma**2)+2)*math.sqrt(math.exp(sigma**2)-1)

'''
Log-triangular distribution
'''
class logtriangular(Distribution):
    @staticmethod
    def random(a,b,c):
        n=uniform.random(a,c)
        if(math.log(n)<=b):
            return a*math.exp(math.sqrt(n*(math.log(b)-math.log(a))*(math.log(c)-math.log(a))))
        return c*math.exp(-math.sqrt(n*(math.log(c)-math.log(b))*(math.log(c)-math.log(a))))
    @staticmethod
    def pdf(a,b,c,x):
        if(a<x and x<=c):
            return 2*math.log(x/a)/(math.log(b/a)*math.log(c/a))
        if(c<x and x<=b):
            return 2*math.log(b/x)/(math.log(b/a)*math.log(b/c))
    @staticmethod
    def median(a,b,c):
        if(math.log(1/2)<=b):
            return a*math.exp(math.sqrt((math.log(b)-math.log(a))/2*(math.log(c)-math.log(a))))
        return c*math.exp(-math.sqrt((math.log(c)-math.log(b))/2*(math.log(c)-math.log(a))))
    @staticmethod
    def ppf(a,b,c,q):
        if(math.log(q)<=b):
            return a*math.exp(math.sqrt(q*(math.log(b)-math.log(a))*(math.log(c)-math.log(a))))
        return c*math.exp(-math.sqrt(q*(math.log(c)-math.log(b))*(math.log(c)-math.log(a))))

'''
Lomax distribution
'''
class lomax(Distribution):
    @staticmethod
    def pdf(lmbda,aa,x):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        if(x<0):
            raise InvalidInputError("x must be non-negative")
        return aa/lmbda*math.pow(1+x/lmbda,-aa-1)
    @staticmethod
    def cdf(lmbda,aa,x):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        if(x<0):
            raise InvalidInputError("x must be non-negative")
        return 1-math.pow(1+x/lmbda,-aa)
    @staticmethod
    def random(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        n=r.random()
        while(n==1):
            n=r.random()
        return lmbda*(math.pow(1/(1-n),1/aa)-1)
    @staticmethod
    def kurtosis(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        if(aa>4):
            return 6*(aa**3+aa**2-6*aa-2)/(aa*(aa-3)*(aa-4))
        return None
    @staticmethod
    def mean(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        if(aa>1):
            return lmbda/(aa-1)
        return None
    @staticmethod
    def median(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        return lmbda*(math.pow(2,1/aa)-1)
    @staticmethod
    def mode(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        return 0
    @staticmethod
    def variance(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        if(aa>2):
            return lmbda**2*aa/((aa-1)**2*(aa-2))
        if(aa>1):
            return float("infinity")
        return None
    @staticmethod
    def stddev(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        if(aa>2):
            return lmbda*math.sqrt(aa/(aa-2))/(aa-1)
        if(aa>1):
            return float("infinity")
        return None
    @staticmethod
    def skewness(lmbda,aa):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        if(aa>3):
            return 2*(1+aa)/(aa-3)*math.sqrt((aa-2)/aa)
        return None
    @staticmethod
    def ppf(lmbda,aa,q):
        if(lmbda<=0 or aa<=0):
            raise InvalidInputError("lambda and aa must be positive")
        return lmbda*(math.pow(1/(1-q),1/aa)-1)

'''
Maxwell-Boltzmann distribution
'''
class maxwellboltzmann(Distribution):
    @staticmethod
    def pdf(aa,x):
        if(aa<=0 or x<=0):
            raise InvalidInputError("all inputs must be positive")
        return math.sqrt(2/math.pi)*(x**2*math.exp(x**2/(2*aa**2)))/(aa**3)
    @staticmethod
    def cdf(aa,x):
        if(aa<=0 or x<=0):
            raise InvalidInputError("all inputs must be positive")
        return math.erf(x/(math.sqrt(2)*aa))-math.sqrt(2/math.pi)*(x*math.exp(x**2/(2*aa**2)))/aa
    @staticmethod
    def random(aa):
        if(aa<=0):
            raise InvalidInputError("aa must be positive")
        r=(-math.log(rg0()))**2
        w1=(rg0())**2
        w2=(rg0())**2
        while(w1+w2>1):
            w1=(rg0())**2
            w2=(rg0())**2
        r2=r-(w1/(w1+w2))*math.log(rg0())
        return aa*math.sqrt(2*r2)
    @staticmethod
    def kurtosis(aa):
        if(aa<=0):
            raise InvalidInputError("aa must be positive")
        return 4*(-96+40*math.pi-3*(math.pi)**2)/((3*math.pi-8)**2)
    @staticmethod
    def mean(aa):
        if(aa<=0):
            raise InvalidInputError("aa must be positive")
        return 2*aa*math.sqrt(2/math.pi)
    @staticmethod
    def mode(aa):
        if(aa<=0):
            raise InvalidInputError("aa must be positive")
        return math.sqrt(2)*aa
    @staticmethod
    def variance(aa):
        if(aa<=0):
            raise InvalidInputError("aa must be positive")
        return aa**2*(3*math.pi-8)/math.pi
    @staticmethod
    def stddev(aa):
        if(aa<=0):
            raise InvalidInputError("aa must be positive")
        return aa*math.sqrt((3*math.pi-8)/math.pi)
    @staticmethod
    def entropy(aa):
        if(aa<=0):
            raise InvalidInputError("aa must be positive")
        return math.log(aa*math.sqrt(2*math.pi))+np.euler_gamma-1/2
    @staticmethod
    def skewness(aa):
        if(aa<=0):
            raise InvalidInputError("aa must be positive")
        return 2*math.sqrt(2)*(16-5*math.pi)/math.pow(3*math.pi-8,3/2)

'''
Marshall-Olkin Esscher transformed Laplace distribution
'''
class marshallolkinesschertransformedlaplace(Distribution):
    @staticmethod
    def random(b,t):
        w=exp.random(1/b)
        return 2*t/(1-t**2)*w+math.sqrt(2/(1-t**2))*math.sqrt(w)*normal.random(0,1)
    @staticmethod
    def pdf(b,t,x):
        l=math.sqrt(b*(1-t**2))
        k=l/(t+math.sqrt(l+t**2))
        if(x<0):
            return l*k/(1+k**2)*math.exp(l*x/k)
        return l*k/(1+k**2)*math.exp(-l*k*x)
    @staticmethod
    def cdf(b,t,x):
        l=math.sqrt(b*(1-t**2))
        k=l/(t+math.sqrt(l+t**2))
        if(x<0):
            return k**2/(1+k**2)*math.exp(l*x/k)
        return 1-1/(1+k**2)*math.exp(-l*k*x)
    @staticmethod
    def mean(b,t):
        l=math.sqrt(b*(1-t**2))
        k=l/(t+math.sqrt(l+t**2))
        return (1-k**2)/(l*k)
    @staticmethod
    def variance(b,t):
        l=math.sqrt(b*(1-t**2))
        k=l/(t+math.sqrt(l+t**2))
        return (1+k**4)/(l**2*k**2)
    @staticmethod
    def stddev(b,t):
        l=math.sqrt(b*(1-t**2))
        k=l/(t+math.sqrt(l+t**2))
        return math.sqrt((1+k**4)/(l*k))

'''
McDonald log-logistic distribution
'''
class mcdonaldloglogistic(Distribution):
    @staticmethod
    def random(a,b,c,aa,bb):
        u=r.random()
        while(u==1):
            u=r.random()
        q1=aa*math.pow((u/(1-u)),1/bb) #Q_alpha,beta
        q2=a/c*math.pow((u/(1-u)),1/b) #Q_{a/c},b
        return q1*math.pow(q2,1/c)/(1-q1*math.pow(q2,1/c))
    @staticmethod
    def pdf(a,b,c,aa,bb):
        return c/sp.beta(a/c,b)*aa/bb*math.pow(x/bb,a*aa-1)*math.pow(1+math.pow(x/bb,aa),-a-1)*math.pow(1-math.pow(1-1/(1+math.pow(x/bb,aa)),c),b-1)
    @staticmethod
    def median(a,b,c,aa,bb):
        q1=aa #Q_alpha,beta
        q2=a/c #Q_{a/c},b
        return q1*math.pow(q2,1/c)/(1-q1*math.pow(q2,1/c))
    @staticmethod
    def ppf(a,b,c,aa,bb,q):
        q1=aa*math.pow((q/(1-q)),1/bb) #Q_alpha,beta
        q2=a/c*math.pow((q/(1-q)),1/b) #Q_{a/c},b
        return q1*math.pow(q2,1/c)/(1-q1*math.pow(q2,1/c))

'''
Meridian distribution*
'''
class meridian(Distribution):
    @staticmethod
    def random(a,s):
        return genpearson7.random(a,s,2,1)

'''
Mielke beta-kappa distribution
'''
class mielkebetakappa(Distribution):
    @staticmethod
    def random(k,theta):
        n=rg0()
        while n==1:
            n=rg0()
        return math.pow((math.pow(n,theta/k))/(1-math.pow(n,theta/k)),1/theta)
    @staticmethod
    def pdf(k,theta,x):
        return k*math.pow(x,k-1)/math.pow(1+math.pow(x,theta),1+k/theta)
    @staticmethod
    def cdf(k,theta,x):
        return math.pow(x,k)/math.pow(1+math.pow(x,theta),k/theta)
    @staticmethod
    def median(k,theta):
        return math.pow((math.pow(1/2,theta/k))/(1-math.pow(1/2,theta/k)),1/theta)
    @staticmethod
    def ppf(k,theta,q):
        return math.pow((math.pow(q,theta/k))/(1-math.pow(q,theta/k)),1/theta)

'''
Minimax distribution
'''
class minimax(Distribution):
    @staticmethod
    def random(aa,bb):
        return math.pow(1-math.pow(1-r.random(),1/aa),1/bb)
    @staticmethod
    def pdf(aa,bb,x):
        return aa*bb*math.pow(x,aa-1)*math.pow(1-math.pow(x,aa),bb-1)
    @staticmethod
    def cdf(aa,bb,x):
        return 1-math.pow(1-math.pow(x,aa),bb)
    @staticmethod
    def median(aa,bb):
        return math.pow(1-math.pow(1/2,1/aa),1/bb)
    @staticmethod
    def ppf(aa,bb,q):
        return math.pow(1-math.pow(1-q,1/aa),1/bb)

'''
Modified Burr III distribution
'''
class modburr3(Distribution):
    @staticmethod
    def random(a,b,g):
        return math.pow((math.pow(rg0(),-g/a)-1)/g,-1/b)
    @staticmethod
    def pdf(a,b,g,x):
        return a*b*math.pow(x,-b-1)*math.pow(1+g*math.pow(x,-b),-a/g-1)
    @staticmethod
    def cdf(a,b,g,x):
        return math.pow(1+g*math.pow(x,-b),-a/g)
    @staticmethod
    def median(a,b,g):
        return math.pow((math.pow(1/2,-g/a)-1)/g,-1/b)
    @staticmethod
    def ppf(a,b,g,q):
        return math.pow((math.pow(q,-g/a)-1)/g,-1/b)

'''
Modified Burr III Burr XII distribution
'''
class modburr3burr12(Distribution):
    @staticmethod
    def random(a,b,c,g,k):
        u=rg0()
        return math.pow(math.pow(math.pow((math.pow(u,-g/a)-1)/g,-1/b)+1,1/k)-1,1/c)
    @staticmethod
    def pdf(a,b,c,g,k,x):
        return a*b*c*k*math.pow(x,c-1)*math.pow(1+math.pow(x,c),k-1)*math.pow(math.pow(1+math.pow(x,c),k)-1,-b-1)*math.pow(1+g*math.pow(math.pow(1+math.pow(x,c),k)-1,-b),-a/g-1)
    @staticmethod
    def cdf(a,b,c,g,k,x):
        return math.pow(1+g*math.pow(math.pow(1+math.pow(x,c),k)-1,-b),-a/g)
    @staticmethod
    def median(a,b,c,g,k):
        return math.pow(math.pow(math.pow((math.pow(1/2,-g/a)-1)/g,-1/b)+1,1/k)-1,1/c)
    @staticmethod
    def ppf(a,b,c,g,k,q):
        return math.pow(math.pow(math.pow((math.pow(q,-g/a)-1)/g,-1/b)+1,1/k)-1,1/c)

'''
Modified Burr III Kumaraswamy distribution
'''
class modburr3kumaraswamy(Distribution):
    @staticmethod
    def random(a,aa,bb,g):
        u=rg0()
        return math.pow(1-math.pow(math.pow((math.pow(u,-g/aa)-1)/g,-1/bb)+1,-1/b),1/a)
    @staticmethod
    def pdf(a,aa,b,bb,g,x):
        r=math.pow(1-math.pow(x,a),-b)-1
        return aa*bb*a*b*math.pow(x,a-1)*math.pow(1-math.pow(x,a),-b-1)*math.pow(r,-bb-1)*math.pow(1+g*math.pow(r,-bb),-aa/g-1)
    @staticmethod
    def cdf(a,aa,b,bb,g,x):
        return math.pow(1+g*math.pow(math.pow(1-math.pow(x,a),-b)-1,-bb),-aa/g)
    @staticmethod
    def median(a,aa,b,bb,g):
        return math.pow(1-math.pow(math.pow((math.pow(1/2,-g/aa)-1)/g,-1/bb)+1,-1/b),1/a)
    @staticmethod
    def ppf(a,aa,b,bb,g,q):
        return math.pow(1-math.pow(math.pow((math.pow(q,-g/aa)-1)/g,-1/bb)+1,-1/b),1/a)

'''
Modified Burr III Weibull distribution
'''
class modburr3weibull(Distribution):
    @staticmethod
    def random(a,b,g,l,n):
        return math.pow(math.log(math.pow(math.pow(rg0(),-g/a)-1,-1/b)/g+1)/l,1/n)
    @staticmethod
    def pdf(a,b,g,l,n,x):
        return a*b*math.pow(math.exp(l*math.pow(x,n)-1),-b-1)*l*v*math.pow(x,n-1)*math.exp(l*math.pow(x,n))*math.pow(1+g*math.pow(math.exp(l*math.pow(x,n))-1,-b),-a/g-1)
    @staticmethod
    def median(a,b,g,l,n):
        return math.pow(math.log(math.pow(math.pow(1/2,-g/a)-1,-1/b)/g+1)/l,1/n)
    @staticmethod
    def ppf(a,b,g,l,n,q):
        return math.pow(math.log(math.pow(math.pow(q,-g/a)-1,-1/b)/g+1)/l,1/n)

'''
Modified extended generalized exponential distribution
'''
class modifextgenexp(Distribution): #modified extended generalized exponential
    @staticmethod
    def random(a,b,e,l):
        p=rg0()
        return -1/l*math.log(1/e*(b-math.pow(p*math.pow(b,a/e)-(p-1)*math.pow(b-e,a/e)),e/a))
    @staticmethod
    def pdf(a,b,e,l,x):
        return a*l*math.pow(b-e*math.exp(-l*x),a/e-1)*math.exp(-l*x)/(math.pow(b,a/e)-math.pow(b-e,a/e))
    @staticmethod
    def cdf(a,b,e,l,x):
        return (math.pow(b-e*math.exp(-l*x),a/e)-math.pow(b-e,a/e))/(math.pow(b,a/e)-math.pow(b-e,a/e))
    @staticmethod
    def median(a,b,e,l):
        return -1/l*math.log((1/e)*(b-math.pow((math.pow(b,a/e)+math.pow(b+e,a/e))/2,e/a)))
    @staticmethod
    def mode(a,b,e,l):
        if(a>e):
            return -1/l*math.log(b/a)
    @staticmethod
    def ppf(a,b,e,l,q):
        return -1/l*math.log(1/e*(b-math.pow(q*math.pow(b,a/e)-(q-1)*math.pow(b-e,a/e)),e/a))

'''
Moffat distribution*
'''
class moffat(Distribution):
    @staticmethod
    def random(s,g):
        return genbetaprime.random(0,s,1,g,2)

'''
Moyal distribution*
'''
class moyal(Distribution):
    @staticmethod
    def random():
        x1=rg0()
        x2=rg0()
        y=math.pi*x1-math.pi/2
        h=x2*0.912
        z=math.tan(y)
        hy=1/math.sqrt(2*math.pi)*1/(math.cos(y)^2)*math.exp(-(math.tan(y)+math.exp(-math.tan(y)))/2)
        while(h>hy):
            x1=rg0()
            x2=rg0()
            y=math.pi*x1-math.pi/2
            h=x2*0.912
            z=math.tan(y)
            hy=1/math.sqrt(2*math.pi)*1/(math.cos(y)^2)*math.exp(-(math.tan(y)+math.exp(-math.tan(y)))/2)
        return z

'''
Multinomial distribution
'''
class multinomial(Distribution): #probs is the array of probabilities (does not have to sum to 1 - rescaling happens inside the function)
    @staticmethod
    def random(probs,n):
        amt=[]
        for _ in range(len(probs)):
            amt.append(0)
        probs2=probs/sum(probs)
        for _ in range(n):
            a=r.random()
            i=0
            while(a>sum(probs2[:i+1])):
                i+=1
            amt[i]+=1
        return amt

'''
Nakagami distribution
'''
class nakagami(Distribution):
    @staticmethod
    def random(m,omega):
        return math.sqrt(omega/(2*m))*chi.random(2*m)
    @staticmethod
    def pdf(m,omega,x):
        return 2*math.pow(m,m)/(math.gamma(m)*math.pow(omega,m))*math.pow(x,2*m-1)*math.exp(-m/omega*x**2)
    @staticmethod
    def mean(m,omega):
        return math.gamma(m+1/2)/math.gamma(m)*math.sqrt(omega/m)
    @staticmethod
    def mode(m,omega):
        return math.sqrt(1/2)*math.sqrt((2*m-1)*omega/m)
    @staticmethod
    def variance(m,omega):
        return omega*(1-1/m*(math.gamma(m+1/2)/math.gamma(m))**2)
    @staticmethod
    def stddev(m,omega):
        return math.sqrt(omega*(1-1/m*(math.gamma(m+1/2)/math.gamma(m))**2))

'''
Negative binomial distribution
'''
class negbin(Distribution):
    @staticmethod
    def random(r,p):
        if(r%1!=0 or r<=0):
            raise InvalidInputError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        num_failure=0
        num_success=0
        while(num_failure<k):
            if(r.random()<=p):
                num_success+=1
            else:
                num_failure+=1
        return num_success
    @staticmethod
    def kurtosis(r,p):
        if(r%1!=0 or r<=0):
            raise InvalidInputError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        return 6/r+(1-p)**2/(p*r)
    @staticmethod
    def mean(r,p):
        if(r%1!=0 or r<=0):
            raise InvalidInputError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        return p*r/(1-p)
    @staticmethod
    def mode(r,p):
        if(r%1!=0 or r<=0):
            raise InvalidInputError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        if(r>1):
            return math.floor(p*(r-1)/(1-p))
        return 0
    @staticmethod
    def variance(r,p):
        if(r%1!=0 or r<=0):
            raise InvalidInputError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        return p*r/((1-p)**2)
    @staticmethod
    def stddev(r,p):
        if(r%1!=0 or r<=0):
            raise InvalidInputError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        return math.sqrt(p*r)/(1-p)
    @staticmethod
    def skewness(r,p):
        if(r%1!=0 or r<=0):
            raise InvalidInputError("r must be a positive integer")
        if(p<=0 or p>=1):
            raise InvalidInputError("p must be between 0 and 1 exclusive")
        return (1+p)/math.sqrt(p*r)

'''
Negative hypergeometric distribution
'''
class neghypergeo(Distribution): #negative hypergeometric
    @staticmethod
    def random(R,N,K):
        num_success=0
        num_failure=0
        total_=0
        while(num_failure<R):
            if(r.random()<=(K-num_success)/(N-total_)):
                num_success+=1
            else:
                num_failure+=1
                total_+=1
        return num_success
    @staticmethod
    def mean(R,N,K):
        return R*K/(N-K+1)
    @staticmethod
    def variance(R,N,K):
        return R*(N+1)*K/((N-K+1)*(N-K+2))*(1-R/(N-K+1))
    @staticmethod
    def stddev(R,N,K):
        return math.sqrt(R*(N+1)*K/((N-K+1)*(N-K+2))*(1-R/(N-K+1)))

'''
Negative multinomial distribution
'''
class negmultinomial(Distribution): #negative multinomial
    @staticmethod
    def random(probs, n):
        amt=[]
        for _ in range(len(probs)):
            amt.append(0)
        probs2=probs/sum(probs)
        while(amt[0]<n):
            a=r.random()
            i=0
            while(a>sum(probs2[:i+1])):
                i+=1
            amt[i]+=1
        return amt

'''
Non-central beta distribution
'''
class noncentralbeta(Distribution):
    @staticmethod
    def random(m,n,mu,sigma):
        return(noncentralchi2.random(mu,sigma,m))/(noncentralchi2.random(mu,sigma,m)+chi2.random(n))

'''
Non-central chi-square distribution
'''
class noncentralchi2(Distribution): #non-central chi-squared
    @staticmethod
    def random(mu,sigma,nu):
        sum_=0
        for _ in range(nu):
            sum_+=normal.random(mu,sigma)
        return sum_

'''
Non-central f distribution
'''
class noncentralf(Distribution):
    @staticmethod
    def random(mu,sigma,m,n):
        return(noncentralchi2.random(mu,sigma,m)/m)/(chi2.random(n)/n)

'''
Non-central t distribution
'''
class noncentralt(Distribution):
    @staticmethod
    def random(mu,nu):
        return (normal.random(0,1)+mu)/math.sqrt(chi2.random(nu)/nu)

'''
Normal distribution
'''
class normal(Distribution):
    @staticmethod
    def pdf(mu,sigma2,x):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return 1/math.sqrt(2*math.pi*sigma2)*math.exp(-(x-mu)**2/(2*sigma2))
    @staticmethod
    def cdf(mu,sigma2,x):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return 1/2*(1+math.erf((x-mu)/math.sqrt(2*sigma2)))
    @staticmethod
    def random(mu,sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        a=rg0()
        b=rg0()
        while(a==1):
            a=rg0()
        while(b==1):
            b=rg0()
        return math.sqrt(-2*math.log(a))*math.cos(2*math.pi*b)*sigma+mu
    @staticmethod
    def kurtosis(mu,sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return 0
    @staticmethod
    def mean(mu,sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return mu
    @staticmethod
    def median(mu,sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return mu
    @staticmethod
    def mode(mu,sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return mu
    @staticmethod
    def variance(mu,sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return sigma2
    @staticmethod
    def stddev(mu,sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return math.sqrt(sigma2)
    @staticmethod
    def entropy(mu,sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return 1/2*math.log(2*math.pi*math.exp(1)*sigma2)
    @staticmethod
    def skewness(mu,sigma2):
        if(sigma2<=0):
            raise InvalidInputError("sigma2 must be positive")
        return 0
    @staticmethod
    def mle(x):
        return {'mu':avg(x),'sigma2':np.var(x)*(len(x)-1)/len(x)}

'''
Normal-gamma distribution*
'''
class normalgamma(Distribution):
    @staticmethod
    def random(aa,bb,mu,lmbda):
        n=gamma.random(aa,bb)
        return normal.random(mu,math.sqrt(1/(lmbda*n)))

'''
Normal-inverse gamma distribution
'''
class norminvgamma(Distribution): #normal-inverse gamma
    @staticmethod
    def random(aa,bb,mu,lmbda):
        n=invgamma.random(aa,bb)
        return normal.random(mu,math.sqrt(n/lmbda))

'''
Odd generalized exponential-Gompertz distribution
'''
class oddgenexpgompertz(Distribution): #odd generalized exponential-gompertz
    @staticmethod
    def random(a,b,c,l):
        return 1/c*math.log(1+c/l*math.log(1-1/a*math.log(1-math.pow(rg0(),1/b))))
    @staticmethod
    def pdf(a,b,c,l,x):
        return a*b*l*math.exp(c*x)*math.exp(l/c*(math.exp(c*x)-1))*math.exp(-a*(math.exp(l/c*(math.exp(c*x)-1))-1))*math.pow(1-math.exp(-a*(math.exp(l/c*(math.exp(c*x)-1))-1)),b-1)
    @staticmethod
    def cdf(a,b,c,l,x):
        return math.pow(1-math.exp(-a*(math.exp(l/c*(math.exp(c*x)-1))-1)),b)
    @staticmethod
    def median(a,b,c,l):
        return 1/c*math.log(1+c/l*math.log(1-1/a*math.log(1-math.pow(1/2,1/b))))
    @staticmethod
    def ppf(a,b,c,l,q):
        return 1/c*math.log(1+c/l*math.log(1-1/a*math.log(1-math.pow(q,1/b))))

'''
Odd generalized exponentiated linear failure rate distribution
'''
class oddgenexplinearfailurerate(Distribution): #odd generalized exponentiated linear failure rate
    @staticmethod
    def random(a,b,aa,bb):
        i=1+math.log(1/math.pow(1-math.pow(rg0(),1/bb),1/aa))
        return (-a+math.sqrt(a*a+2*b*math.log(i)))/b
    @staticmethod
    def pdf(a,b,aa,bb,x):
        aa*bb*(a+b*x)*math.exp(a*x+b/2*x**2)*math.exp(-aa*(math.exp(a*x+b/2*x**2-1)))*math.pow(1-math.exp(-aa*(math.exp(a*x+b/2*x**2-1))),bb-1)
    @staticmethod
    def cdf(a,b,aa,bb,x):
        return math.pow(1-math.exp(-aa*(math.exp(a*x+b/2*x**2-1))),bb)
    @staticmethod
    def median(a,b,aa,bb):
        i=1+math.log(1/math.pow(1-math.pow(1/2,1/bb),1/aa))
        return (-a+math.sqrt(a*a+2*b*math.log(i)))/b
    @staticmethod
    def ppf(a,b,aa,bb,q):
        i=1+math.log(1/math.pow(1-math.pow(q,1/bb),1/aa))
        return (-a+math.sqrt(a*a+2*b*math.log(i)))/b

'''
Odd generalized exponential log-logistic distribution
'''
class oddgenexploglog(Distribution): #odd generalized exponential log-logistic
    @staticmethod
    def random(g,l,s,t):
        return s*math.pow(-l*math.log(1-math.pow(rg0(),1/g)),1/t)
    @staticmethod
    def pdf(g,l,s,t,x):
        return g*t/(l*s)*math.pow(x/s,t-1)*math.pow(1-math.exp(-1/l*math.pow(x/s,t)),g-1)*math.exp(-1/l*math.pow(x/s,t))
    @staticmethod
    def cdf(g,l,s,t,x):
        return math.pow(1-math.exp(-1/l*math.pow(x/s,t)),g)
    @staticmethod
    def median(g,l,s,t):
        return s*math.pow(-l*math.log(1-math.pow(rg0(),1/g)),1/t)
    @staticmethod
    def mode(g,l,s,t):
        return s*math.pow((1-t)/t*l,1/t)
    @staticmethod
    def ppf(g,l,s,t,q):
        return s*math.pow(-l*math.log(1-math.pow(q,1/g)),1/t)

'''
Odds generalized exponential-exponential distribution
'''
class oddsgenexpexp(Distribution): #odds generalized exponential-exponential
    @staticmethod
    def pdf(l,t,x):
        return l*t*math.exp(x*t)*math.exp(-l*(math.exp(x*t)-1))
    @staticmethod
    def cdf(l,t,x):
        return 1-math.exp(-l*(math.exp(t*x)-1))
    @staticmethod
    def random(l,t):
        return math.log(1-math.log(1-rg0())/l)/t
    @staticmethod
    def median(l,t):
        return math.log(1-math.log(1/2)/l)/t
    @staticmethod
    def ppf(l,t,q):
        return math.log(1-math.log(1-q)/l)/t

'''
Paralogistic distribution*
'''
class paralogistic(Distribution):
    @staticmethod
    def random(bb):
        return genbetaprime.random(0,1,1,bb,bb)

'''
Pareto distribution
'''
class pareto(Distribution):
    @staticmethod
    def pdf(alpha,xm,x):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        if(x<xm):
            raise InvalidInputError("x must be greater than or equal to xm")
        return alpha*math.pow(xm,alpha)/math.pow(x,alpha+1)
    @staticmethod
    def cdf(alpha,xm,x):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        if(x<xm):
            raise InvalidInputError("x must be greater than or equal to xm")
        return 1-math.pow(xm/x,alpha)
    @staticmethod
    def random(alpha,xm):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        return xm/math.pow(rg0(),1/alpha)
    @staticmethod
    def kurtosis(alpha,xm):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        if(alpha>4):
            return 6*(alpha**3+alpha**2-6*alpha-2)/(alpha*(alpha-3)*(alpha-4))
        return None
    @staticmethod
    def mean(alpha,xm):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        if(alpha<=1):
            return float("infinity")
        return alpha*xm/(alpha-1)
    @staticmethod
    def median(alpha,xm):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        return xm*math.pow(2,1/alpha)
    @staticmethod
    def mode(alpha,xm):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        return xm
    @staticmethod
    def variance(alpha,xm):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        if(alpha<=2):
            return float("infinity")
        return (xm**2)*alpha/((alpha-1)**2*(alpha-2))
    @staticmethod
    def stddev(alpha,xm):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        if(alpha<=2):
            return float("infinity")
        return xm*math.sqrt(alpha/(alpha-2))/(alpha-1)
    @staticmethod
    def entropy(alpha,xm):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        return math.log((xm/alpha)*math.exp(1+1/alpha))
    @staticmethod
    def skewness(alpha,xm):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        if(alpha>3):
            return 2*(1+alpha)/(alpha-3)*math.sqrt((alpha-2)/alpha)
        return None
    @staticmethod
    def ppf(alpha,xm,q):
        if(alpha<=0 or xm<=0):
            raise InvalidInputError("alpha and xm must be positive")
        return xm/math.pow(q,1/alpha)

'''
Pareto II variant (Shifted Pareto) distribution
'''
class pareto2var(Distribution): #shifted pareto
    @staticmethod
    def random(a,b,c):
        return a*(1/(math.pow(rg0(),b))-1)+c
    @staticmethod
    def pdf(a,b,c,x):
        return b/a*math.pow(a/(x+a-c),b+1)
    @staticmethod
    def cdf(a,b,c,x):
        return 1-math.pow(a/(y+a-c),b)
    @staticmethod
    def kurtosis(a,b,c):
        return 3*a**4*b*(3*b**2+b+2)/((b-4)*(b-3)*(b-2)*(b-1)**4)
    @staticmethod
    def mean(a,b,c):
        return a/(b-1)+c
    @staticmethod
    def median(a,b,c):
        return a*(math.pow(2,1/b)-1)+c
    @staticmethod
    def mode(a,b,c):
        return a
    @staticmethod
    def variance(a,b,c):
        return a**2*b/((b-2)*(b-1)**2)
    @staticmethod
    def stddev(a,b,c):
        return math.sqrt(a**2*b/((b-2)*(b-1)**2))
    @staticmethod
    def skewness(a,b,c):
        return 2*a**3*b*(b+1)/((b-3)*(b-2)*(b-1)**3)
    @staticmethod
    def ppf(a,b,c,q):
        return a*(1/(math.pow(q,b))-1)+c

'''
Pareto III distribution*
'''
class pareto3(Distribution):
    @staticmethod
    def random(mu,sigma,gmma):
        return fellerpareto.random(mu,sigma,gmma,1,1)

'''
Pareto IV distribution*
'''
class pareto4(Distribution):
    @staticmethod
    def random(mu,sigma,gmma,a):
        return fellerpareto.random(mu,sigma,gmma,1,a)

'''
Pearson III distribution
'''
class pearson3(Distribution):
    @staticmethod
    def random(skew):
        return st.pearson3.rvs(skew)

'''
Pearson VII distribution*
'''
class pearson7(Distribution):
    @staticmethod
    def random(s,m):
        return normal.random(0,s)/math.sqrt(gamma.random(1/2,m-1/2))

'''
Planck distribution
'''
class planck(Distribution):
    @staticmethod
    def random(l):
        return math.ceil(-math.log(1-rg0())/l-1)
    @staticmethod
    def pdf(l,x):
        if(x*l>=0):
            return (1-math.exp(-l))*math.exp(-l*x)
    @staticmethod
    def cdf(l,x):
        if(x*l>=0):
            return 1-math.exp(-l*(math.floor(x)+1))
    @staticmethod
    def median(l):
        return math.ceil(-math.log(1/2)/l-1)
    @staticmethod
    def ppf(l,q):
        return math.ceil(-math.log(1-q)/l-1)

'''
Poisson distribution
'''
class poisson(Distribution):
    @staticmethod
    def pdf(lmbda,k):
        if(lmbda<=0):
            raise InvalidInputError("lambda must be positive")
        if(k%1!=0 or k<0):
            raise InvalidInputError("k must be a non-negative integer")
        return math.pow(lmbda,k)*math.exp(-lmbda)/math.factorial(k)
    @staticmethod
    def random(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lambda must be positive")
        l=math.exp(-lmbda)
        k=0
        p=1
        while(p>l):
            k+=1
            p*=r.random()
        return int(k)
    @staticmethod
    def kurtosis(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lambda must be positive")
        return 1/lmbda
    @staticmethod
    def mean(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lambda must be positive")
        return lmbda
    @staticmethod
    def median(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lambda must be positive")
        return math.floor(lmbda+1/3-0.02/lmbda)
    @staticmethod
    def mode(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lambda must be positive")
        return (math.floor(lmbda)-1,math.floor(lmbda))
    @staticmethod
    def variance(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lambda must be positive")
        return lmbda
    @staticmethod
    def stddev(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lambda must be positive")
        return math.sqrt(lmbda)
    @staticmethod
    def skewness(lmbda):
        if(lmbda<=0):
            raise InvalidInputError("lambda must be positive")
        return math.sqrt(1/lmbda)

'''
Poisson exponential distribution
'''
class poissonexp(Distribution): #poisson exponential
    @staticmethod
    def random(l,t):
        u=rg0()
        return -math.log(-math.log(1-(1-math.exp(-t))*u)/t)/l
    @staticmethod
    def pdf(l,t,x):
        return t*l*math.exp(-l*x-t*math.exp(-l*x))/(1-math.exp(-t))
    @staticmethod
    def cdf(l,t,x):
        return 1-(1-math.exp(-t*math.exp(-l*x)))/(1-math.exp(-t))
    @staticmethod
    def median(l,t):
        return -math.log(-math.log(1-(1-math.exp(-t))/2)/t)/l
    @staticmethod
    def ppf(l,t,q):
        return -math.log(-math.log(1-(1-math.exp(-t))*q)/t)/l

'''
Polya distribution*
'''
class polya(Distribution):
    @staticmethod
    def random(a,b):
        return poisson.random(gamma.random(a,b))

'''
Positive negative binomial distribution*
'''
class posnegbin(Distribution): #positive negative binomial (returns only >0)
    @staticmethod
    def random(k,p):
        n=negbin.random(k,p)
        while(n==0):
            n=negbin.random(k,p)
        return n

'''
Power distribution*
'''
class power(Distribution):
    @staticmethod
    def random(c,a,b):
        n=rg0()
        return a+b*math.pow(n,1/c)
    @staticmethod
    def median(c,a,b):
        return a+b*math.pow(1/2,1/c)
    @staticmethod
    def ppf(c,a,b,q):
        return a+b*math.pow(q,1/c)

'''
Power-log normal distribution
'''
class powerlognorm(Distribution): #power-log normal
    @staticmethod
    def random(c,s):
        return st.powerlognorm.rvs(c,s)

'''
Power-normal distribution
'''
class powernorm(Distribution): #power-normal
    @staticmethod
    def random(c):
        return st.powernorm.rvs(c)

'''
Prentice distribution*
'''
class prentice(Distribution):
    @staticmethod
    def random(z,l,a,g):
        return z+l*math.log(gamma.random(g,1)/gamma.random(a,1))

'''
Q-exponential distribution
'''
class qexp(Distribution): #q-exponential
    @staticmethod
    def random(q1,lmbda):
        return (-(1/(2-q1))*qlog(r.random(),1/(2-q1)))/lmbda
    @staticmethod
    def pdf(q1,lmbda,x):
        return (2-q1)*lmbda*qexp(-lmbda*x,q1)
    @staticmethod
    def cdf(q1,lmbda,x):
        return 1-qexp(-lmbda*x/(1/(2-q1)),1/(2-q1))
    @staticmethod
    def kurtosis(q,lmbda):
        if(q1<6/5):
            return 6*(-4*q1**3+17*q1**2-20*q1+6)/((q1-2)*(4*q1-5)*(5*q1-6))
        return None
    @staticmethod
    def mean(q1,lmbda):
        if(q1<3/2):
            return 1/(lmbda*(3-2*q1))
        return None
    @staticmethod
    def median(q1,lmbda):
        qp=1/(2-q1)
        return -qp*qlog(1/2,qp)/lmbda
    @staticmethod
    def mode(q1,lmbda):
        return 0
    @staticmethod
    def variance(q1,lmbda):
        if(q1<4/3):
            return (q1-2)/((2*q1-3)**2*(3*q1-4)*lmbda**2)
        return None
    @staticmethod
    def stddev(q1,lmbda):
        if(q1<4/3):
            return math.sqrt((q1-2)/((2*q1-3)**2*(3*q1-4)*lmbda**2))
        return None
    @staticmethod
    def skewness(q1,lmbda):
        if(q1<5/4):
            return 2/(5-4*q1)*math.sqrt((3*q1-4)/(q1-2))
        return None
    @staticmethod
    def ppf(q1,lmbda,q):
        return (-(1/(2-q1))*qlog(q,1/(2-q1)))/lmbda

'''
Q-Gaussian distribution
'''
class qgaussian(Distribution):
    @staticmethod
    def random(bb,q1,mu):
        return mu+(math.sqrt(-2*qlog((1+q1)/(3-q1),rg0()))*math.cos(2*math.pi*rg0()))/math.sqrt(bb*(3-q1))
    @staticmethod
    def median(bb,q1,mu):
        return mu+(math.sqrt(-2*qlog((1+q1)/(3-q1),1/2))*math.cos(2*math.pi/2))/math.sqrt(bb*(3-q1))
    @staticmethod
    def ppf(bb,q1,mu,q):
        return mu+(math.sqrt(-2*qlog((1+q1)/(3-q1),q))*math.cos(2*math.pi*q))/math.sqrt(bb*(3-q1))

'''
Rademacher distribution
'''
class rademacher(Distribution):
    @staticmethod
    def pdf(k):
        if(k!=-1 or k!=1):
            raise InvalidInputError("k must be -1 or 1")
        return 1/2
    @staticmethod
    def cdf(k):
        if(k!=-1 or k!=1):
            raise InvalidInputError("k must be -1 or 1")
        return 1/2+(1+k)/4
    @staticmethod
    def random():
        n=r.random()
        if(n<=0.5):
            return -1
        return 1
    @staticmethod
    def kurtosis():
        return -2
    @staticmethod
    def mean():
        return 0
    @staticmethod
    def median():
        return 0
    @staticmethod
    def mode():
        return None
    @staticmethod
    def variance():
        return 1
    @staticmethod
    def stddev():
        return 1
    @staticmethod
    def entropy():
        return math.log(2)
    @staticmethod
    def skewness():
        return 0
    @staticmethod
    def ppf(q):
        if(q<=0.5):
            return -1
        return 1

'''
Raised cosine distribution
'''
class raisedcosine(Distribution):
    @staticmethod
    def pdf(mu,s,x):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        if(math.abs(x-mu)>s):
            raise InvalidInputError("x must be between mu-s and mu+s")
        return 1/(2*s)*(1+math.cos((x-mu)/s*math.pi))
    @staticmethod
    def cdf(mu,s,x):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        if(math.abs(x-mu)>s):
            raise InvalidInputError("x must be between mu-s and mu+s")
        return 1/2*(1+(x-mu)/s+1/math.pi*math.sin((x-mu)/s*math.pi))
    @staticmethod
    def kurtosis(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return 6/5*(90-math.pi**4)/((math.pi**2-6)**2)
    @staticmethod
    def mean(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return mu
    @staticmethod
    def median(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return mu
    @staticmethod
    def mode(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return mu
    @staticmethod
    def variance(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return s**2*(1/3-2/math.pi**2)
    @staticmethod
    def stddev(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return s*math.sqrt(1/3-2/math.pi**2)
    @staticmethod
    def skewness(mu,s):
        if(s<=0):
            raise InvalidInputError("s must be positive")
        return 0

'''
Rayleigh distribution
'''
class rayleigh(Distribution):
    @staticmethod
    def pdf(sigma,x):
        if(sigma<=0 or x<0):
            raise InvalidInputError("sigma must be positive and x must be non-negative")
        return x/(sigma**2)*math.exp(-x**2/(2*sigma**2))
    @staticmethod
    def cdf(sigma,x):
        if(sigma<=0 or x<0):
            raise InvalidInputError("sigma must be positive and x must be non-negative")
        return 1-math.exp(-x**2/(2*sigma**2))
    @staticmethod
    def random(sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return sigma*math.sqrt(-2*math.log(rg0()))
    @staticmethod
    def kurtosis(sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return -(6*math.pi**2-24*math.pi+16)/((4-math.pi)**2)
    @staticmethod
    def mean(sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return sigma*math.sqrt(math.pi/2)
    @staticmethod
    def median(sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return sigma*math.sqrt(2*math.log(2))
    @staticmethod
    def mode(sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return sigma
    @staticmethod
    def variance(sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return (4-math.pi)/2*sigma**2
    @staticmethod
    def stddev(sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return math.sqrt((4-math.pi)/2)*sigma
    @staticmethod
    def entropy(sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return 1+math.log(sigma/math.sqrt(2))+np.euler_gamma/2
    @staticmethod
    def skewness(sigma):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return 2*math.sqrt(math.pi)*(math.pi-3)/math.pow(4-math.pi,3/2)
    @staticmethod
    def ppf(sigma,q):
        if(sigma<=0):
            raise InvalidInputError("sigma must be positive")
        return sigma*math.sqrt(-2*math.log(q))

'''
R distribution
'''
class rdist(Distribution):
    @staticmethod
    def random(c):
        return st.rdist.rvs(c)

'''
Reciprocal distribution
'''
class reciprocal(Distribution):
    @staticmethod
    def pdf(a,b,x):
        if(a<=0 or a>=b):
            raise InvalidInputError("a must be positive and less than b")
        if(x<a or x>b):
            raise InvalidInputError("x must be between a and b inclusive")
        return 1/(x*(math.log(b)-math.log(a)))
    @staticmethod
    def cdf(a,b,x):
        if(a<=0 or a>=b):
            raise InvalidInputError("a must be positive and less than b")
        if(x<a or x>b):
            raise InvalidInputError("x must be between a and b inclusive")
        return (math.log(x)-math.log(a))/(math.log(b)-math.log(a))
    @staticmethod
    def random(a,b):
        if(a<=0 or a>=b):
            raise InvalidInputError("a must be positive and less than b")
        return math.exp(r.random()*(math.log(b)-math.log(a))+math.log(a))
    @staticmethod
    def mean(a,b):
        if(a<=0 or a>=b):
            raise InvalidInputError("a must be positive and less than b")
        return (b-a)/(math.log(b)-math.log(a))
    @staticmethod
    def ppf(a,b,q):
        if(a<=0 or a>=b):
            raise InvalidInputError("a must be positive and less than b")
        return math.exp(q*(math.log(b)-math.log(a))+math.log(a))

'''
Reciprocal-inverse Gaussian distribution
'''
class recipinvgauss(Distribution): #reciprocal-inverse gaussian
    @staticmethod
    def random(mu):
        return st.recipinvgauss.rvs(mu)
    @staticmethod
    def pdf(mu,x):
        return 1/math.sqrt(2*math.pi*x)*math.exp(-(1-mu*x)**2/(2*x*mu**2))

'''
Rectified normal distribution*
'''
class rectifiednormal(Distribution):
    @staticmethod
    def random(mu,sigma):
        n=normal.random(mu,sigma)
        if(n<0):
            return 0
        return n

'''
Reflected power distribution*
'''
class reflectedpower(Distribution):
    @staticmethod
    def random(c):
        return 1-math.pow(1-r.random(),1/c)
    @staticmethod
    def median(c):
        return 1-math.pow(1/2,1/c)
    @staticmethod
    def random(c,q):
        return 1-math.pow(1-q,1/c)

'''
Reflected Topp-Leone distribution
'''
class reflectedtoppleone(Distribution):
    @staticmethod
    def random(aa,bb):
        n=rg0()
        if aa==1:
            return 1-math.pow(1-n,1/bb)
        if(0<aa and aa<1):
            return 1-(aa+math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1-n,1/bb)))/(2*(aa-1))
        return 1-(aa-math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1-n,1/bb)))/(2*(aa-1))
    @staticmethod
    def median(aa,bb):
        if aa==1:
            return 1-math.pow(1/2,1/bb)
        if(0<aa and aa<1):
            return 1-(aa+math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1/2,1/bb)))/(2*(aa-1))
        return 1-(aa-math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1/2,1/bb)))/(2*(aa-1))
    @staticmethod
    def ppf(aa,bb):
        if aa==1:
            return 1-math.pow(1-q,1/bb)
        if(0<aa and aa<1):
            return 1-(aa+math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1-q,1/bb)))/(2*(aa-1))
        return 1-(aa-math.sqrt(math.pow(aa,2)-4*(aa-1)*math.pow(1-q,1/bb)))/(2*(aa-1))

'''
Reversed Burr II distribution*
'''
class reversedburr2(Distribution):
    @staticmethod
    def random(g):
        return prentice.random(0,-1,1,g)

'''
Rice distribution
'''
class rice(Distribution):
    @staticmethod
    def random(nu,sigma):
        return sigma*math.sqrt(chi2.random(2*poisson.random(nu**2/(2*sigma**2))+2))

'''
Scaled chi-square distribution*
'''
class scaledchi2(Distribution): #scaled chi-squared
    @staticmethod
    def random(s,k):
        return gamma.random(2*(s**2),k/2.0)

'''
Scaled inverse chi-square distribution
'''
class scaledinvchi2(Distribution): #scaled inverse chi-squared
    @staticmethod
    def random(nu,t2):
        return invchi2(nu)*t2*nu
    @staticmethod
    def pdf(nu,t2,x):
        return math.pow(t2*nu/2,nu/2)/math.gamma(nu/2)*math.exp((-nu*t2/(2*x)))/math.pow(x,1+nu/2)
    @staticmethod
    def kurtosis(nu,t2):
        if(nu>8):
            return 12*(5*nu-22)/((nu-6)*(nu-8))
    @staticmethod
    def mean(nu,t2):
        if(nu>2):
            return nu*t2/(nu-2)
    @staticmethod
    def mode(nu,t2):
        return nu*t2/(nu+2)
    @staticmethod
    def variance(nu,t2):
        if(nu>4):
            return 2*nu**2*t2**2/((nu-2)**2*(nu-4))
    @staticmethod
    def stddev(nu,t2):
        if(nu>4):
            return 2*nu*t2/math.sqrt((nu-2)**2*(nu-4))
    @staticmethod
    def skewness(nu,t2):
        if(nu>6):
            return 4/(nu-6)*math.sqrt(2*(nu-4))

'''
Shifted Gompertz distribution
'''
class shiftedgompertz(Distribution):
    @staticmethod
    def random(b,eta):
        return max(exp.random(b),gumbel.random(b,eta))
    @staticmethod
    def pdf(b,eta,x):
        return b*math.exp(-b*x)*math.exp(-eta*math.exp(-b*x))*(1+eta*(1-math.exp(-b*x)))
    @staticmethod
    def cdf(b,eta,x):
        return (1-math.exp(-b*x))*math.exp(-eta*math.exp(-b*x))
    @staticmethod
    def mode(b,eta):
        if(0<eta and eta<=0.5):
            return 0
        z=(3+eta-math.sqrt(eta**2+2*eta+5))/(2*eta)
        return (-1/b)*math.log(z)

'''
Shifted log-logistic distribution
'''
class shiftedloglogistic(Distribution):
    @staticmethod
    def random(xi,mu,sigma):
        n=rg0()
        return (math.pow(1/n-2,-xi)*(xi*mu*math.pow(1/n-2,xi)+sigma))/xi
    @staticmethod
    def pdf(xi,mu,sigma,x):
        z=(x-mu)/sigma
        return math.pow(1+xi*z,-1/xi-1)/(sigma*(1+math.pow(1+xi*z,-1/xi))**2)
    @staticmethod
    def cdf(xi,mu,sigma,x):
        z=(x-mu)/sigma
        return 1/(1+math.pow(1+xi*z,-1/xi))
    @staticmethod
    def mean(xi,mu,sigma):
        return mu+sigma/xi*(math.pi*xi*math.csc(math.pi*xi)-1)
    @staticmethod
    def median(xi,mu,sigma):
        return mu
    @staticmethod
    def mode(xi,mu,sigma):
        return mu+sigma/xi*(math.pow((1-xi)/(1+xi),xi)-1)
    @staticmethod
    def variance(xi,mu,sigma):
        aa=math.pi*xi
        return sigma**2/xi**2*(2*aa*math.csc(2*aa)-(aa*math.csc(aa))**2)
    @staticmethod
    def stddev(xi,mu,sigma):
        return sigma/xi*math.sqrt(2*aa*math.csc(2*aa)-(aa*math.csc(aa))**2)
    @staticmethod
    def ppf(xi,mu,sigma,q):
        return (math.pow(1/n-2,-xi)*(xi*mu*math.pow(1/q-2,xi)+sigma))/xi

'''
Skellam distribution
'''
class skellam(Distribution):
    @staticmethod
    def random(mu1,mu2):
        return poisson.random(mu1)-poisson.random(mu2)
    @staticmethod
    def kurtosis(mu1,mu2):
        return 1/(mu1+mu2)
    @staticmethod
    def mean(mu1,mu2):
        return mu1-mu2
    @staticmethod
    def variance(mu1,mu2):
        return mu1+mu2
    @staticmethod
    def stddev(mu1,mu2):
        return math.sqrt(mu1+mu2)
    @staticmethod
    def skewness(mu1,mu2):
        return (mu1-mu2)/math.pow(mu1+mu2,3/2)

'''
Skew Laplace distribution*
'''
class skewlaplace(Distribution):
    @staticmethod
    def random(a,b,c):
        n=r.random()
        if(n<a):
            return a-b*math.log(b/(n*(b+c)))
        return a-b*math.log(-((n-1)*(b+c))/b)
    @staticmethod
    def median(a,b,c):
        if(1/2<a):
            return a-b*math.log(b/((b+c)/2))
        return a-b*math.log((b+c)/(2*b))
    @staticmethod
    def ppf(a,b,c,q):
        if(q<a):
            return a-b*math.log(b/(q*(b+c)))
            return a-b*math.log(-((q-1)*(b+c))/b)

'''
Skew logistic distribution*
'''
class skewlogistic(Distribution):
    @staticmethod
    def random(a,b,c):
        n=r.random()
        return math.log(math.pow(1-n,-1/c))*(a*math.pow(math.log(1-n),1/c)-b)
    @staticmethod
    def median(a,b,c):
        return math.log(math.pow(-1/2,-1/c))*(a*math.pow(math.log(1/2),1/c)-b)
    @staticmethod
    def ppf(a,b,c,q):
        return math.log(math.pow(1-q,-1/c))*(a*math.pow(math.log(1-q),1/c)-b)

'''
Skew normal distribution
'''
class skewnormal(Distribution):
    @staticmethod
    def random(a):
        return st.skewnorm.rvs(a)

'''
Slash distribution
'''
class slash(Distribution):
    @staticmethod
    def random():
        n=rg0()
        while(n==1):
            n=rg0()
        return normal.random(0,1)/n
    @staticmethod
    def pdf():
        return (st.norm.pdf(0)-st.norm.pdf(x))/x**2
    @staticmethod
    def cdf():
        if(x!=0):
            return st.norm.cdf(x)-(st.norm.pdf(0)-st.norm.pdf(x))/x
        return 1/2
    @staticmethod
    def median():
        return 0
    @staticmethod
    def mode():
        return 0

'''
Slope distribution*
'''
class slope(Distribution):
    @staticmethod
    def random(aa):
        if aa==1:
            return 1
        return (-aa+math.sqrt(math.pow(aa,2)+4*r.random()*(1-aa)))/(2*(1-aa))
    @staticmethod
    def median(aa):
        if aa==1:
            return 1
        return (-aa+math.sqrt(math.pow(aa,2)+2*(1-aa)))/(2*(1-aa))
    @staticmethod
    def ppf(aa,q):
        if aa==1:
            return 1
        return (-aa+math.sqrt(math.pow(aa,2)+4*q*(1-aa)))/(2*(1-aa))

'''
Standard beta prime distribution*
'''
class stdbetaprime(Distribution): #standard beta prime
    @staticmethod
    def random(a,b):
        c=gamma.random(b,1)
        while(c==0):
            c=gamma.random(b,1)
        return gamma.random(a,1)/c

'''
Standard Prentice distribution*
'''
class stdprentice(Distribution): #standart prentice
    @staticmethod
    def random(a,g):
        return math.log(stdbetaprime.random(a,g))

'''
Suzuki distribution
'''
class suzuki(Distribution):
    @staticmethod
    def random(mu,nu):
        return rayleigh.random(lognormal.random(0,mu,nu))
    @staticmethod
    def kurtosis(mu,nu):
        return (32*math.exp(6*nu**2)+24*math.exp(nu**2)*math.pi-24*math.exp(3*nu**2)*math.pi-3*math.pi**2)/((4*math.exp(nu**2)-math.pi)**2)
    @staticmethod
    def mean(mu,nu):
        return math.exp(mu+nu**2/2)*math.sqrt(math.pi/2)
    @staticmethod
    def variance(mu,nu):
        return math.exp(2*mu+nu**2)*(2*math.exp(nu**2)-math.pi/2)
    @staticmethod
    def stddev(mu,nu):
        return math.sqrt(math.exp(2*mu+nu**2)*(2*math.exp(nu**2)-math.pi/2))
    @staticmethod
    def skewness(mu,nu):
        return 2*math.sqrt(pi)*(-6*math.exp(nu**2)+3*math.exp(3*nu**2)+math.pi)/math.pow(4*math.exp(nu**2)-math.pi,3/2)

'''
Symmetric Prentice distribution*
'''
class symprentice(Distribution): #symmetric prentice
    @staticmethod
    def random(l,a):
        return prentice.random(0,l,a,a)

'''
T distribution
'''
class t(Distribution):
    @staticmethod
    def random(nu):
        return math.sqrt(f.random(1,nu))
    @staticmethod
    def pdf(nu,x):
        return math.gamma((nu+1)/2)/(math.sqrt(nu*math.pi)*math.gamma(nu/2))*math.pow(1+x**2/nu,-(nu+1)/2)
    @staticmethod
    def kurtosis(nu):
        if(nu>4):
            return 6/(nu-4)
        if(2<nu and nu<=4):
            return float("infinity")
        return None
    @staticmethod
    def mean(nu):
        if(nu>1):
            return 0
        return None
    @staticmethod
    def median(nu):
        return 0
    @staticmethod
    def mode(nu):
        return 0
    @staticmethod
    def variance(nu):
        if(nu>2):
            return nu/(nu-2)
        if(1<nu and nu<=2):
            return float("infinity")
        return None
    @staticmethod
    def stddev(nu):
        if(nu>2):
            return math.sqrt(nu/(nu-2))
        if(1<nu and nu<=2):
            return float("infinity")
        return None
    @staticmethod
    def skewness(nu):
        if(nu>3):
            return 0
        return None

'''
Topp-Leone Distribution
'''
class toppleone(Distribution):
    @staticmethod
    def random(bb):
        return 1-math.sqrt(1-math.pow(r.random(),1/bb))
    @staticmethod
    def pdf(bb,x):
        return bb*(2-2*x)*math.pow(2*x-x**2,bb-1)
    @staticmethod
    def cdf(bb,x):
        return math.pow(2*x-x**2,bb)
    @staticmethod
    def median(bb):
        return 1-math.sqrt(1-math.pow(1/2,1/bb))
    @staticmethod
    def ppf(bb,q):
        return 1-math.sqrt(1-math.pow(q,1/bb))

'''
Transmuted complementary Weibull-geometric distribution
'''
class transmutedcompweibullgeo(Distribution): #transmuted complementary weibull geometric
    @staticmethod
    def random(a,b,d,g):
        q=rg0()
        return math.pow(math.log((2*a*a-2*a*q*(a-1)-a*(1+d)+a*math.sqrt(1+d*(d-4*q+2)))/(2*a*a*(1-q))),1/b)/g
    @staticmethod
    def pdf(a,b,d,g,x):
        return a*b*g*math.pow(g*x,b-1)*math.exp(-math.pow(g*x,b))*(a*(1-d)-(a-a*d-d-1)*math.exp(-math.pow(g*x,b)))/math.pow(a+(1-a)*math.exp(-math.pow(g*x,b)),3)
    @staticmethod
    def median(a,b,d,g):
        return math.pow(math.log((2*a**2-a*(a-1)-a*(1+d)+a*math.sqrt(1+d**2))/(a**2)),1/b)/g
    @staticmethod
    def ppf(a,b,d,g,q):
        return math.pow(math.log((2*a*a-2*a*q*(a-1)-a*(1+d)+a*math.sqrt(1+d*(d-4*q+2)))/(2*a*a*(1-q))),1/b)/g

'''
Transmuted exponentiated exponential distribution
'''
class transmutedexponentiatedexp(Distribution): #transmuted exponentiated exponential
    @staticmethod
    def random(a,g,l):
        i=((1+l)-math.sqrt((1+l)**2-4*l*rg0()))/(2*l)
        return -math.log(1-math.pow(i,1/a))/g
    @staticmethod
    def median(a,g,l):
        i=((1+l)-math.sqrt((1+l)**2-2*l))/(2*l)
        return -math.log(1-math.pow(i,1/a))/g
    @staticmethod
    def random(a,g,l,q):
        i=((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l)
        return -math.log(1-math.pow(i,1/a))/g

'''
Transmuted exponentiated Frechet distribution
'''
class transmutedexponentiatedfrechet(Distribution):
    @staticmethod
    def random(a,b,l,t):
        return t*math.pow(-math.log(1-math.pow(((l-1)+math.sqrt((l+1)**2)-4*l*rg0())/2*l,1/a)),-1/b)
    @staticmethod
    def pdf(a,b,l,t,x):
        qx=1-math.exp(-math.pow(t/x,b))
        return a*b*math.pow(t,b)*math.pow(x,-1-b)*math.exp(-math.pow(t/x,b))*math.pow(qx,a-1)*((1-l)+2*l*math.pow(qx,a))
    @staticmethod
    def cdf(a,b,l,t,x):
        return (1-math.pow(1-math.exp(-math.pow(t/x,b)),a))*(1+l*math.pow(1-math.exp(-math.pow(t/x,b)),a))
    @staticmethod
    def median(a,b,l,t):
        return t*math.pow(-math.log(1-math.pow(((l-1)+math.sqrt(l+1**2))/2*l,1/a)),-1/b)
    @staticmethod
    def ppf(a,b,l,t,q):
        return t*math.pow(-math.log(1-math.pow(((l-1)+math.sqrt((l+1)**2)-4*l*q)/2*l,1/a)),-1/b)

'''
Transmuted exponentiated Lomax distribution
'''
class transmutedexponentiatedlomax(Distribution): #transmuted exponentiated lomax
    @staticmethod
    def random(a,g,l,t):
        u=rg0()
        i=(1+l-math.sqrt((1+l)**2-4*l*u))/(2*l)
        return (math.pow(1-math.pow(i,1/a),-1/t)-1)/g
    @staticmethod
    def pdf(a,g,l,t,x):
        return a*t*g*math.pow(1-math.pow(1+g*x,-t),a-1)/(math.pow(1+g*x,t+1))*(1+l-2*l*math.pow(1-math.pow(1+g*x,-t),a))
    @staticmethod
    def cdf(a,g,l,t,x):
        return math.pow(1-math.pow(1+g*x,-t),a)*(1+l-l*math.pow(1-math.pow(1+g*x,-t),a))
    @staticmethod
    def mean(a,g,l,t):
        return a*(1+l)/g*(sp.beta(1-1/t,a)-1/a)-2*a*l/g*(sp.beta(1-1/t,2*a)-1/(2*a))
    @staticmethod
    def median(a,g,l,t):
        i=(1+l-math.sqrt((1+l)**2-2*l))/(2*l)
        return (math.pow(1-math.pow(i,1/a),-1/t)-1)/g
    @staticmethod
    def ppf(a,g,l,t,q):
        i=(1+l-math.sqrt((1+l)**2-4*l*q))/(2*l)
        return (math.pow(1-math.pow(i,1/a),-1/t)-1)/g

'''
Transmuted generalized linear failure rate distribution
'''
class transmutedgenlinearfailurerate(Distribution): #transmuted generalized linear failure rate
    @staticmethod
    def random(a,g,l,t):
        i=((1+l)-math.sqrt((1+l)**2-4*l*rg0()))/(2*l)
        return (-t+math.sqrt(t**2-4*g*math.log(1-math.pow(i,1/a))))/(2*g)
    @staticmethod
    def median(a,g,l,t):
        i=((1+l)-math.sqrt((1+l)**2-2*l))/(2*l)
        return (-t+math.sqrt(t**2-4*g*math.log(1-math.pow(i,1/a))))/(2*g)
    @staticmethod
    def ppf(a,g,l,t,q):
        i=((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l)
        return (-t+math.sqrt(t**2-4*g*math.log(1-math.pow(i,1/a))))/(2*g)

'''
Transmuted generalized Rayleigh distribution
'''
class transmutedgenrayleigh(Distribution): #transmuted generalized rayleigh
    @staticmethod
    def random(a,b,l):
        i=((1+l)-math.sqrt((1+l)**2-4*l*rg0()))/(2*l)
        return math.sqrt(-math.log(1-math.pow(i,1/a))/b)
    @staticmethod
    def pdf(a,b,l,x):
        return 2*a*b**2*x*math.exp(-(b*x)**2)*math.pow(1-math.exp(-(b*x)**2),a-1)*(1+l-2*l*math.pow(1-math.exp(-(b*x)**2),a))
    @staticmethod
    def cdf(a,b,l,x):
        return math.pow(1-math.exp(-(b*x)**2),a)*(1+l-2*l*math.pow(1-math.exp(-(b*x)**2),a))
    @staticmethod
    def median(a,b,l):
        i=((1+l)-math.sqrt(1+l**2)/(2*l))
        return math.sqrt(-math.log(1-math.pow(i,1/a))/b)
    @staticmethod
    def ppf(a,b,l,q):
        i=((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l)
        return math.sqrt(-math.log(1-math.pow(i,1/a))/b)

'''
Transmuted inverse exponential distribution
'''
class transmutedinvexp(Distribution): #transmuted inverse exponential
    @staticmethod
    def random(a,l):
        u=rg0()
        i=2*l/((1+l)-math.sqrt((1+l)**2-4*l*u))
        return a/math.log(i)

'''
Transmuted Kumaraswamy distribution
'''
class transmutedkumaraswamy(Distribution):
    @staticmethod
    def random(a,l,t):
        q=rg0()
        return math.pow(1-math.pow(1-((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l),1/t),1/a)
    @staticmethod
    def pdf(a,l,t,x):
        return a*t*math.pow(x,a-1)*math.pow(1-math.pow(x,a),t-1)*(1-l+2*l*math.pow(1-math.pow(x,a),t))
    @staticmethod
    def cdf(a,l,t,x):
        return (1-math.pow(1-math.pow(x,a),t))*(1+l*math.pow(1-math.pow(x,a),t))
    @staticmethod
    def median(a,l,t):
        return math.pow(1-math.pow(1-((1+l)-math.sqrt(1+l**2))/(2*l),1/t),1/a)
    @staticmethod
    def ppf(a,l,t,q):
        return math.pow(1-math.pow(1-((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l),1/t),1/a)

'''
Transmuted modified inverse Rayleigh distribution
'''
class transmutedmodifiedinvrayleigh(Distribution): #transmuted modified inverse rayleigh
    @staticmethod
    def random(a,b,l):
        return 2*b/(-a+math.sqrt(a**2-4*b*math.log(((1+l)-math.sqrt((1+l)**2)-4*l*rg0())/(2*l))))
    @staticmethod
    def pdf(a,b,l,x):
        return (a+2*b/x)*(1/x**2)*math.exp(-a/x-b/x**2)*(1+l-2*l*math.exp(-a/x-b/x**2))
    @staticmethod
    def cdf(a,b,l,x):
        return math.exp(-a/x-b/x**2)*(1+l-l*math.exp(-a/x-b/x**2))
    @staticmethod
    def median(a,b,l):
        return 2*b/(-a+math.sqrt(a**2-4*b*math.log(((1+l)-math.sqrt(1+l**2))/(2*l))))
    @staticmethod
    def ppf(a,b,l,q):
        return 2*b/(-a+math.sqrt(a**2-4*b*math.log(((1+l)-math.sqrt((1+l)**2)-4*l*q)/(2*l))))

'''
Transmuted Weibull-Lomax distribution
'''
class transmutedweibulllomax(Distribution):
    @staticmethod
    def random(a,b,aa,bb,l):
        u=rg0()
        D=0
        if(l==0):
            D=u
        else:
            D=((1+l)-math.sqrt((1+l)**2-4*l*u))/(2*l)
        return bb*(math.pow(math.pow(math.log(math.pow(1-D,-1/a)),1/b)+1,1/aa)-1)
    @staticmethod
    def pdf(a,b,aa,bb,l,x):
        return a*b*aa/bb*math.pow(1+x/bb,b*a-1)*math.pow(1-math.pow(1+x/bb,-aa),b-1)*math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b))*((1+l)-2*l*(1-math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b))))
    @staticmethod
    def cdf(a,b,aa,bb,l,x):
        return (1-math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b)))*((1+l)-l*(1-math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b))))
    @staticmethod
    def median(a,b,aa,bb,l):
        D=0
        if(l==0):
           D=1/2
        else:
           D=((1+l)-math.sqrt(1+l**2))/(2*l)
        return bb*(math.pow(math.pow(math.log(math.pow(1-D,-1/a)),1/b)+1,1/aa)-1)
    @staticmethod
    def ppf(a,b,aa,bb,l,q):
        D=0
        if(l==0):
            D=q
        else:
            D=((1+l)-math.sqrt((1+l)**2-4*l*q))/(2*l)
        return bb*(math.pow(math.pow(math.log(math.pow(1-D,-1/a)),1/b)+1,1/aa)-1)

'''
Triangular distribution
'''
class triangular(Distribution):
    @staticmethod
    def pdf(a,b,c,x):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        if(x>a or b>x):
            raise InvalidInputError("x must be between a and b, inclusive")
        if(x<c):
            return 2*(x-a)/((b-a)*(c-a))
        if(x==c):
            return 2/(b-a)
        if(c<x):
            return 2*(b-x)/((b-a)*(b-c))
    @staticmethod
    def cdf(a,b,c,x):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        if(x>a or b>x):
            raise InvalidInputError("x must be between a and b, inclusive")
        if(x<=c):
            return (x-a)**2/((b-a)*(c-a))
        if(c<x and x<b):
            return 1-(b-x)**2/((b-a)*(b-c))
        if(x==b):
            return 1
    @staticmethod
    def random(a,b,c):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        n=r.random()
        f=(c-a)/(b-a)
        if(n<=f):
            return a+math.sqrt(n*(b-a)*(c-a))
        return b-math.sqrt((1-n)*(b-a)*(c-a))
    @staticmethod
    def kurtosis(a,b,c):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        return -3/5
    @staticmethod
    def mean(a,b,c):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        return (a+b+c)/3
    @staticmethod
    def median(a,b,c):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        if(c>=(a+b)/2):
            return a+math.sqrt((b-a)*(c-a)/2)
        return b-math.sqrt((b-a)*(b-c)/2)
    @staticmethod
    def mode(a,b,c):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        return c
    @staticmethod
    def variance(a,b,c):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        return (a**2+b**2+c**2-a*b-a*c-b*c)/18
    @staticmethod
    def stddev(a,b,c):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        return math.sqrt((a**2+b**2+c**2-a*b-a*c-b*c)/18)
    @staticmethod
    def entropy(a,b,c):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        return 1/2+math.log((b-a)/2)
    @staticmethod
    def skewness(a,b,c):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        return math.sqrt(2)*(a+b-2*c)*(2*a-b-c)*(a-2*b+c)/(5*math.pow(a**2+b**2+c**2-a*b-a*c-b*c,3/2))
    @staticmethod
    def ppf(a,b,c,q):
        if(b>=a or a>c or c>b):
            raise InvalidInputError("b must be greater than c, which must be greater than a")
        f=(c-a)/(b-a)
        if(q<=f):
            return a+math.sqrt(q*(b-a)*(c-a))
        return b-math.sqrt((1-q)*(b-a)*(c-a))

'''
Truncated exponential distribution
'''
class truncexp(Distribution): #truncated exponential
    @staticmethod
    def random(b):
        return st.truncexp.rvs(b)

'''
Truncated normal distribution (truncated on both ends)*
'''
class truncnormb(Distribution): #truncated normal (truncated on both ends)
    @staticmethod
    def random(a,b,mu,sigma):
        n=normal.random(mu,sigma**2)
        while (n<a or n>b):
            n=normal.random(mu,sigma**2)
        return n

'''
Truncated normal distribution (truncated on the left)*
'''
class truncnorml(Distribution): #left-truncated normal #truncated on the left
    @staticmethod
    def random(a,mu,sigma):
        n=normal.random(mu,sigma**2)
        while n<a:
            n=normal.random(mu,sigma**2)
        return n

'''
Truncated normal distribution (truncated on the right)*
'''
class truncnormr(Distribution): #right-truncated normal #truncated on the right
    @staticmethod
    def random(b,mu,sigma):
        n=normal.random(mu,sigma**2)
        while n>b:
            n=normal.random(mu,sigma**2)
        return n

'''
Tsallis-q distribution*
'''
class tsallisq(Distribution):
    @staticmethod
    def random(lmbda,q1):
        x=rg0()
        if(q1==1):
            return -math.log(1-x)/lmbda
        return -(1-math.pow(1-x,1/(-2+q1))+math.pow(1-x,1/(-2+q1))*x)/((-1+q1)*lmbda)
    @staticmethod
    def median(lmbda,q1):
        if(q1==1):
            return -math.log(1/2)/lmbda
        return -(1-math.pow(1/2,1/(-2+q1))+math.pow(1/2,1/(-2+q1))/2)/((-1+q1)*lmbda)
    @staticmethod
    def random(lmbda,q1):
        if(q1==1):
            return -math.log(1-q)/lmbda
        return -(1-math.pow(1-q,1/(-2+q1))+math.pow(1-q,1/(-2+q1))*q)/((-1+q1)*lmbda)

'''
Tukey-lambda distribution*
'''
class tukeylambda(Distribution):
    @staticmethod
    def random(lmbda):
        p=r.random()
        if(lmbda==0):
            while(p==1):
                p=r.random()
            return math.log(p/(1-p))
        return 1/lmbda*(math.pow(p,lmbda)-math.pow(1-p,lmbda))
    @staticmethod
    def median(lmbda):
        return 0
    @staticmethod
    def ppf(lmbda,q):
        if(lmbda==0):
            return math.log(q/(1-q))
        return 1/lmbda*(math.pow(q,lmbda)-math.pow(1-q,lmbda))

'''
Two-sided power distribution*
'''
class twosidedpower(Distribution):
    @staticmethod
    def random(d,theta):
        n=rg0()
        if(n<=theta):
            return theta*math.pow(n/theta,1/d)
        return 1-(1-theta)*math.pow((1-n)/(1-d),1/d)
    @staticmethod
    def median(d,theta):
        if(1/2<=theta):
            return theta*math.pow(1/(2*theta),1/d)
        return 1-(1-theta)*math.pow(1/(2*(1-d)),1/d)
    @staticmethod
    def ppf(d,theta,q):
        if(q<=theta):
            return theta*math.pow(q/theta,1/d)
            return 1-(1-theta)*math.pow((1-q)/(1-d),1/d)

'''
Uniform distribution
'''
class uniform(Distribution):
    @staticmethod
    def pdf(a,b,x):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        if(x<a or x>b):
            raise InvalidInputError("x must be between a and b inclusive")
        return 1/(b-a)
    @staticmethod
    def cdf(a,b,x):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        if(x<a or x>b):
            raise InvalidInputError("x must be between a and b inclusive")
        return (x-a)/(b-a)
    @staticmethod
    def random(a,b):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        return ((b-a)*r.random()+a)
    @staticmethod
    def kurtosis(a,b):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        return -6/5
    @staticmethod
    def mean(a,b):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        return 1/2*(a+b)
    @staticmethod
    def median(a,b):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        return 1/2*(a+b)
    @staticmethod
    def variance(a,b):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        return 1/12*((b-a)**2)
    @staticmethod
    def stddev(a,b):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        return 1/math.sqrt(12)*(b-a)
    @staticmethod
    def entropy(a,b):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        return math.log(b-a)
    @staticmethod
    def skewness(a,b):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        return 0
    @staticmethod
    def ppf(a,b,q):
        if(b<=a):
            raise InvalidInputError("b must be greater than a")
        return ((b-a)*q+a)

'''
Uniform product distribution*
'''
class uniformproduct(Distribution):
    @staticmethod
    def random(n):
        prod=1
        for _ in range(n):
            prod*=r.random()
        return prod

'''
Unit gamma distribution*
'''
class unitgamma(Distribution):
    @staticmethod
    def random(a,b):
        return math.exp(gamma.random(-1/b,a))

'''
U power distribution
'''
class upower(Distribution):
    @staticmethod
    def random(k):
        return math.pow((2*r.random()-1),1/(2*k+1))
    @staticmethod
    def pdf(k,x):
        return (2*k+1)/2*math.pow(x,2*k)
    @staticmethod
    def cdf(k,x):
        return 1/2*(1+math.pow(x,2*k+1))
    @staticmethod
    def kurtosis(k):
        return (2*k+3)**2/((2*k+5)*(2*k+1))
    @staticmethod
    def mean(k):
        return 0
    @staticmethod
    def median(k):
        return 0
    @staticmethod
    def variance(k):
        return (2*k+1)/(2*k+3)
    @staticmethod
    def stddev(k):
        return math.sqrt((2*k+1)/(2*k+3))
    @staticmethod
    def skewness(k):
        return 0
    @staticmethod
    def ppf(k,q):
        return math.pow((2*q-1),1/(2*k+1))

'''
Upsilon distribution
'''
class upsilon(Distribution):
    @staticmethod
    def random(ts, nus):
        y=0
        for i in range(len(ts)):
            y+=ts[i]*math.sqrt(chi2.random(nus[i]))
        return y

'''
U quadratic distribution*
'''
class uquad(Distribution):
    @staticmethod
    def random(a,b):
        c=12/(math.pow((b-a),3))
        d=(a+b)/2
        return 3*r.random()/c-math.pow(math.pow((d-c),3),1/3.0)+d
    @staticmethod
    def median(a,b):
        c=12/(math.pow((b-a),3))
        d=(a+b)/2
        return 3*(1/2)/c-math.pow(math.pow((d-c),3),1/3.0)+d
    @staticmethod
    def ppf(a,b,q):
        c=12/(math.pow((b-a),3))
        d=(a+b)/2
        return 3*q/c-math.pow(math.pow((d-c),3),1/3.0)+d

'''
Voigt distribution*
'''
class voigt(Distribution):
    @staticmethod
    def random(a,s,sigma):
        return normal.random(0,sigma)+cauchy.random(a,s)

'''
Wakeby distribution
'''
class wakeby(Distribution):
    @staticmethod
    def random(a,b,c,d,m):
        n=r.random()
        #return -a*pow(1-n,b)+c*pow(1-n,-d)+e
        return m+a*(1-math.pow(1-n,b))-c*(1-math.pow(1-n,-d))
    @staticmethod
    def median(a,b,c,d,m):
        return m+a*(1-math.pow(1/2,b))-c*(1-math.pow(1/2,-d))
    @staticmethod
    def random(a,b,c,d,m,q):
        return m+a*(1-math.pow(1-q,b))-c*(1-math.pow(1-q,-d))

'''
Wedge distribution*
'''
class wedge(Distribution):
    @staticmethod
    def random(a,s):
        return power.random(a,s,2)

'''
Weibull distribution
'''
class weibull(Distribution):
    @staticmethod
    def pdf(k,lmbda,x):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        if(x<0):
            raise InvalidInputError("x must be non-negative")
        return k/lmbda*math.pow(x/lmbda,k-1)*math.exp(-math.pow(x/lmbda,k))
    @staticmethod
    def cdf(k,lmbda,x):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        if(x<0):
            raise InvalidInputError("x must be non-negative")
        return 1-math.exp(-math.pow(x/lmbda,k))
    @staticmethod
    def random(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        return lmbda*math.pow((-math.log(rg0())),1/k)
    @staticmethod
    def kurtosis(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        g1=math.gamma(1+1/k)
        g2=math.gamma(1+2/k)
        g3=math.gamma(1+3/k)
        g4=math.gamma(1+4/k)
        return (-6*g1**4+12*g1**2*g2-3*g2**2-4*g1*g3+g4)/((g2-g1**2)**2)
    @staticmethod
    def mean(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        return lmbda*math.gamma(1+1/k)
    @staticmethod
    def median(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        return lmbda*math.pow(math.log(2),1/k)
    @staticmethod
    def mode(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        if(k>1):
            return lmbda*math.pow((k-1)/k,1/k)
        return 0
    @staticmethod
    def variance(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        return lmbda**2*(math.gamma(1+2/k)-(math.gamma(1+1/k)**2))
    @staticmethod
    def stddev(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        return lmbda*math.sqrt(math.gamma(1+2/k)-(math.gamma(1+1/k)**2))
    @staticmethod
    def entropy(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        return np.euler_gamma*(1-1/k)+math.log(lmbda/k)+1
    @staticmethod
    def skewness(k,lmbda):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        return (math.gamma(1+3/k)*lmbda**3-3*weibull.mean(k,lmbda)*weibull.variance(k,lmbda)-weibull.mean(k,lmbda)**3)/math.pow(weibull.variance,3/2)
    @staticmethod
    def ppf(k,lmbda,q):
        if(lmbda<=0 or k<=0):
            raise InvalidInputError("k and lambda must be positive")
        return lmbda*math.pow((-math.log(q)),1/k)

'''
Weibull-Burr XII distribution*
'''
class weibullburr12(Distribution):
    @staticmethod
    def random(a,b,c,k,s):
        return s*math.pow(math.pow(math.pow((-math.log(1-rg0())/a),1/b)+1,1/k)-1,1/c)
    @staticmethod
    def median(a,b,c,k,s):
        return s*math.pow(math.pow(math.pow((-math.log(1/2)/a),1/b)+1,1/k)-1,1/c)
    @staticmethod
    def ppf(a,b,c,k,s,q):
        return s*math.pow(math.pow(math.pow((-math.log(1-q)/a),1/b)+1,1/k)-1,1/c)

'''
Weibull-Frechet distribution
'''
class weibullfrechet(Distribution):
    @staticmethod
    def random(a,b,aa,bb):
        u=rg0()
        return aa*math.pow(math.log(1+math.pow(-math.log(1-u)/a,-1/b)),-1/bb)
    @staticmethod
    def pdf(a,b,aa,bb,x):
        return a*b*bb*math.pow(aa,bb)*math.pow(x,-bb-1)*math.exp(-b*math.pow(aa/x,bb))*math.pow(1-math.exp(-math.pow(aa/x,bb))-b-1)*math.exp(-a*math.pow(math.exp(math.pow(aa/x,bb)-1),-b))
    @staticmethod
    def cdf(a,b,aa,bb,x):
        return 1-math.exp(-a*math.pow(math.exp(math.pow(aa/x,bb)-1),-b))
    @staticmethod
    def median(a,b,aa,bb):
        return aa*math.pow(math.log(1+math.pow(-math.log(1/2)/a,-1/b)),-1/bb)
    @staticmethod
    def ppf(a,b,aa,bb,q):
        return aa*math.pow(math.log(1+math.pow(-math.log(1-q)/a,-1/b)),-1/bb)

'''
Weibull-generalized exponential distribution
'''
class weibullgenexp(Distribution): #weibull-generalized exponential
    @staticmethod
    def random(a,b,l):
        return 1/l*math.log(1+math.pow(1/a*math.log(1-rg0()),1/b))
    @staticmethod
    def pdf(a,b,l,x):
        return a*b*l*math.exp(l*x)*math.pow(math.exp(l*x)-1,b-1)*math.exp(-a*math.pow(math.exp(l*x)-1,b))
    @staticmethod
    def cdf(a,b,l,x):
        return 1-math.exp(-a*math.pow(math.exp(l*x)-1,b))
    @staticmethod
    def median(a,b,l):
        return 1/l*math.log(1+math.pow(1/a*math.log(1/2),1/b))
    @staticmethod
    def ppf(a,b,l,q):
        return 1/l*math.log(1+math.pow(1/a*math.log(1-q),1/b))

'''
Weibull-Lomax distribution
'''
class weibulllomax(Distribution):
    @staticmethod
    def random(a,b,aa,bb):
        return bb*(math.pow(math.pow(-math.log(1-rg0())/a,1/b)+1,1/aa)-1)
    @staticmethod
    def pdf(a,b,aa,bb,x):
           return a*b*aa/bb*math.pow(1+x/bb,b*aa-1)*math.pow(1-math.pow(1+x/bb,-aa),b-1)*math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b))
    @staticmethod
    def cdf(a,b,aa,bb,x):
        return 1-math.exp(-a*math.pow(math.pow(1+x/bb,aa)-1,b))
    @staticmethod
    def median(a,b,aa,bb):
        return bb*(math.pow(math.pow(-math.log(1/2)/a,1/b)+1,1/aa)-1)
    @staticmethod
    def ppf(a,b,aa,bb,q):
        return bb*(math.pow(math.pow(-math.log(1-q)/a,1/b)+1,1/aa)-1)

'''
Weibull-uniform distribution
'''
class weibulluniform(Distribution):
    @staticmethod
    def random(a,b,t):
        return t/(1+math.pow(-math.log(1-rg0())/a,-1/b))
    @staticmethod
    def pdf(a,b,t,x):
        return t*a*b/((t-x)**2)*math.pow(x/(t-x),b-1)*math.exp(-a*math.pow(x/(t-x),b))
    @staticmethod
    def cdf(a,b,t,x):
        return 1-math.exp(-a*math.pow(x/(t-x),b))
    @staticmethod
    def median(a,b,t):
        return t/(1+math.pow(-math.log(1/2)/a,-1/b))
    @staticmethod
    def ppf(a,b,t,q):
        return t/(1+math.pow(-math.log(1-q)/a,-1/b))

'''
Weibull-Weibull distribution*
'''
class weibullweibull(Distribution):
    @staticmethod
    def random(a,b,g,l):
        return math.pow(math.log(math.pow(-math.log(1-u)/a,1/b)+1)/l,1/g)
    @staticmethod
    def median(a,b,g,l):
        return math.pow(math.log(math.pow(-math.log(1/2)/a,1/b)+1)/l,1/g)
    @staticmethod
    def ppf(a,b,g,l,q):
        return math.pow(math.log(math.pow(-math.log(1-q)/a,1/b)+1)/l,1/g)

'''
Wein distribution*
'''
class wein(Distribution):
    @staticmethod
    def random(t):
        return gamma.random(t,4)

'''
Xgamma distribution
'''
class xgamma(Distribution):
    @staticmethod
    def random(t):
        u=rg0()
        v=exp.random(t)
        w=gamma.random(3,t)
        if(u<=t/(1+t)):
            return v
        return w
    @staticmethod
    def pdf(t,x):
        return t**2/(1+t)*(1+t/2*x**2)*math.exp(-t*x)
    @staticmethod
    def cdf(t,x):
        return 1-(1+t+t*x+(t*x)**2/2)/(1+t)*math.exp(-t*x)
    @staticmethod
    def mode(t):
        if(0<t and t<1/2):
           return (1+math.sqrt(1-2*t))/t
        return 0

'''
Yule-Simon distribution
'''
class yulesimon(Distribution):
    @staticmethod
    def random(rho):
        a=r.random()
        i=0
        while True:
            b=1-i*sp.beta(i,a+1)
            if(b>a):
                return i
    @staticmethod
    def pdf(rho,x):
        return rho*sp.beta(x,rho+1)
    @staticmethod
    def cdf(rho,x):
        return 1-x*sp.beta(x,rho+1)
    @staticmethod
    def kurtosis(rho):
        if(rho>4):
           return rho+3+(11*rho**3-49*rho-22)/((rho-4)*(rho-3)*rho)
    @staticmethod
    def mean(rho):
        if(rho>1):
            return rho/(rho-1)
    @staticmethod
    def mode(rho):
        return 1
    @staticmethod
    def variance(rho):
        if(rho>2):
            return rho**2/((rho-1)**2*(rho-2))
    @staticmethod
    def stddev(rho):
        if(rho>2):
            return rho/math.sqrt((rho-1)**2*(rho-2))
    @staticmethod
    def skewness(rho):
        if(rho>3):
            return (rho+1)**2*math.sqrt(rho-2)/((rho-3)*rho)

'''
Zero-inflated negative binomial-generalized exponential distribution*
'''
class zeroinflatednegbingenexp(Distribution): #zero-inflated negative binomial-generalized exponential
    @staticmethod
    def random(a,b,f,r):
        u=rg0()
        l=-1/b*math.log(1-math.pow(u,1/a))
        y=rnegbin(r,math.exp(-l))
        uu=rg0()
        if(uu>f):
            return y
        return x

'''
Zero-truncated Poisson distribution*
'''
class zerotruncpoisson(Distribution): #zero-truncated poisson
    @staticmethod
    def random(lmbda):
        k=1
        t=math.exp(-lmbda)/((1-math.exp(-lmbda))*lmbda)
        s=t
        u=r.random()
        while(s<u):
            k+=1
            t=t*lmbda/k
            s=s+t
        return k

'''
Zipf distribution
'''
class zipf(Distribution):
    @staticmethod
    def random(a):
        return st.zipf.rvs(a)
