### **ACM数论专题**
#### **1、常用数论算法**
```C++
void get_phi()  
{  
    for (int i=2;i<maxn;i++)  
        phi[i]=0;  
    phi[1]=1;  
    for (int i=2;i<maxn;i++)   
        if (!phi[i])  
            for (int j=i;j<maxn;j+=i)  
            {  
                if (!phi[j]) phi[j]=j;  
                phi[j]=phi[j]/i*(i-1);  
            }  
}  
void get_mu()  
{  
    int cnt=0;  
    memset(vis,false,sizeof(vis));  
    mu[1]=1;  
    for (int i=2;i<maxn;i++)  
    {  
        if (!vis[i])  
        {  
            prim[++cnt]=i;  
            mu[i]=-1;  
        }  
        for (int j=1;j<=cnt && prim[j]*i<maxn;j++)  
        {  
            vis[prim[j]*i]=1;  
            if (i%prim[j]==0) break;  
            else mu[i*prim[j]]=-mu[i];  
        }  
    }  
} 
void get_inv(int n,int mod)  
{  
    inv[1]=1;  
    for(int i=2;i<=n;i++)  
        inv[i]=(ll)(mod-mod/i)*inv[mod%i]%mod;  
}  

void get_fac()  
{  
    fac[1]=1;  
    for (int i=2;i<maxn;i++)   
        fac[i]=fac[i-1]*i%mod;  
    facinv[maxn-1]=inv(fac[maxn-1]);  
    for (int i=maxn-2;i>=0;i--)  
        facinv[i]=facinv[i+1]*(i+1)%mod;  
} 
void get_phi()  
{  
    int i,j,k=0;  
    for (i=2;i<maxn;i++)  
    {  
        if (!unprime[i])  
        {  
            prime[k++]=i;  
            phi[i]=i-1;  
        }  
        for (j=0;j<k && prime[j]*i<maxn;j++)  
        {  
            unprime[prime[j]*i]=true;  
             if (i%prime[j])
                phi[prime[j]*i]=phi[i]*(prime[j]-1);  
            else  
            {  
                phi[prime[j]*i]=phi[i]*prime[j];  
                break;  
            }  
        }  
    }  
} 
```
#### **2、中国剩余定理**
```
//求整数x和y，使得ax+by=d，且abs(x)+abs(y)最小  
//其中d=gcd(a,b)   
void gcd(ll a,ll b,ll &d,ll &x,ll &y)   
{  
    if (!b)   
    {  
        d=a,x=1,y=0;  
    }  
    else  
    {  
        gcd(b,a%b,d,y,x);  
        y-=x*(a/b);  
    }  
}  
//n个方程：x=a[i] (mod m[i])  
ll china(int n,int *a,int *m)  
{  
    ll M=1,d,y,x=0;  
    for (int i=0;i<n;i++)  
        M*=m[i];  
    for (int i=0;i<n;i++)  
    {  
        ll w=M/m[i];  
        gcd(m[i],w,d,d,y);  
        x=(x+y*w*a[i])%M;  
    }  
    return (x+M)%M;  
}  
```
**拓展中国剩余定理（可以处理m[i]不互质情况）**
```C++
ll gcd(ll a,ll b,ll &x,ll &y)
{
    if (!b)
    {
        x=1,y=0;
        return a;
    }
    else
    {
        ll ret=gcd(b,a%b,y,x);
        y-=x*(a/b);
        return ret;
    }
}
ll china(ll *m,ll *r,ll n)
{
    if (!n) return 0;
    ll M=m[0],R=r[0],x,y,d;
    for (int i=0;i<n;i++)
    {
        d=gcd(M,m[i],x,y);
        if ((r[i]-R)%d) return -1;
        x=(r[i]-R)/d*x%(m[i]/d);
        R+=x*M;
        M=M/d*m[i];
        R%=M;
    }
    return R>=0?R:R+M;
}
```
#### **4、大步小步算法**
```C++
int log_mod(int a,int b)  
{  
    int m,v,e=1,i;  
    m=(int)sqrt(mod);  
    v=inv(pow_mod(a,m));  
    map <int,int> x;  
    x[1]=0;  
    for (int i=1;i<m;i++)  
    {  
        e=mul_mod(e,a);  
        if (!x.count(e)) x[e]=i;  
    }  
    for (int i=0;i<m;i++)  
    {  
        if (x.count(b)) return i*m+x[b];  
        b=mul_mod(b,v);  
    }  
    return -1;  
}  
```
#### **5、莫比乌斯反演**
（1）莫比乌斯反演公式
**约数的莫比乌斯反演：**
若：$f(n)=\sum_{d|n}^{}g(d) $
则：$g(n)=\sum_{d|n}\mu(d)f(\frac{n}{d})$
**倍数的莫比乌斯反演：**
若：$f(n)=\sum_{n|d}^{}g(d) $
则：$g(n)=\sum_{n|d}\mu(\frac{d}{n})f(d)$
（2）莫比乌斯函数
$ \mu(x)=\left\{
\begin{aligned}
1 & & x=1 \\
0 & & x存在平方因子 \\
-1 & & x有奇数个质因子 \\
1 & & x偶奇数个质因子
\end{aligned}
\right.
$
#### **6、FFT和NTT算法**
（1）FFT
```C++
#include <bits/stdc++.h>
using namespace std;
typedef complex <double> cd;
typedef long long ll;
const int maxl=500000+5;
const double PI=acos(-1);
const int mod=998244353;
ll a[maxl],b[maxl],c[maxl];
cd p[maxl],q[maxl];
int rev[maxl];
void getrev(int bit)
{
    for (int i=0;i<(1<<bit);i++)
	{
        rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    }
}
void fft(cd* a,int n,int dft)
{
    for (int i=0;i<n;i++)
	{
        if (i<rev[i])
            swap(a[i],a[rev[i]]);
    }
    for (int step=1;step<n;step<<=1)
	{
        cd wn=exp(cd(0,dft*PI/step));
        for (int j=0;j<n;j+=step<<1)
		{
            cd wnk(1,0);
            for (int k=j;k<j+step;k++)
			{
                cd x=a[k];
                cd y=wnk*a[k+step];
                a[k]=x+y;
                a[k+step]=x-y;
                wnk*=wn;
            }
        }
    }
    if (dft==-1)
	{
        for (int i=0;i<n;i++)
            a[i]/=n;
    }
}
/*
	求多项式a和多项式b的乘积 
*/
void solvefft(int n,int m)
{
	int bit=1,s=2;
	while ((1<<bit)<n+m+1) bit++,s<<=1;
	getrev(bit);
	memset(p,0,sizeof(p));
	memset(q,0,sizeof(q));
	for (int i=0;i<=n;i++)
		p[i]=a[i];
	for (int i=0;i<=m;i++)
		q[i]=b[i];
	fft(p,s,1);
	fft(q,s,1);
	for (int i=0;i<s;i++) p[i]*=q[i];
	fft(p,s,-1);
	for (int i=0;i<=n+m;i++) 
		c[i]=(ll)(p[i].real()+0.5)%mod;	
}
int main()
{
	int n,m;
	scanf("%d%d",&n,&m);
	for (int i=0;i<=n;i++)
		scanf("%lld",&a[i]);
	for (int i=0;i<=m;i++)
		scanf("%lld",&b[i]);
	solvefft(n,m);
	for (int i=0;i<=n+m;i++)
		printf("%lld ",c[i]);
	printf("\n");
	return 0;
}
```
（2）NTT
```C++
#include <bits/stdc++.h>
using namespace std;
typedef complex<double> cd;
typedef long long ll;
const int maxl=300000+5;
const double PI=acos(-1);
const int mod=998244353;
int a[maxl],b[maxl],c[maxl];
ll A[maxl],B[maxl];
const int g=3;  //原根
ll quick_mod(ll a,ll b)
{
    ll ans=1;
    for(;b;b/=2)
    {
        if(b&1)
            ans=ans*a%mod;
        a=a*a%mod;
    }
    return ans;
}
int rev(int x,int r)  //蝴蝶操作
{
    int ans=0;
    for(int i=0; i<r; i++)
    {
        if(x&(1<<i))
        {
            ans+=1<<(r-i-1);
        }
    }
    return ans;
}
void NTT(int n,ll A[],int on) // 长度为N (2的次数) 
{
    int r=0;
    for(;; r++)
    {
        if((1<<r)==n)
            break;
    }
    for(int i=0; i<n; i++)
    {
        int tmp=rev(i,r);
        if(i<tmp)
            swap(A[i],A[tmp]);
    }
    for(int s=1; s<=r; s++)
    {
        int m=1<<s;
        ll wn=quick_mod(g,(mod-1)/m);
        for(int k=0; k<n; k+=m)
        {
            ll  w=1;
            for(int j=0; j<m/2; j++)
            {
                ll t,u;
                t=w*(A[k+j+m/2]%mod)%mod;
                u=A[k+j]%mod;
                A[k+j]=(u+t)%mod;
                A[k+j+m/2]=((u-t)%mod+mod)%mod;
                w=w*wn%mod;
            }
        }
    }
    if(on==-1)
    {
        for(int i=1;i<n/2;i++)
            swap(A[i],A[n-i]);
        ll inv=quick_mod(n,mod-2);
        for(int i=0;i<n;i++)
            A[i]=A[i]%mod*inv%mod;
    }
}
/*
	求多项式a和多项式b的乘积 
*/
void solventt(int n,int m)
{
	int bit=1,s=2;
	while ((1<<bit)<n+m+1) bit++,s<<=1;
	memset(A,0,sizeof(A));
	memset(B,0,sizeof(B));
	for (int i=0;i<=n;i++)
		A[i]=a[i];
	for (int i=0;i<=m;i++)
		B[i]=b[i];
	NTT(s,A,1);
    NTT(s,B,1);
    for(int i=0;i<s;i++)
        A[i]=A[i]*B[i]%mod;
    NTT(s,A,-1);
	for (int i=0;i<=n+m;i++) 
		c[i]=(c[i]+A[i])%mod;	
}
int main()
{
	int n,m;
	scanf("%d%d",&n,&m);
	for (int i=0;i<=n;i++)
		scanf("%d",&a[i]);
	for (int i=0;i<=m;i++)
		scanf("%d",&b[i]);
	solventt(n,m);
	for (int i=0;i<=n+m;i++)
		printf("%d ",c[i]);
	printf("\n");
	return 0;
}
```
#### **7、杜教筛**
（1）模板1：求$ \sum_{i=1}^{n} \varphi(i),\mu(i) $。
```C++
#include<bits/stdc++.h>
#define maxx 2500005
using namespace std;
int prime[maxx],cnt;
bool p[maxx];
long long mu[maxx],phi[maxx];
map<long long ,long long>M,P;//用map来记忆化
void init()
{
    p[1]=mu[1]=phi[1]=1;
    for(int i=2;i<maxx;i++)
    {
        if(!p[i])
        {
            prime[cnt++]=i;
            mu[i]=-1;
            phi[i]=i-1;
        }
        for(int j=0;j<cnt&&i*prime[j]<N;j++)
        {
            p[i*prime[j]]=true;
            if(i%prime[j])
            {
                mu[i*prime[j]]=-mu[i];
                phi[i*prime[j]]=phi[i]*phi[prime[j]];
            }
            else
            {
                mu[i*prime[j]]=0;
                phi[i*prime[j]]=phi[i]*prime[j];
                break;
            }
        }
    }
    /*for(int i=1;i<100;i++)
        cout<<i<<" "<<phi[i]<<endl;*/
    for(int i=1;i<N;i++)
        mu[i]+=mu[i-1],phi[i]+=phi[i-1];
}
long long Phi(long long x)
{
    if(x<N)return phi[x];
    if(P[x])return P[x];
    long long res=0;
    for(long long i=2,last;i<=x;i=last+1)
    {
        last=x/(x/i);
        res+=(last-i+1)*Phi(x/i);
    }
    return P[x]=x*(x+1)/2-res;
}
long long Mu(long long x)
{
    if(x<N)return mu[x];
    if(M[x])return M[x];
    long long res=0;
    for(long long i=2,last;i<=x;i=last+1)
    {
        last=x/(x/i);
        res+=(last-i+1)*Mu(x/i);
    }
    return M[x]=1-res;
}
int main()
{
    init();
    int t;
    cin>>t;
    while(t--)
    {
        int n;
        scanf("%d",&n);
        printf("%lld %lld\n",Phi(n),Mu(n));
    }
    return 0;
}
```
（2）模板2：求$ \sum_{i=1}^{n} \mu(i) $。
```C++
#include <bits/stdc++.h>  
typedef long long ll;  
const int INF=0x3f3f3f3f;  
const ll Maxn=1e7*2;  
const ll mod=2333333;  
using namespace std;  
int pcnt=0,prime[1300000]; // 质数  
int mu[Maxn]; // 莫比乌斯函数值  
bool vis[Maxn];  
void init()   
{  
    pcnt=0;  
    mu[1]=1;  
    for (int i=2;i<Maxn;i++)   
    {  
        if (vis[i]==0)   
        {  
            mu[i]=-1;  
            prime[++pcnt]=i;  
        }  
        for (int j=1;j<=pcnt && i*prime[j]<Maxn;j++)   
        {  
            vis[i*prime[j]]=1;  
            if (i%prime[j]!=0)  
                mu[i*prime[j]]=-mu[i];  
            else   
            {  
                mu[i*prime[j]]=0;  
                break;  
            }  
        }  
    }  
    for (int i=2;i<Maxn;i++) mu[i]+=mu[i-1];   //  函数前缀和  
}  
struct Hash   
{  
    long long key;  
    int value,next;  
} node[mod];  
int cnt=0,Head[mod]={0};  
  
void Insert(long long N,int v)   
{  // 记忆  
    int ha=N%mod;  
    node[++cnt].key=N;  
    node[cnt].value=v;  
    node[cnt].next=Head[ha];  
    Head[ha]=cnt;  
}  
int mu_pre(long long N)   
{  
    if (N<Maxn) return mu[N];  
    int ha=N%mod;  
    for (int i=Head[ha];i!=0;i=node[i].next)   
    {    
        if (node[i].key==N)   // 如果已经计算过  
            return node[i].value;  
    }  
    int ans=0;  
    for (long long i=2,j;i<=N;i=j+1)   
    {  
        j=N/(N/i);   // 分块加速  
        ans+=(j-i+1)*mu_pre(N/i);  
    }  
    Insert(N,1-ans);  // 做记忆  
    return 1-ans;  
}
```
（3）狄利克雷卷积
定义：两个数论函数f和g的卷积为$\left(f^{*} g\right)(n)=\sum_{d|n} f(g) \cdot g\left(\frac{n}{d}\right)$。
性质：满足交换律，结合律，加法分配律。
莫比尤斯函数：$\sum_{d|n} \mu(d)=[n=1]$，所以有$ \mu^{*} I=\varepsilon $，这在狄利克雷卷积中是一个很常用的恒等式。
欧拉函数：$ \sum_{d|n} \varphi(d)=n $，所以有$ \varphi^{*} I=i d $。
推论：$ \varphi^{*} I=i d \rightarrow \varphi^{*} I^{*} \mu=i d^{*} \mu \rightarrow \varphi^{*} \varepsilon=i d^{*} \mu \rightarrow \varphi=\mu^{*} i d $，所以$ \varphi(i)=\sum_{d | n} \mu(d) \cdot \frac{n}{d} $。
（4）杜教筛求，其中f是一个积性函数。
1、构造两个积性函数g,h，使得$h = f^{*}g $。
2、$ \sum_{i=1}^{n} h(i)=\sum_{i=1}^{n} \sum_{d|i} g(d) \cdot f\left(\frac{i}{d}\right)=\sum_{d=1}^{n} g(d) \sum_{i=1}^{\lfloor\frac{n}{d}\rfloor} f(i), \sum_{i=1}^{n} h(i)=\sum_{d=1}^{n} g(d) \cdot S\left(\left\lfloor\frac{n}{d}\right\rfloor\right) $。
3、$ \sum_{i=1}^{n} h(i)=g(1) \cdot S(n)+\sum_{d=2}^{n} g(d) \cdot S\left(\lfloor\frac{n}{d}\rfloor\right) \rightarrow g(1) \cdot S(n)=\sum_{i=1}^{n} h(i)-\sum_{d=2}^{n} g(d) \cdot S\left(\lfloor \frac{n}{d}\rfloor\right) $。
4、经各种分析，只要当你的h(i)的前缀和很好求，能在较短的时间内求出，那么当我们对后面的式子进行整除分块时，求S(n)的复杂度为$ O\left(n^{\frac{2}{3}}\right) $。
（5）例题1：求 $  S(n)=\sum_{i=1}^{n} \mu(i)$。
分析：构造$ g(i)=1, \quad h(i)=[i=1] $，所以$ S(n)=1-\sum_{d=2}^{n} S\left(\left\lfloor\frac{n}{d}\right\rfloor\right)$。
（6）例题2：求$ S(n)=\sum_{i=1}^{n} \varphi(i)$。
分析：构造$ g(i)=1, \quad h(i)=i$，所以$S(n)=\frac{i^{*}(i+1)}{2}-\sum_{d=2}^{n} S\left(\left\lfloor\frac{n}{d}\right\rfloor\right) $。
（7）例题3：求$S(n)=\sum_{i=1}^{n} i \cdot \varphi(i)$。
分析：构造$g(i)=i, h(i)=i^{2}$
所以$\left(f^{*} g\right)(n)=\sum_{d n} f(d) \cdot g\left(\frac{n}{d}\right)=\sum_{d | n} d \cdot \varphi(d) \cdot \frac{n}{d}=\sum_{d | n} d \cdot \varphi(d) \cdot n=n \cdot \sum_{d n} d \cdot \varphi(d)=n^{2}$
所以$S(n)=\frac{n \cdot(n+1) \cdot(2 n+1)}{6}-\sum_{d=2}^{n} d \cdot S\left(\lfloor\frac{n}{d}\rfloor\right)$。
```C++
#include <bits/stdc++.h>  
typedef long long ll;  
const int INF=0x3f3f3f3f;  
const ll Maxn=1e6*2;  
const ll mod=1e9+7;  
using namespace std;  
int pcnt=0,prime[1300000]; // 质数  
ll phi[Maxn];
ll inv6,inv2;  
bool vis[Maxn];  
void _init()   
{  
    pcnt=0;  
    phi[1]=1;  
    for (int i=2;i<Maxn;i++)   
    {  
        if (vis[i]==0)   
        {  
            phi[i]=i-1;  
            prime[++pcnt]=i;  
        }  
        for (int j=1;j<=pcnt && i*prime[j]<Maxn;j++)   
        {  
            vis[i*prime[j]]=1;  
            if (i%prime[j]!=0)  
                phi[i*prime[j]]=phi[i]*(prime[j]-1);  
            else   
            {  
                phi[i*prime[j]]=phi[i]*prime[j];  
                break;  
            }  
        }  
    }  
for (int i=2;i<Maxn;i++) phi[i]=(phi[i-1]+1ll*i*phi[i])%mod;   
//  函数前缀和  
}  
struct Hash   
{  
    ll key;  
    ll value,next;  
};
int cnt=0;  
map<int,Hash> node;
map<int,ll> Head;
  
void Insert(long long N,int v)   
{  // 记忆  
    int ha=N%mod;  
    node[++cnt].key=N;  
    node[cnt].value=v;  
    node[cnt].next=Head[ha];  
    Head[ha]=cnt;  
}  
ll sum_h(ll N) //因题而异
{
	return 1ll*N%mod*(N+1)%mod*(2*N+1)%mod*inv6; 
}
ll sum_g(ll N)  //因题而异
{
	return N%mod*(N+1)%mod*inv2%mod;
}
ll get_ans(ll N)   
{  
    if (N<Maxn) return phi[N];  
    int ha=N%mod;  
    for (int i=Head[ha];i!=0;i=node[i].next)   
    {    
        if (node[i].key==N)   // 如果已经计算过  
            return node[i].value;  
    }  
    ll ans=0;  
    for (ll i=2,j;i<=N;i=j+1)   
    {  
        j=N/(N/i);   // 分块加速 
// 求当d取i-j的贡献，所以需要求i-j的和 
        ans=(ans+1ll*(sum_g(j)-sum_g(i-1)+mod)%mod*get_ans(N/i))%mod; 
// 因题而异 
    }  
    Insert(N,(sum_h(N)-ans+mod)%mod);  // 因题而异  
    return (sum_h(N)-ans+mod)%mod;  // 因题而异
}
ll pow_mod(ll a,ll p)
{
	if (p==0) return 1;
	ll ret=pow_mod(a,p/2);
	ret=ret*ret%mod;
	if (p%2==1) ret=ret*a%mod;
	return ret;
}
int main()
{
	inv6=pow_mod(6,mod-2);
	inv2=pow_mod(2,mod-2);
	_init();
	ll n;
	int t,a,b;
	cin>>t;
	while (t--)
	{
		cin>>n>>a>>b;
		cout<<(get_ans(n)-1)*inv2%mod<<endl;
	}
	return 0;
}
```
#### **8、线性筛（$b-a<10^{6}$）**
```C++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
ll phi[2000000+5],num[2000000+5];
bool prime[4000000+5];
void deal(ll a,ll b,ll c)
{
	for (ll k=(a/c)*c;k<=b;k+=c)
	{
		if (k>=a)
		{
			phi[k-a]=phi[k-a]/c*(c-1);
			while (num[k-a]%c==0)
				num[k-a]/=c;
		}	
	}
}
void solve(ll a,ll b)
{
	memset(prime,true,sizeof(prime));
	for (ll k=a;k<=b;k++)
		phi[k-a]=k;
	for (ll k=a;k<=b;k++)
		num[k-a]=k;
	deal(a,b,2);
	for (int i=3;i<2000000;i+=2)
	{
		if (prime[i])
		{
			deal(a,b,i);
			for (int j=2*i;j<2000000;j+=i)
				prime[j]=false;
		}
	}
	for (ll k=a;k<=b;k++)
	{
		if (num[k-a]!=1)
			phi[k-a]=phi[k-a]/(num[k-a])*(num[k-a]-1);
	}
}
int main()
{
	ll a,b;
	scanf("%lld%lld",&a,&b);
	solve(a,b);
	ll ans=0;
	for (ll k=a;k<=b;k++)
		ans+=phi[k-a];
	printf("%lld\n",ans);
	return 0;
}
```
#### **9、SG定理**
（1）必胜点和必败点
P点：必败点，换而言之，就是谁处于此位置，则在双方操作正确的情况下必败。
N点：必胜点，处于此情况下，双方操作均正确的情况下必胜。
（2）必胜点和必败点的性质：
1、所有终结点是 必败点 P 。（我们以此为基本前提进行推理，换句话说，我们以此为假设）
2、从任何必胜点N 操作，至少有一种方式可以进入必败点 P。
3、无论如何操作，必败点P 都只能进入 必胜点 N。
（3）Sprague-Grundy定理（SG定理）：
游戏和的SG函数等于各个游戏SG函数的Nim和。这样就可以将每一个子游戏分而治之，从而简化了问题。而Bouton定理就是Sprague-Grundy定理在Nim游戏中的直接应用，因为单堆的Nim游戏 SG函数满足$ SG(x)=x $。
（4）SG函数
首先定义mex(minimal excludant)运算，这是施加于一个集合的运算，表示最小的不属于这个集合的非负整数。例如$ mex\{ 0,1,2,4\} =3、mex\{ 2,3,5\} =0、mex\{\}=0 $。
对于任意状态x，定义 $SG(x)=mex(S)$,其中S是x后继状态的SG函数值的集合。如x有三个后继状态分别为 $SG(a),SG(b),SG(c)$，那么$SG=mex\{SG(a),SG(b),SG(c)\}$。 这样 集合S的终态必然是空集，所以SG函数的终态为 SG(x)=0,当且仅当x为必败点P时。
注：无论游戏规则是“无法操作”算先手赢还是先手输，SG函数的终态x一定满足SG(x)=0。
（5）Anti-SG游戏定义
1、决策集合为空的操作者胜。 
2、其余规则与SG游戏一致。
（6）SJ定理
对于任意一个Anti-SG游戏，如果定义所有子游戏的SG值为0时游戏结束，先手必胜的条件： 
1、游戏的SG值为0且所有子游戏SG值均不超过1。 
2、游戏的SG值不为0且至少一个子游戏SG值超过1。
（7）模板代码
```C++
#include <bits/stdc++.h>
using namespace std;
int sg[5005];
bool mex[5005];
void init()
{
	sg[0]=0,sg[1]=0,sg[2]=1;
	for (int i=3;i<=5000;i++)
	{
		memset(mex,false,sizeof(mex));
		for (int j=0;j<=i-2;j++)
			mex[sg[j]^sg[i-2-j]]=true;
		int k=0;
		while (mex[k]) k++;
		sg[i]=k;
	}
}
int main()
{
	init();
	int t,n;
	cin>>t;
	while (t--)
	{
		cin>>n;
		if (sg[n]!=0)
			cout<<"First"<<endl;
		else cout<<"Second"<<endl;
	}
	return 0;
}
```
（8）问题描述：A和B玩一个画点小游戏，初始有一个整点三角形，两个人轮流选某个三角形内部的一个整点（不能在边界上），然后将这个点与该三角形三个顶点相连，构成三个新的三角形，一直重复如此操作，直到一方无法操作者负。A想知道在最优策略下先手必胜还是后手必胜。
数据范围：所有点的坐标小于等于18。
分析：
简单而又套路的 SG 函数，由于分割后的三角形都是独立的，所以就可以看作不同的	堆，然后就是裸的 SG 函数 + 记忆化搜索了。复杂度是$ O(n^8) $，常数非常小。
当然，记忆化还是需要优雅地记忆化，不然会 TLE。
可优化的有几个地方：
三角形平移都相同。（只需要存四维，必要。会改善复杂度，变成$ O(n^6) $。
1、三角形三点可以有序存。
2、三角形旋转 90 度相同。
3、三角形对称相同。
```C++
#include <bits/stdc++.h>  
using namespace std;  
struct Point  
{  
    int x,y;  
    Point(int x=0,int y=0):x(x),y(y){}  
    bool operator < (const Point &rhs) const  
    {  
        return x<rhs.x || (x==rhs.x && y<rhs.y);  
    }  
};  
typedef Point Vector;  
Vector operator + (Vector a,Vector b)  
{  
    return Vector(a.x+b.x,a.y+b.y);  
}  
Vector operator - (Vector a,Vector b)  
{  
    return Vector(a.x-b.x,a.y-b.y);  
}  
Vector operator * (Vector a,int p)  
{  
    return Vector (a.x*p,a.y*p);  
}  
int Cross(Vector a,Vector b)  
{  
    return a.x*b.y-a.y*b.x;  
}  
int Dot(Vector a,Vector b)  
{  
    return a.x*b.x+a.y*b.y;  
}  
int Area2(Point a,Point b,Point c)  
{  
    return abs(Cross(b-a,c-a));  
}  
int min(int a,int b,int c)  
{  
    return min(a,min(b,c));  
}  
int max(int a,int b,int c)  
{  
    return max(a,max(b,c));  
}  
Point Rotate(Point a)  
{  
    return Point(-a.y,a.x);  
}  
Point Symmetric_x(Point a)  
{  
    return Point(a.x,-a.y);  
}  
Point Symmetric_y(Point a)  
{  
    return Point(-a.x,a.y);  
}  
bool vis[19][19][19][19][19][19];  
int sg[19][19][19][19][19][19];  
int SG(Point a,Point b,Point c)  
{  
    Point point[3];  
    point[0]=a,point[1]=b,point[2]=c;  
    sort(point,point+3);  
    a=point[0],b=point[1],c=point[2];  
    for (int i=0;i<4;i++)  
    {  
        a=Rotate(a),b=Rotate(b),c=Rotate(c);  
        int left=min(a.x,b.x,c.x);  
        int down=min(a.y,b.y,c.y);  
        Vector tmp(left,down);  
        a=a-tmp,b=b-tmp,c=c-tmp;  
        if (vis[a.x][a.y][b.x][b.y][c.x][c.y]) return sg[a.x][a.y][b.x][b.y][c.x][c.y];  
    }  
    for (int i=0;i<2;i++)  
    {  
        a=Symmetric_x(a),b=Symmetric_x(b),c=Symmetric_x(c);  
        int left=min(a.x,b.x,c.x);  
        int down=min(a.y,b.y,c.y);  
        Vector tmp(left,down);  
        a=a-tmp,b=b-tmp,c=c-tmp;  
        if (vis[a.x][a.y][b.x][b.y][c.x][c.y]) return sg[a.x][a.y][b.x][b.y][c.x][c.y];  
    }  
    for (int i=0;i<2;i++)  
    {  
        a=Symmetric_y(a),b=Symmetric_y(b),c=Symmetric_y(c);  
        int left=min(a.x,b.x,c.x);  
        int down=min(a.y,b.y,c.y);  
        Vector tmp(left,down);  
        a=a-tmp,b=b-tmp,c=c-tmp;  
        if (vis[a.x][a.y][b.x][b.y][c.x][c.y]) return sg[a.x][a.y][b.x][b.y][c.x][c.y];  
    }  
    bool mex[512];  
    memset(mex,false,sizeof(mex));  
    vis[a.x][a.y][b.x][b.y][c.x][c.y]=true;  
    int &ans=sg[a.x][a.y][b.x][b.y][c.x][c.y];  
    int area=Area2(a,b,c);  
    int minx=min(a.x,b.x,c.x),maxx=max(a.x,b.x,c.x);  
    int miny=min(a.y,b.y,c.y),maxy=max(a.y,b.y,c.y);  
    for (int x=minx;x<=maxx;x++)   
        for (int y=miny;y<=maxy;y++)  
        {  
            Point p=Point(x,y);  
            int a1=Area2(p,a,b),a2=Area2(p,a,c),a3=Area2(p,b,c);  
            if (a1==0 || a2==0 || a3==0 || a1+a2+a3!=area) continue;  
            mex[SG(p,a,b)^SG(p,a,c)^SG(p,b,c)]=true;  
        }  
    int i;  
    for (i=0;;i++)  
        if (!mex[i])  
            break;  
    return ans=i;  
}  
int main()  
{  
    int t;  
    Point a,b,c;  
    memset(vis,false,sizeof(vis));  
    scanf("%d",&t);  
    while (t--)  
    {  
        scanf("%d%d%d%d%d%d",&a.x,&a.y,&b.x,&b.y,&c.x,&c.y);  
        if (SG(a,b,c)==0) printf("Second\n");  
        else printf("First\n");  
    }  
    return 0;  
}  
```
#### **10、佩尔方程**
（1）定义：不定方程$x^{2}-d y^{2}=1$ （d为正整数）称为佩尔（Pell）方程：
（2）佩尔方程求解
1、当d为完全平方数时，佩尔方程只有两解$(-1,0),(1,0)$。
2、当d不是完全平方数时，佩尔方程有无穷多解$\left(x_{1}, y_{1}\right)$
（3）当d不是完全平方数时，佩尔方程求解过程
1、暴力枚举找到一组特解。
2、求出通解
通解的矩阵形式$\left(\begin{array}{c}{\mathrm{x}_{\mathrm{n}}} \\ {\mathrm{y}_{\mathrm{n}}}\end{array}\right)=\left(\begin{array}{cc}{x_{1}} & {d y_{1}} \\ {y_{1}} & {x_{1}}\end{array}\right)^{n-1}\left(\begin{array}{l}{x_{1}} \\ {y_{1}}\end{array}\right)$，此处容易想到矩阵快速幂。
通解的递推关系式1：$\left(\begin{array}{c}{\mathrm{x}_{\mathrm{n}}} \\ {\mathrm{y}_{\mathrm{n}}}\end{array}\right)=\left(\begin{array}{c}{x_{1} x_{n-1}+d y_{1} y_{n-1}} \\ {y_{1} x_{n-1}+x_{1} y_{n-1}}\end{array}\right)$。
通解的递推关系式2：$\left(\begin{array}{c}{\mathrm{x}_{n}} \\ {y_{n}}\end{array}\right)=\left(\begin{array}{l}{2 x_{1} x_{n-1}-x_{n-2}} \\ {2 x_{1} y_{n-1}-y_{n-2}}\end{array}\right)$。
（4）算法实现：Java大整数类，C++大整数类。
#### **11、斐波那契数列性质**
（1）斐波那契递推式：
$\begin{array}{l}{F(1)=1, F(0)=0} \\ {F(n)=F(n-1)+F(n-2)}\end{array}$
斐波那契通项公式：
$F(n)=\frac{1}{\sqrt{5}}\left[\left(\frac{1+\sqrt{5}}{2}\right)^{n}-\left(\frac{1-\sqrt{5}}{2}\right)^{n}\right]$
（2）斐波那契和矩阵的关系：
我们可以构造矩阵$\left[\begin{array}{ll}{1} & {1} \\ {1} & {0}\end{array}\right]$和矩阵$\left[\begin{array}{l}{F(1)} \\ {F(0)}\end{array}\right]$
二者乘积为：
$\left[\begin{array}{cc}{1} & {1} \\ {1} & {0}\end{array}\right] *\left[\begin{array}{c}{F(1)} \\ {F(0)}\end{array}\right]=\left[\begin{array}{c}{F(1)+F(0)} \\ {F(1)}\end{array}\right]=\left[\begin{array}{c}{F(2)} \\ {F(1)}\end{array}\right]$
所以要求F(n),即求$\left[\begin{array}{ll}{1} & {1} \\ {1} & {0}\end{array}\right]^{n-1}$我们可以利用快速矩阵幂。就可以在log(n)复杂度中解决了。
（3）关于斐波那契的一些恒等式：
$\begin{array}{l}{1 : \mathrm{F}(1)+\mathrm{F}(2)+\mathrm{F}(3) \ldots+\mathrm{F}(\mathrm{n})=\mathrm{F}(\mathrm{n}+2)-1} \\ {2 : \mathrm{F}(1)^{2}+\mathrm{F}(2)^{2}+\mathrm{F}(3)^{2} \ldots+\mathrm{F}(\mathrm{n})^{2}=F(n) F(n+1)} \\ {3 : \mathrm{F}(1)+\mathrm{F}(3)+\mathrm{F}(5)+\ldots \mathrm{F}(2 \mathrm{n}-1)=\mathrm{F}(2 \mathrm{n})} \\ {4 : \mathrm{F}(2)+\mathrm{F}(4)+\mathrm{F}(6)+\ldots \mathrm{F}(2 \mathrm{n})=\mathrm{F}(2 \mathrm{n}+1)-1} \\ {5 : \mathrm{F}(\mathrm{n})=\mathrm{F}(\mathrm{m}) \mathrm{F}(\mathrm{n}-\mathrm{m}+1)+\mathrm{F}(\mathrm{m}-1) \mathrm{F}(\mathrm{n}-\mathrm{m}) \quad \text { ps: } \mathrm{n}>=\mathrm{m}} \\ {6 : \mathrm{F}(\mathrm{n}-1) \mathrm{F}(\mathrm{n}+1)=\mathrm{F}(\mathrm{n})^{2}+(-1)^{n}}\end{array}$
（4）斐波那契的数论相关：
性质1：$\operatorname{gcd}(F(n), F(m))=F(\operatorname{gcd}(n, m))$
性质2：$n|m \Leftrightarrow F(n)| F(m)$
#### **12、矩阵快速幂**
```C++
#define maxm 30
using namespace std;
struct matrix
{
    long long a[maxm][maxm];
    int row,col;
    matrix():row(maxm),col(maxm)
	{
		memset(a,0,sizeof(a));
	}
    matrix(int x, int y):row(x),col(y)
	{
		memset(a,0,sizeof(a));
	}
    long long* operator [](int x)
	{
		return a[x];	
	}
    matrix operator *(matrix x)
	{
        matrix tmp(col,x.row);
        for (int i=0;i<row;i++)
            for(int j=0;j<col;j++) 
				if(a[i][j])   //稀疏矩阵优化
                	for(int k=0;k<x.col;k++) 
						if (x[j][k])
						{
                    		tmp[i][k]+=a[i][j]*x[j][k];
                    		tmp[i][k]%=mod;
                		}
        return tmp;
    }
    void operator *=(matrix x)
	{
		*this=*this*x;
	}
    matrix operator ^(long long x)
	{
        matrix ret(row,col);
        for (int i=0;i<col;i++) ret[i][i]=1;  //单位矩阵 
        matrix tmp=*this;
        for (;x>0;x>>=1,tmp*=tmp)
			if (x&1) ret*=tmp;
        return ret;
    }
}; 
```
#### **13、类欧几里得的算法**
$\begin{array}{l}{f(a, b, c, n)=\sum_{i=0}^{n}\left\lfloor\frac{a i+b}{c}\right\rfloor} \\ {g(a, b, c, n)=\sum_{i=0}^{n}i\left\lfloor\frac{a i+b}{c}\right\rfloor} \\ {h(a, b, c, n)=\sum_{i=0}^{n}\left\lfloor\frac{a i+b}{c}\right\rfloor^{2}}\end{array}$
```C++
#include <bits/stdc++.h>
using namespace std;
const int mo=1e9+7,inv2=500000004,inv6=166666668;
typedef long long LL;
int a,b,c,l,r;
struct data
{
    int f,g,h;
};
data calc(int a,int b,int c,LL n)
{
    data tmp;
    if (!a)
    {
        tmp.f=tmp.g=tmp.h=0;
        return tmp;
    }
    if (a>=c || b>=c)
    {
        tmp=calc(a%c,b%c,c,n);
        n%=mo;
        tmp.h=(tmp.h+n*(n+1)%mo*(2*n+1)%mo*inv6%mo*(a/c)%mo*(a/c)%mo
        +(n+1)*(b/c)%mo*(b/c)%mo
        +(LL)2*(a/c)*tmp.g%mo
        +(LL)2*(b/c)*tmp.f%mo
        +n*(n+1)%mo*(a/c)%mo*(b/c))%mo;
        tmp.f=(tmp.f+n*(n+1)/2%mo*(a/c)+(n+1)*(b/c))%mo;
        tmp.g=(tmp.g+n*(n+1)%mo*(2*n+1)%mo*inv6%mo*(a/c)
		+n*(n+1)/2%mo*(b/c))%mo;
        return tmp;
    }
    LL m=((LL)a*n+b)/c;
    data nxt=calc(c,c-b-1,a,m-1);
    n%=mo; m%=mo;
    tmp.f=((n*m-nxt.f)%mo+mo)%mo;
    tmp.g=(LL)((n*(n+1)%mo*m-nxt.f-nxt.h)%mo+mo)*inv2%mo;
    tmp.h=((m*(m+1)%mo*n-(LL)2*(nxt.g+nxt.f)%mo-tmp.f)%mo+mo)%mo;
    return tmp;
}
int main()
{
    scanf("%d%d%d%d%d",&a,&c,&b,&l,&r);
    printf("%d\n",(calc(a,b,c,r).g-calc(a,b,c,l-1).g+mo)%mo);
    return 0;
}
```
（1）求$\sum_{k=1}^{N}((k M) \& M) \bmod \left(10^{9}+7\right)$。
分析：考虑kM第i位是1的有多少个，这个问题等价于$\left\lfloor\frac{k M}{2^{i}}\right\rfloor \equiv 1(\bmod 2)$。
所以$\left.\sum_{k=1}^{N}((k M) \& M)=\sum_{i=0}^{63} 2^{i}\left(\sum_{k=1}^{N} \lfloor\frac{k M}{2^{i}}\right\rfloor-\left\lfloor\frac{k M}{2^{i+1}}\right\rfloor * 2\right) |\left(m \& 2^{i}=2^{i}\right)$。
该式子可以直接用类欧几里得算法求解。
（2）注意a,b,c,d均为正数
（3）应用：求最小的x,y使得$\frac{a}{b}<\frac{y}{x}<\frac{d}{c}$成立。
分析：
若$ \frac{a}{b}>1 $，则等价于 $\frac{a%b}{b}<\frac{y}{x}-\lfloor\frac{a}{b}\rfloor < \frac{d}{c}-\lfloor\frac{a}{b}\rfloor$。
否则若$ \frac{a}{b}<1 $，则此时解出$x=y=1$。
否则有且等号不同时取到，则等价于，递归求解即可。
```C++
void f(ll a,ll b,ll c,ll d,ll& x,ll& y)
{
	if (b>a)
	{
		ll tmp=b/a;
		f(a,b%a,c,d-c*tmp,x,y);
		y+=x*tmp;
		return;
	}
	else if (d>c)
	{
		y=1;
		x=1;
		return;
	}
	else
	{
		f(d,c,b,a,y,x);
	}
}
```
#### **14、齐次线性递推式：**
问题描述：已知$F_n=\sum_{i=1}^{k}f_i*F_{n-i}$，且已知$ \{F_n\}$的前k项，求$F_n$。
注意：代码中传的参数为数列的前n项，一般若时k阶线性递推式，则传前$ pk $项$p<10 $。
```C++
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cassert>
using namespace std;
#define rep(i,a,n) for (int i=a;i<n;i++)
#define per(i,a,n) for (int i=n-1;i>=a;i--)
#define pb push_back
#define mp make_pair
#define all(x) (x).begin(),(x).end()
#define fi first
#define se second
#define SZ(x) ((int)(x).size())
typedef vector<int> VI;
typedef long long ll;
typedef pair<int,int> PII;
const ll mod=1000000007;
ll powmod(ll a,ll b) {ll res=1;a%=mod; assert(b>=0); for(;b;b>>=1){if(b&1)res=res*a%mod;a=a*a%mod;}return res;}
 
int _,n;
namespace linear_seq {
    const int N=10010;
    ll res[N],base[N],_c[N],_md[N];
 
    vector<int> Md;
    void mul(ll *a,ll *b,int k) {
        rep(i,0,k+k) _c[i]=0;
        rep(i,0,k) if (a[i]) rep(j,0,k) _c[i+j]=(_c[i+j]+a[i]*b[j])%mod;
        for (int i=k+k-1;i>=k;i--) if (_c[i])
            rep(j,0,SZ(Md)) _c[i-k+Md[j]]=(_c[i-k+Md[j]]-_c[i]*_md[Md[j]])%mod;
        rep(i,0,k) a[i]=_c[i];
    }
    int solve(ll n,VI a,VI b) {  /// a 系数 b 初值 b[n+1]=a[0]*b[n]+...
        ll ans=0,pnt=0;
        int k=SZ(a);
        assert(SZ(a)==SZ(b));
        rep(i,0,k) _md[k-1-i]=-a[i];_md[k]=1;
        Md.clear();
        rep(i,0,k) if (_md[i]!=0) Md.push_back(i);
        rep(i,0,k) res[i]=base[i]=0;
        res[0]=1;
        while ((1ll<<pnt)<=n) pnt++;
        for (int p=pnt;p>=0;p--) {
            mul(res,res,k);
            if ((n>>p)&1) {
                for (int i=k-1;i>=0;i--) res[i+1]=res[i];res[0]=0;
                rep(j,0,SZ(Md)) res[Md[j]]=(res[Md[j]]-res[k]*_md[Md[j]])%mod;
            }
        }
        rep(i,0,k) ans=(ans+res[i]*b[i])%mod;
        if (ans<0) ans+=mod;
        return ans;
    }
    VI BM(VI s) {
        VI C(1,1),B(1,1);
        int L=0,m=1,b=1;
        rep(n,0,SZ(s)) {
            ll d=0;
            rep(i,0,L+1) d=(d+(ll)C[i]*s[n-i])%mod;
            if (d==0) ++m;
            else if (2*L<=n) {
                VI T=C;
                ll c=mod-d*powmod(b,mod-2)%mod;
                while (SZ(C)<SZ(B)+m) C.pb(0);
                rep(i,0,SZ(B)) C[i+m]=(C[i+m]+c*B[i])%mod;
                L=n+1-L; B=T; b=d; m=1;
            } else {
                ll c=mod-d*powmod(b,mod-2)%mod;
                while (SZ(C)<SZ(B)+m) C.pb(0);
                rep(i,0,SZ(B)) C[i+m]=(C[i+m]+c*B[i])%mod;
                ++m;
            }
        }
        return C;
    }
    int gao(VI a,ll n) {
        VI c=BM(a);
        c.erase(c.begin());
        rep(i,0,SZ(c)) c[i]=(mod-c[i])%mod;
        return solve(n,c,VI(a.begin(),a.begin()+SZ(c)));
    }
};
 
int main() {
    for (scanf("%d",&_);_;_--) {
        scanf("%d",&n);  printf("%d\n",linear_seq::gao(VI{2,24,96,416,1536,5504,18944,64000,212992,702464},n-1));
//在gao后面尽可能的给出前n项，给的越多，结果越对
    }
}
```
#### **15、Xor专题**
（1）XOR基本性质
。
，则在正整数范围内是一个双射函数。
（2）XOR线性基
例题：
求。
分析：
该问题可以等价于一个异或线性方程组。
设为异或线性方程组的一组基，则方案数为，即。
考虑如果不在线性方程组的基中，则若取了，则方案数为。
而不在基中的数共有个，所以这部分的贡献为。
对于在基中的数共个，对剩下的数求新的基，若插入新的基后不改变基的大小，表示可以被剩下的数异或线性表出，设为剩下的数的异或线性基中自由列数，则的贡献为。
优化：
注意到，所以个数求次线性基中有大量重复计算，可以先求出个数的线性基，然后每次只插入剩下的个数。
代码：
```C++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
struct XorBase
{
	static const int Max=63;
	//pivit为线性基大小，free为自由列个数 
	int pivit,free;
	ll b[Max+5];
	XorBase()
	{
		memset(b,0,sizeof(b));
		pivit=free=0;
	}
	//若n构成新的基，则return true，否则return false 
	bool insert(ll n)
	{
	    for (int j=Max;j>=0;j--) 
		{
	        if ((n>>j)&1) 
			{
	            if(b[j]) n^=b[j];
	            else 
				{
	                b[j]=n;
	                for (int k=j-1;k>=0;k--) 
						if (b[k] && b[j]>>k&1) 
							b[j]^=b[k];
	                for (int k=j+1;k<=Max;k++) 
						if (b[k]>>j&1) 
							b[k]^=b[j];
					if (n)
					{
						pivit++;
		            	return true;	
					}
					else 
					{
						free++;
						return false;
					}
	            }
	        }
	    }
	    free++; 
	    return false;
	}
    bool check(int x)
	{
		for (int j=Max;j>=0;j--)
		{
			if ((x>>j)&1)
			{
				if (a[j]) x^=a[j];
				else return false;
			}
		}
		return true;
	}
	//线性基求交 
	friend XorBase intersection(const XorBase &a, const XorBase &b)
    {
        XorBase ans,c=b,d=b;
        for (int i=0;i<=Max;i++)
        {
            int x=a.b[i];
            if (!x) continue;
            int T=0;
            for (int j=i;j>=0;j--)
            {
                if ((x>>j)&1)
                    if (c.b[j]) 
					{ 
						x^=c.b[j]; 
						T^=d.b[j]; 
					}
                    else break;
            }
            if (!x) ans.b[i]=T;
            else 
			{ 
				c[j]=x; 
				d[j]=T; 
			}
        }
        return ans;
    }
};
const int maxn=1e5+5;
const ll mod=1e9+7;
ll a[maxn];
ll pow2[maxn];
int main()
{
	int n;
	pow2[0]=1;
	for (int i=1;i<maxn;i++)
		pow2[i]=pow2[i-1]*2%mod;
	while (~scanf("%d",&n))
	{
		XorBase base1,base2;
		base1.init();
		base2.init();
		vector <int> piv;
		for (int i=0;i<n;i++)
		{
			scanf("%lld",&a[i]);
			if (!base1.insert(a[i]))
				base2.insert(a[i]);
			else piv.push_back(i);
		}
		ll ans=0;
		if (n!=base1.pivit) 
            ans=1ll*(base1.free)*pow2[base1.free-1]%mod;
		for (int i=0;i<piv.size();i++)
		{
			XorBase base3=base2;
			for (int j=0;j<piv.size();j++)
			{
				if (i==j) continue;
				base3.insert(a[piv[j]]);
			}
			if (!base3.insert(a[piv[i]]))
				ans=(ans+pow2[base3.free-1])%mod;
		}
		printf("%lld\n",ans);
	}
	return 0;
}
```
（3）线性基区间查询
分析：考虑维护每个点的前缀线性基,线性基里将靠右的数字尽可能放高位,就是存一个额外存一个位置 p,表示这个位上的数的位置,从高位到低位扫,如果当前位置大于这个位上的位置那么交换,然后就得到了一个靠右的数字尽可能在高位的线性基
然后对于询问[l,r]在r的前缀线性基里找,只在位置大于等于 l的位更新答案。
```C++
struct XorBase
{
	static const int Max=30;
	int b[Max+5],pos[Max+5];
	void insert(int n,int p)
	{
	    for (int i=Max;i>=0;i--) 
		{
	        if ((n>>i)&1) 
			{
	            if(b[i])
				{
					if (pos[i]<p)
					{
						swap(n,b[i]);
						swap(p,pos[i]);	
					}
					n^=b[i];
				} 
	            else 
				{
	                b[i]=n;
	                pos[i]=p;
					break;
	            }
	        }
	    }
	}
	int query(int l)
	{
		int ans=0;
		for (int i=Max;i>=0;i--)
		{
			if (pos[i]>=l && (ans^b[i])>ans)
				ans^=b[i];
		}
		return ans;
	}
}b[maxn];
int main()
{
	int t,m,n,x,l,r;
	scanf("%d",&t);
	while (t--)
	{
		scanf("%d%d",&n,&m);
		for (int i=1;i<=n;i++)
		{
			scanf("%d",&x);
			b[i]=b[i-1];
			b[i].insert(x,i);
		}
		for (int i=0;i<m;i++)
		{
			scanf("%d%d",&l,&r);
			int ans=b[r].query(l);
			printf("%d\n",ans);
		}
	}
	return 0;
}
```
（4）Trie树
问题描述：给定一个数组a和一个数x，求。
```C++
struct trie
{
	int ch[maxn][2];
	int cnt[maxn],val[maxn];
	int sz;
	const static int k=31;
	trie()
	{
		memset(ch[0],0,sizeof(ch[0]));
		cnt[0]=0;
		sz=1;
	}
	void init()
	{
		memset(ch[0],0,sizeof(ch[0]));
		cnt[0]=0;
		sz=1;
	}
	int idx(int x,int i)
	{
		return (x>>i)&1;
	}
	void insert(int x)
	{
		int u=0;
		for (int i=k;i>=0;i--)
		{
			int c=idx(x,i);
			if (!ch[u][c])
			{
				memset(ch[sz],0,sizeof(ch[sz]));
				cnt[sz]=0;
				ch[u][c]=sz++;
			}
			u=ch[u][c];
			cnt[u]++;
		}
		val[u]=x;
	}
	void remove(int x)
	{
		int u=0;
		for (int i=k;i>=0;i--)
		{
			int c=idx(x,i);
			if (!ch[u][c])
				assert(false);
			u=ch[u][c];
			cnt[u]--;
		}
	}
	int find(int x)
	{
		int u=0;
		for (int i=k;i>=0;i--)
		{
			int c=idx(x,i);
			if (!ch[u][c] || cnt[ch[u][c]]==0)
			{
				assert(ch[u][c^1] && cnt[ch[u][c^1]]>0);
				u=ch[u][c^1];
			}
			else u=ch[u][c];
		}
		return val[u];
	}
};
```
（5）FWT
。
。
。
```C++
#include <iostream>
#include <cstdio>
using namespace std;
const int MOD = 998244353;
const int N = 20;
long long inv2;
long long inv(long long a,long long m)
{
	if (a==1) return 1;
	return inv(m%a,m)*(m-m/a)%m;
}
void FWT_or(int* a,int len,int opt)
{
	for (int i=1;i<len;i<<=1)
		for (int p=i<<1,j=0;j<len;j+=p)
			for (int k=0;k<i;k++)
				if (opt==1) a[i+j+k]=(a[j+k]+a[i+j+k])%MOD;
				else a[i+j+k]=(a[i+j+k]+MOD-a[j+k])%MOD;
}
void FWT_and(int* a,int len,int opt)
{
	for (int i=1;i<len;i<<=1)
		for (int p=i<<1,j=0;j<len;j+=p)
			for (int k=0;k<i;k++)
				if (op ==1) a[j+k]=(a[j+k]+a[i+j+k])%MOD;
				else a[j+k]=(a[j+k]+MOD-a[i+j+k])%MOD;
}
void FWT_xor(int* a,int len,int opt)
{
	for (int i=1;i<len;i<<=1)
		for (int p=i<<1,j=0;j<len;j+=p)
			for (int k=0;k<i;k++)
			{
				int X=a[j+k],Y=a[i+j+k];
				a[j+k]=(X+Y)%MOD;a[i+j+k]=(X+MOD-Y)%MOD;
				if (opt == -1)
					a[j+k]=1ll*a[j+k]*inv2%MOD,a[i+j+k]=1ll*a[i+j+k]*inv2%MOD;
			}
}
int a[1<<N],b[1<<N],c[1<<N];
int main()
{
	int n;
	inv2=inv(2,MOD);
	scanf("%d",&n);
	for (int i=0;i<(1<<n);i++)
		scanf("%d",&a[i]);
	for (int i=0;i<(1<<n);i++)
		scanf("%d",&b[i]);
	FWT_or(a,1<<n,1);
	FWT_or(b,1<< n,1);
	for (int i=0;i<(1<<n);i++)
		c[i]=1ll*a[i]*b[i]%MOD;
	FWT_or(c,1<<n,-1);
	for (int i=0;i<(1<<n);i++)
		printf("%d%c",c[i],i==(1<<n)-1?'\n':' ');
	FWT_or(a,1<<n,-1);
	FWT_or(b,1<<n,-1);
	FWT_and(a,1<<n,1);
	FWT_and(b,1<<n,1);
	for (int i=0;i<(1<<n);i++)
		c[i]=1ll*a[i]*b[i]%MOD;
	FWT_and(c,1<<n,-1);
	for (int i=0;i<(1<<n);i++)
		printf("%d%c",c[i],i==(1<<n)-1?'\n':' ');
	FWT_and(a,1<<n,-1);
	FWT_and(b,1<<n,-1);
	FWT_xor(a,1<<n,1);
	FWT_xor(b,1<<n,1);
	for (int i=0;i<(1<<n);i++)
		c[i]=1ll*a[i]*b[i]%MOD;
	FWT_xor(c,1<<n,-1);
	for (int i=0;i<(1<<n);i++)
		printf("%d%c",c[i],i==(1<<n)-1?'\n': ' ');
	FWT_xor(a,1<<n,-1);
	FWT_xor(b,1<<n,-1);
	return 0;
}
```
题意：
给定两个序列 A[0...2^m-1], B[0...2^m-1]
求 C[0...2^m-1]  ，满足：，m <= 19。
分析：
看C[k]的形式与集合卷积的形式接近，故转化式子时主要向普通的集合卷积式方向靠
与三种位运算都相关的结论是 ： i^j + i&j = i|j
设 x = i^j, y = i|j，则显然 k = y-x，且 k 与 x 互成关于 y 的补集，即 k = x^y
再来关心给定(x,y)，符合 x = i^j, y = i|j的(i,j)对的数目
注意到相同的位 i&j 是确定的，x = i^j 是i和j不同的位的数目，这部分谁是 0 谁是 1 不固定，故(i,j)对的数目为 2^bits(x)
此时重写原式： 。
设 
由于 [k == x^y]，第二个条件 [k == y-x] 等价于 bits(k) == bits(y) - bits(x)
所以。
将 A,B,C三个数组按 bits 划分：
。
最后按不同的维度(bits)做 FWT即可。
```C++
#include <iostream>
#include <cstring>
#include <cstdio>
using namespace std;
const int MOD=998244353;
const int N=1<<20;
int rev2;
long long inv(long long a,long long m)
{
	if (a==1) return 1;
	return inv(m%a,m)*(m-m/a)%m;
}
void FWT(int a[],int n) 
{
	for (int d=1;d<n;d<<=1)
		for (int m=d<<1,i=0;i<n;i+=m)
			for (int j=0;j<d;j++)
			{
				int x=a[i+j],y=a[i+j+d];
				a[i+j]=(x+y)%MOD;
				a[i+j+d]=(x-y+MOD) %MOD;
			}
}
void UFWT(int a[],int n) 
{
	for (int d=1;d<n;d<<=1)
		for (int m=d<<1,i=0;i<n;i+=m)
			for (int j=0;j<d;j++)
			{
				int x=a[i+j],y=a[i+j+d];
				a[i+j]=1LL*(x+y)*rev2%MOD;
				a[i+j+d]=(1LL*(x-y)*rev2%MOD+MOD)%MOD;
			}
}
int A[20][N],B[20][N],C[20][N],bit[N];
int main()
{
	int n,x;
	scanf("%d",&n);
	rev2=inv(2, MOD);
	bit[0]=0;
	for (int i=1;i<N;i++)
	{
		bit[i]=bit[i>>1]+(i & 1);
	}
	for (int i=0;i<(1<<n);i++)
	{
		scanf("%d",&x);
		A[bit[i]][i]=1ll*x*(1<<bit[i])%MOD;
	}
	for (int i=0;i<(1<<n);i++)
	{
		scanf("%d",&x);
		B[bit[i]][i]=x;
	}
	for (int i=0;i<=n;i++)
	{
		FWT(A[i],1<<n);
		FWT(B[i],1<<n);
	}
	for (int i=0;i<=n;i++)
		for (int j=i;j<=n;j++)
			for (int k=0;k<(1<<n);k++)
				C[j-i][k]=(C[j-i][k]+1ll*A[i][k]*B[j][k]%MOD)%MOD;
	for (int i=0;i<=n;i++)
		UFWT(C[i],1<<n);
	long long ans=0,base=1;
	for (int i=0;i<(1<<n);i++)
	{
		ans=(ans+C[bit[i]][i]*base%MOD)%MOD;
		base=base*1526%MOD;
	}
	printf("%lld\n",ans);
	return 0;
}
```
#### **16、Miller_Rabin判大素数**
```C++
#include <bits/stdc++.h>
using namespace std;
typedef __int128_t ll;
typedef long long LL;
const int N = 1e5 + 7;
const int times = 10;
ll fast_mod(ll a,ll b,ll mod)//计算2^q的过程
{    
    ll res = 0;    
    while(b)
    {       
        if(b & 1) res = res + a;        
        a <<= 1;        
        if(a >= mod) a -= mod;        
        if(res >= mod) res -= mod;        
        b >>= 1;    
    }    
    return res;
}
ll fast_pow_mod(ll a,ll b,ll mod)//快速幂算出a^m
{    
    ll res = 1;    
    while(b)
    {        
        if(b & 1) res = (res * a) % mod;        
        a = (a * a) % mod;        
        b >>= 1;    
    }    
    return res;
}
bool check(ll a,ll m,ll p,ll n)//对于每次随机的a进行测试
{    
    ll temp = fast_pow_mod(a,m,n),ret = temp;    
    for(int i = 0;i < p;++i)
    {        
        ret = fast_mod(temp,temp,n);        
        if(ret == 1 && temp != n - 1 && temp != 1) return true;        
        temp = ret;    
    }    
    return ret != 1;
}
bool Miller_Pabin(ll n)//Miller测试的主体结构
{   
    if(n < 2) return false;    
    if(n == 2) return true;    
    if(n & 1 == 0) return false;//对于偶数的优化    
    ll p = 0,x = n - 1;//p为Miller测试的q，x为Miller测试的m    
    while(x & 1 == 0)
    {        
        x >>= 1;        
        p++;    
    }    
    srand(time(NULL));    
    for(int i = 0;i < times;++i)
    {        
        ll o = rand() % (n - 1) + 1;//o就是Miller测试的底数a        
        if(check(o,x,p,n)) return false;   
    }    
    return true;
}
16、Java版exgcd
static ArrayList<BigInteger> gcd(BigInteger a,BigInteger b)
    {
        BigInteger ans;
        ArrayList<BigInteger> result=new ArrayList<BigInteger>();
        if (b.equals(new BigInteger("0")))
        {
            result.add(a);
            result.add(new BigInteger(String.valueOf(1)));
            result.add(new BigInteger(String.valueOf(0)));
            return result;
        }
        else
        {
            ArrayList<BigInteger> tmp=gcd(b,a.mod(b));
            ans=tmp.get(0);
            result.add(ans);
            result.add(tmp.get(2));
            result.add(tmp.get(1).subtract(tmp.get(2)
.multiply(a.divide(b))));
            return result;
        }
}
//d=result[0],x=result[1],y=result[2]
```
#### **17、计算$ a^{b^{b^{...^b}}} $**
分析：利用拓展欧拉定理。
$ a^c=\left\{
\begin{aligned}
& a^{c\ mod\ phi(m)} &  ,gcd(a,m)=1 \\
& a^{c} &  ,gcd(a,m)\neq 1,c<phi(m) \\
& a^{c\ mod\ phi(m)+phi(m)} & ,gcd(a,m)\neq 1,c \geq phi(m) 
\end{aligned}
\right.
$
```C++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const int maxn=1e7+5;
bool unprime[maxn];
int prime[maxn],phi[maxn];
void get_phi()
{
	int i,j,k=0;
	memset(unprime,false,sizeof(unprime));
	for (int i=2;i<maxn;i++)
	{
		if (!unprime[i])
		{
			prime[k++]=i;
			phi[i]=i-1;
		}
		for (j=0;j<k && prime[j]*i<maxn;j++)
		{
			unprime[prime[j]*i]=true;
			if (i%prime[j])
				phi[prime[j]*i]=phi[i]*(prime[j]-1);
			else
			{
				phi[prime[j]*i]=phi[i]*prime[j];
				break;
			} 
		}
	}
} 
int q,a,b,p;
ll k;
ll pow(ll a,ll p,ll mod)
{
	if (p==0) return 1;
	ll ret=pow(a,p/2,mod);
	ret=ret*ret%mod;
	if (p%2==1) ret=ret*a%mod;
	return ret; 
}
ll pow2(ll a,ll q)
{
	if (a==-1) return -1;
	if (q==0) return 1;
	ll ret=pow2(a,q/2);
	if (ret==-1) return -1;
	ret=ret*ret;
	if (ret>p) return -1;
	if (q%2==1) 
	{
		ret=ret*a;
		if (ret>p) return -1;
	}
	return ret;
}
ll cal(int dep)
{
	ll tmp=1;
	for (int i=0;i<k-dep+1;i++)
	{
		tmp=pow2(b,tmp);
		if (tmp==-1) return -1;
	}
	return tmp;
}
pair<ll,ll> dfs(int dep,int mod)
{
	if (mod==1) return make_pair(0,cal(dep));
	if (dep==k+1) return make_pair(1,1);
	if (dep==0)
	{
		pair<ll,ll> ret=dfs(dep+1,phi[mod]);
		ll t=ret.second,c=ret.first,cur;
		if (t!=-1)
			cur=pow2(a,t);
		else cur=-1;
		if (t==-1 || t>=phi[mod]) return make_pair(pow(a,c%phi[mod]+phi[mod],mod),cur);
		else return make_pair(pow(a,c%phi[mod],mod),cur);	
	}
	else
	{
		pair<ll,ll> ret=dfs(dep+1,phi[mod]);
		ll t=ret.second,c=ret.first,cur;
		if (t!=-1)
			cur=pow2(b,t);
		else cur=-1;
		if (t==-1 || t>=phi[mod]) return make_pair(pow(b,c%phi[mod]+phi[mod],mod),cur);
		else return make_pair(pow(b,c%phi[mod],mod),cur);	
	} 
}
int main()
{
	get_phi();	
	scanf("%d",&q);
	while (q--)
	{
		scanf("%d%d%lld%d",&a,&b,&k,&p);
		ll ans=dfs(0,p).first;
		printf("%lld\n",ans);
	}
	return 0;
}
```
#### **18、拉格朗日插值**
#### **19、数论小结**
