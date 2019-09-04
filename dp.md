<script type="text/javascript" async src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML"> </script>

## ACM动态规划专题
### 1、数位dp
#### （1）模板代码
``` C++
int dfs(int i, int s, bool e)   
{  
    if (i==-1) return s==target_s;  
    if (!e && ~f[i][s]) return f[i][s];  
    int res = 0;  
    int u = e?num[i]:9;  
    for (int d = first?1:0; d <= u; ++d)  
        res += dfs(i-1, new_s(s, d), e&&d==u);  
    return e?res:f[i][s]=res;  
}  
```
#### 其中：
f为记忆化数组；
i为当前处理串的第i位（权重表示法，也即后面剩下i+1位待填数）；
s为之前数字的状态（如果要求后面的数满足什么状态，也可以再记一个目标状态t之类，for的时候枚举下t）；
e表示之前的数是否是上界的前缀（即后面的数能否任意填）。
for循环枚举数字时，要注意是否能枚举0，以及0对于状态的影响，有的题目前导0和中间的0是等价的，但有的不是，对于后者可以在dfs时再加一个状态变量z，表示前面是否全部是前导0，也可以看是否是首位，然后外面统计时候枚举一下位数。It depends.
#### （2）题意：统计一个范围内数的个数，要求该数能被各位上的数整除。范围2^64。
分析：dp之前，为2520打个表，把LCM给离散化Hash。
```C++
#include <bits/stdc++.h>  
using namespace std;  
typedef long long ll;  
ll dp[20][2520][50];   
int __gcd(int a,int b)  
{  
    return b==0?a:__gcd(b,a%b);  
}  
int __lcm(int a,int b)  
{  
    if (a==0) return b;  
    if (b==0) return a;  
    return a/__gcd(a,b)*b;   
}  
int Hash[2550];  
void get_hash()  
{  
    int cnt=0;  
    for (int i=1;i<=2520;i++)  
        if (2520%i==0) Hash[i]=++cnt;  
}  
int x[20];   
ll dfs(int i,int num,int lcm,bool e)  
{  
    //e=false表示剩下i+1为可以随便填  
    if (i==-1) return num%lcm==0;  
    if (!e && dp[i][num][Hash[lcm]]>=0) return dp[i][num][Hash[lcm]];  
    ll ret=0;  
    int u=e?x[i]:9;  
    for (int d=0;d<=u;d++)  
        ret+=dfs(i-1,(num*10+d)%2520,__lcm(lcm,d),e && d==u);  
    if (!e) dp[i][num][Hash[lcm]]=ret;  
    return ret;  
}  
ll cal(ll n)  
{  
    int len=0;  
    for (ll p=n;p;p/=10) x[len++]=p%10;  
    return dfs(len-1,0,1,true);  
}  
int main()  
{  
    get_hash();  
    memset(dp,-1,sizeof(dp));  
    int t;  
    ll l,r;  
    scanf("%d",&t);  
    while (t--)  
    {  
        scanf("%I64d %I64d",&l,&r);  
        printf("%I64d\n",cal(r)-cal(l-1));  
    }  
    return 0;  
}  
```
#### （3）题意：给你一个区间，问你在这个区间内最长递增子序列长度恰为K的数有多少个。
分析：回忆LIS复杂度为nlogn的算法，用g[i]表示d值为i的最小状态编号，此题中，由于g[i]的取值为0-9，所以可以状态压缩保存。
``` C++
#include <bits/stdc++.h>  
using namespace std;  
typedef long long ll;  
int bitcount(int x)  
{  
    return x==0?0:bitcount(x>>1)+(x&1);  
}  
int newstate(int state,int k)  
{  
    /* 
        在原来的状态下插入k 
        相当于在数列最后增加k，更新g数组   
        分两种情况，原来是1,3,5,6插入4后，更新为1,3,4,6 
        原来是1,3,5,6,插入7后，更新为1,3,5,6,7  
    */   
    for (int i=k;i<=9;i++)  
    {  
        if (state&(1<<i))  
            return (state^(1<<i))|(1<<k);  
    }  
    return state|(1<<k);  
}  
int x[20],k;   
ll dp[20][1<<10][11];  
ll dfs(int i,int state,bool e,bool z)  
{  
    //e=false表示剩下i+1为可以随便填  
    if (i==-1) return bitcount(state)==k;  
    if (!e && dp[i][state][k]>=0) return dp[i][state][k];  
    ll ret=0;  
    int u=e?x[i]:9;  
    for (int d=0;d<=u;d++)  
        ret+=dfs(i-1,(z && d==0)?0:newstate(state,d),e && d==u,z && d==0);  
    if (!e) dp[i][state][k]=ret;  
    return ret;  
}  
ll cal(ll n)  
{  
    int len=0;  
    for (ll p=n;p;p/=10) x[len++]=p%10;  
    return dfs(len-1,0,true,true);  
}  
int main()  
{  
    memset(dp,-1,sizeof(dp));  
    int t,cas=0;  
    ll l,r;  
    scanf("%d",&t);  
    while (t--)  
    {  
        scanf("%I64d %I64d %d",&l,&r,&k);  
        printf("Case #%d: %I64d\n",++cas,cal(r)-cal(l-1));  
    }  
    return 0;  
}
```
#### （4）题意：求出现的数字，所有偶数出现奇数次，所有奇数出现偶数次的个数。
分析：注意前导0的处理。
```C++
#include <bits/stdc++.h>  
using namespace std;  
typedef long long ll;  
typedef unsigned long long ull;  
int x[20];  
ll dp[20][1024][32];  
//dp[i][j][k]表示长度为i，数字状态为j，偶数位状态为k，含有前导的个数   
bool vis[20][1024][32];  
bool check(int state,int evenstate)  
{  
    for (int i=0;i<=9;i++)  
    {  
        if (i%2==1 && (state&(1<<i))) return false;  
        if (i%2==0 && (evenstate&(1<<(i/2))) && !(state&(1<<i))) return false;  
    }  
    return true;  
}  
ull dfs(int i,int state,int evenstate,bool e,bool z)  
{  
    if (i==-1) return check(state,evenstate);  
    if (!e && !z && dp[i][state][evenstate]>=0)   
        return dp[i][state][evenstate];  
    ull ret=0;  
    int u=e?x[i]:9;  
    for (int d=0;d<=u;d++)  
    {  
        if (d==0 && z) ret+=dfs(i-1,state,evenstate,e && d==u,true);  
        else if (d==0 && !z) ret+=dfs(i-1,state^1,evenstate|1,e && d==u,false);  
        else if (d%2==1) ret+=dfs(i-1,state^(1<<d),evenstate,e && d==u,false);  
        else if (d%2==0) ret+=dfs(i-1,state^(1<<d),evenstate|(1<<(d/2)),e && d==u,false);   
    }  
    if (!e && !z) dp[i][state][evenstate]=ret;   
    return ret;  
}  
ull cal(ull n)  
{  
    int len=0;  
    for (ull p=n;p;p/=10) x[len++]=p%10;  
    return dfs(len-1,0,0,true,true);  
}  
int main()  
{  
    int t;  
    ull l,r;  
    memset(dp,-1,sizeof(dp));  
    memset(vis,false,sizeof(vis));  
    scanf("%d",&t);  
    while (t--)  
    {  
        scanf("%llu %llu",&l,&r);  
        printf("%llu\n",cal(r)-cal(l-1));  
    }  
    return 0;  
}  
```
#### （5）题意：统计区间中不含7，且每一位数加起来的和不是7的倍数，且这个数不是7的倍数的所有数的平方和。
分析：数位dp+统计平方和，先统计个数，再统计和，最后统计平方和。
```C++
#include <bits/stdc++.h>  
#define MOD 1000000007  
using namespace std;  
typedef long long ll;  
ll p[20],dpcnt[20][7][7],dpsum[20][7][7],dpsqr[20][7][7];  
void get_pow10()  
{  
    p[0]=1;  
    for (int i=1;i<=18;i++)  
        p[i]=p[i-1]*10%MOD;  
}  
int x[20];  
ll dfs1(int i,int sum,int mod,bool e)  
{  
    if (i==-1) return sum!=0 && mod!=0;  
    if (!e && dpcnt[i][sum][mod]>=0) return dpcnt[i][sum][mod];  
    ll ret=0;  
    int u=e?x[i]:9;  
    for (int d=0;d<=u;d++)  
    {  
        if (d==7) continue;  
        ret=(ret+dfs1(i-1,(sum+d)%7,(mod*10+d)%7,e && d==u))%MOD;  
    }  
    if (!e) dpcnt[i][sum][mod]=ret;  
    return ret;  
}  
ll dfs2(int i,int sum,int mod,ll num,bool e)  
{  
    if (i==-1)  
    {  
        if (sum!=0 && mod!=0) return num;  
        else return 0;  
    }  
    if (!e && dpsum[i][sum][mod]>=0)   
    {  
        ll ret=(num*p[i+1])%MOD;  
        ret=(ret*dpcnt[i][sum][mod])%MOD;  
        ret=(dpsum[i][sum][mod]+ret)%MOD;  
        return ret;  
    }  
    ll ret=0;  
    int u=e?x[i]:9;  
    for (int d=0;d<=u;d++)  
    {  
        if (d==7) continue;  
        ret=(ret+dfs2(i-1,(sum+d)%7,(mod*10+d)%7,(num*10+d)%MOD,e && d==u))%MOD;  
    }  
    if (!e)   
    {  
        ll tmp=(num*p[i+1])%MOD;  
        tmp=(tmp*dpcnt[i][sum][mod])%MOD;  
        tmp=((ret-tmp)%MOD+MOD)%MOD;  
        dpsum[i][sum][mod]=tmp;  
    }  
    return ret;  
}  
ll dfs3(int i,int sum,int mod,ll num,bool e)  
{  
    if (i==-1)  
    {  
        if (sum!=0 && mod!=0) return (num*num)%MOD;  
        else return 0;  
    }  
    if (!e && dpsqr[i][sum][mod]>=0)   
    {  
        ll ret=(num*p[i+1])%MOD;  
        ll tmp=(ret*ret)%MOD;  
        ll ttmp=(((2*ret)%MOD)*dpsum[i][sum][mod])%MOD;  
        ll tttmp=(tmp*dpcnt[i][sum][mod])%MOD;  
        ret=(dpsqr[i][sum][mod]+ttmp+tttmp)%MOD;  
        return ret;  
    }  
    ll ret=0;  
    int u=e?x[i]:9;  
    for (int d=0;d<=u;d++)  
    {  
        if (d==7) continue;  
        ret=(ret+dfs3(i-1,(sum+d)%7,(mod*10+d)%7,(num*10+d)%MOD,e && d==u))%MOD;  
    }  
    if (!e)   
    {  
        ll tmp=(num*p[i+1])%MOD;  
        ll ttmp=(tmp*tmp)%MOD;  
        ll tttmp=(((2*tmp)%MOD)*dpsum[i][sum][mod])%MOD;  
        tmp=(ttmp*dpcnt[i][sum][mod])%MOD;  
        tmp=((ret-tmp-tttmp)%MOD+MOD)%MOD;  
        dpsqr[i][sum][mod]=tmp;  
    }  
    return ret;  
}  
ll cal(ll n)  
{  
    int len=0;  
    for (ll p=n;p;p/=10) x[len++]=p%10;  
    dfs1(len-1,0,0,true);  
    dfs2(len-1,0,0,0,true);  
    return dfs3(len-1,0,0,0,true);  
}  
int main()  
{  
    get_pow10();  
    memset(dpcnt,-1,sizeof(dpcnt));  
    memset(dpsum,-1,sizeof(dpsum));  
    memset(dpsqr,-1,sizeof(dpsqr));  
    int t;  
    ll l,r;  
    scanf("%d",&t);  
    while (t--)  
    {  
        scanf("%lld %lld",&l,&r);  
        printf("%lld\n",(cal(r)-cal(l-1)+MOD)%MOD);  
    }  
    return 0;  
}  
```
#### （6）求$ x\in[1,a],y\in[a,b],x\&y>z $的方案数。
```C++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
ll dp1[35][2][2][2][2][2];
ll x,y,z;
int a[35],b[35],c[35];
ll dfs1(int len,bool limc,bool lima,bool limb,bool leax,bool leay)
{
    if (len==-1)
    {
        if (limc || leax || leay) return 0;
        return 1;
    }
if (dp1[len][limc][lima][limb][leax][leay]!=-1) 
return dp1[len][limc][lima][limb][leax][leay];
    int lx=(lima)?a[len]:1;
    int ly=(limb)?b[len]:1;
    ll ret=0;
    for (int i=0;i<=lx;i++)
        for (int j=0;j<=ly;j++)
        {
            if (!limc) ret+=dfs1(len-1,false,lima && i==a[len],limb && j==b[len],leax && i==0,leay && j==0);
            else
            {
                if (c[len]==1 && (i&j)==1) ret+=dfs1(len-1,true,lima && i==a[len],limb && j==b[len],leax && i==0,leay && j==0);
                if (c[len]==0 && (i&j)==1) ret+=dfs1(len-1,false,lima && i==a[len],limb && j==b[len],leax && i==0,leay && j==0);
                if (c[len]==0 && (i&j)==0) ret+=dfs1(len-1,true,lima && i==a[len],limb && j==b[len],leax && i==0,leay && j==0);
            }
        }
    dp1[len][limc][lima][limb][leax][leay]=ret;
    return ret;
}
ll solve()
{
    int lena=0;
    while (x)
    {
        a[lena++]=x%2;
        x/=2;
    }
    int lenb=0;
    while (y)
    {
        b[lenb++]=y%2;
        y/=2;
    }
    int lenc=0;
    while (z)
    {
        c[lenc++]=z%2;
        z/=2;
    }
    int len=max(lena,max(lenb,lenc));
    while (lena<len)
        a[lena++]=0;
    while (lenb<len)
        b[lenb++]=0;
    while (lenc<len)
        c[lenc++]=0;
    memset(dp1,-1,sizeof(dp1));
    return dfs1(len-1,true,true,true,true,true);
}
int main()
{
    int t;
    cin>>t;
    while (t--)
    {
        cin>>x>>y>>z;
        cout<<solve()<<endl;
    }
    return 0;
}
```
#### 注意：转移前需要判断合法性，所有dfs中的状态变量都可以记忆化。
### 2、轮廓线dp
#### 题意:在$ n \times m $的棋盘中放置$ 1 \times 1 $或$ 1 \times 2 $的骨牌，且$ 1 \times 1 $的骨牌数在$ c $到$ d $之间。输入$ n \times m$的棋盘的初始状态，求填满棋盘的方法数。
```C++
#include <iostream>  
#include <cstring>  
#include <algorithm>  
using namespace std;  
const int maxn=11;  
const int mod=1000000000+7;  
typedef long long ll;  
ll dp[2][1<<maxn][30];  
int n,m,cur;  
void update(int a,int b,int p1,int p2)  
{  
    if (b&(1<<m))  
        dp[cur][b^(1<<m)][p2]=(dp[cur][b^(1<<m)][p2]+dp[1-cur][a][p1])%mod;  
}  
char s[105][15];  
int main()  
{  
    int c,d;  
    while (scanf("%d%d%d%d",&n,&m,&c,&d)==4)  
    {  
        for (int i=0;i<n;i++)  
            scanf("%s",s[i]);  
        memset(dp,0,sizeof(dp));  
        dp[0][(1<<m)-1][0]=1;  //相当于第0行上面一行全填满  
        cur=0;  
        for (int i=0;i<n;i++)  
            for (int j=0;j<m;j++)  
            {  
                cur^=1;  
                memset(dp[cur],0,sizeof(dp[cur]));  
                if (s[i][j]=='1')  
                {  
                    for (int p=0;p<=d;p++)  
                        for (int k=0;k<(1<<m);k++)  
                        {  
                            update(k,k<<1,p,p);  //不放  
                            if (i && !(k&(1<<(m-1))))   
                                update(k,(k<<1)^(1<<m)^1,p,p);   
                                //竖放1*2  
                            if (j && !(k&1))   
                                update(k,(k<<1)^3,p,p);   
                                //横放1*2  
                            update(k,(k<<1)^1,p,p+1);    
                            //放1*1   
                        }  
                }  
                else  
                {  
                    for (int p=0;p<=d;p++)  
                        for (int k=0;k<(1<<m);k++)  
                        {  
                            update(k,(k<<1)^1,p,p);     
                            //相当于放1*1   
                        }  
                 }   
            }  
        ll ret=0;  
        for (int p=c;p<=d;p++)  
            ret=(ret+dp[cur][(1<<m)-1][p])%mod;  
        printf("%lld\n",ret);   
    }  
    return 0;  
}  
```
### 3、树形dp
#### 题意：给定一棵树，每个节点有一个价值 （ $ v_i(v_i = a_i - b_i $)），两人轮流操作，若一个人选一个节点$ u $，则另一个人只能在其儿子节点中选一个节点，先手目标是最大化所有价值之和，后手目标是最小化价值之和。根可以由先手确定，问在左右情况下，最后所有走过的节点的价值之和为多少。
分析：设根为1（即下手选择1号节点为根），$ f[u][0] $ 为以$ u $ 为根的子树中先手到达该节点的最优解， $ f[u][1] $为以$ u $为根的子树中后手到达该节点的最优解，由 minmax准则，有$ f[u][0]=max{f[v][1]}+v_u , f[u][1]=min{f[v][0]} + v_u $。
dfs一次后$ f[1][0] $为以1为根的答案。
再做一次dfs就可以得到所有情况下的答案。
```C++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const int maxn=1e5+5;
vector <int> g[maxn];
ll f[maxn][2],up[maxn][2],val[maxn];
int a[maxn],b[maxn];
int son[maxn];
void dfs(int u,int fa)
{
    f[u][0]=f[u][1]=val[u];
    ll t1=1e18,t2=-1e18;
    for (int i=0;i<g[u].size();i++)
    {
        int v=g[u][i];
        if (v==fa) continue;
        dfs(v,u);
        son[u]++;
        t1=min(t1,f[v][1]);
        t2=max(t2,f[v][0]);
    }
    if (son[u]!=0) f[u][0]+=t1;
    if (son[u]!=0) f[u][1]+=t2;
}
ll ans;
void solve(int u,int fa) 
{
    vector <int> s;
    ll temp;
    if (fa==0)
    	temp=f[u][0];
    else if (son[u]==0)
        temp=up[u][0];
    else temp=min(up[u][0],f[u][0]);
    ans=max(ans,temp);
    if (son[u]==1)
    {
        for (int i=0;i<g[u].size();i++)
        {
            int v=g[u][i];
            if (v==fa) continue;
            up[v][0]=(fa==0?0:up[u][1])+val[v];
            up[v][1]=(fa==0?0:up[u][0])+val[v];
            solve(v,u);
        }
    }
    else if (son[u]>=2)
    {
        priority_queue< pair<ll,int> ,vector < pair<ll,int> >, less< pair<ll,int> > > pq;
        for (int i=0;i<g[u].size();i++)
        {
            int v=g[u][i];
            if (v==fa) continue;
            pq.push(make_pair(f[v][0],v));
        }
        pair<ll,int> maxv=pq.top();
        pq.pop();
        pair<ll,int> semiv=pq.top();
        for (int i=0;i<g[u].size();i++)
        {
            int v=g[u][i];
            if (v==fa) continue;
            if (v==maxv.second)
                up[v][0]=max(semiv.first+val[u],up[u][1])+val[v];    
            else up[v][0]=max(maxv.first+val[u],up[u][1])+val[v];
        }
        priority_queue < pair<ll,int> ,vector < pair<ll,int> >, greater< pair<ll,int> > > qp;
        for (int i=0;i<g[u].size();i++)
        {
            int v=g[u][i];
            if (v==fa) continue;
            qp.push(make_pair(f[v][1],v));
        }
        pair<ll,int> minv=qp.top();
        qp.pop();
        semiv=qp.top();
        for (int i=0;i<g[u].size();i++)
        {
            int v=g[u][i];
            if (v==fa) continue;
            if (v==minv.second)
                up[v][1]=min(semiv.first+val[u],up[u][0])+val[v];
            else up[v][1]=min(minv.first+val[u],up[u][0])+val[v];
        }
        for (int i=0;i<g[u].size();i++)
        {
            int v=g[u][i];
            if (v==fa) continue;
            solve(v,u);
        }
    }
}
int main()
{
    int t,n,u,v;
    scanf("%d",&t);
    while (t--)
    {
        scanf("%d",&n);
        for (int i=1;i<=n;i++)
            scanf("%d",&a[i]);
        for (int i=1;i<=n;i++)
            scanf("%d",&b[i]);
        for (int i=1;i<=n;i++)
            val[i]=a[i]-b[i];
        for (int i=1;i<=n;i++)
            g[i].clear();
        for (int i=1;i<n;i++)
        {
            scanf("%d%d",&u,&v);
            g[u].push_back(v);
            g[v].push_back(u);
        }
        for (int i=1;i<=n;i++)
            son[i]=0;
        for (int i=0;i<=n;i++)
            f[i][0]=1e18;
        for (int i=0;i<=n;i++)
            f[i][1]=-1e18;
        dfs(1,0);
        up[1][0]=1e18;
        up[1][1]=-1e18;
        ans=-1e18;
        solve(1,0);
        printf("%lld\n",ans);
    }
    return 0;
}
```
#### 注意：各种细节的讨论和初始化条件。
### 4、斜率优化dp
模板：考虑dp转移式$ dp[i]=min(dp[j]+(s[i]-s[j])^2+m) $ 。
当i固定时，若j比k优，则有 $ dp[j]+(s[i]-s[j])^2+m < dp[k]+(s[i]-s[k])^2+m $，将这个式子化简得$ \frac{dp[j]+sum[j]^2-(dp[k]+sum[k]^2)}{2(sum[j]-sum[k])}<sum[i] $，
令$ f[k]=dp[k]+sum[k]^2,g[k]=2*sum[k] $，则有$ \frac{f[j]-f[k]}{g[j]-g[k]}<sum[i] $，则不等式右边是一个常数，因此可以用斜率优化dp。
```C++
//这种写法比较简单，但是需要保证k(sum[i])单调递增
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
ll f[500005],s[500005];
int n,m,Q[500005];
double Slope(int j,int k) 
{ //求斜率 
    return double((f[j]+s[j]*s[j])-(f[k]+s[k]*s[k]))/(2*s[j]-2*s[k]);
}
//注意都考虑精度问题可以写成如下版本
ll getup(int j,int k) 
{
	return (f[j]+s[j]*s[j])-(f[k]+s[k]*s[k]);
}
ll getdown(int j,int k)
{
	return (2*s[j]-2*s[k]);
}
int main() 
{
    while(scanf("%d%d",&n,&m)!=EOF) 
	{
        for(int i=1;i<=n;i++) scanf("%lld",&s[i]);
        int Left=1,Right=1;
        Q[1]=0;
        f[0]=0;
        for(int i=1;i<=n;i++) 
		{
            while (Left<Right && Slope(Q[Left],Q[Left+1])<=s[i])
            //while (Left<Right && getup(Q[Left,Q[Left+1])<=
//s[i]*getdown(Q[Left],Q[left+1]) 
				Left++; 
				//维护队首（删除非最优决策） 
            int Front=Q[Left];
            f[i]=f[Front]+(s[i]-s[Front])*(s[i]-s[Front])+m; 
			//计算当前f 
            while (Left<Right && 
Slope(Q[Right-1],Q[Right])>=Slope(Q[Right],i))
            //while (Left<Right && 
//getup(Q[Right],i)*getdown(Q[Right-1],Q[Right])<=
//getup(Q[Right-1],Q[Right])*getdown(Q[Right],i))
				Right--; 
				//维护队尾（维护下凸包性质） 
            Q[++Right]=i; //入队 
        }
        printf("%lld\n",f[n]);
    }
    return 0;
}
```
例题：$ dp[i][j]=dp[k][j-1]+prewh[i]-prewh[k]-h[k+1] \times (prew[i]-prew[k]) $ 。
（$ prew[i],prewh[i] $单调递增。）
分析：考虑从小到大枚举j，当i固定时若$ dp[k_i][j]<dp[k_2][j] $，令
$ f[k]=dp[k][j-1]-prewh[k]+h[k+1] \times prew[k],g[k]=h[k+1] $，
则有 $ \frac{f[k_2]-f[k_1]}{g[k_1]-g[k_2]}<prew[i] $ 。
```C++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const int maxn=5e3+5;
typedef __int128_t LL;
struct node 
{
	int w,h;
	bool operator < (const node &rhs) const
	{
		return h<rhs.h;
	}
}p[maxn];
ll dp[maxn][maxn];
ll prewh[maxn],prew[maxn];
LL getup(int k1,int k2,int j)
{
	LL fk1=dp[k1][j-1]-prewh[k1]+1ll*p[k1+1].h*prew[k1];
	LL fk2=dp[k2][j-1]-prewh[k2]+1ll*p[k2+1].h*prew[k2];
	return fk2-fk1;
}
LL getdown(int k1,int k2)
{
	LL gk1=p[k1+1].h;
	LL gk2=p[k2+1].h;
	return gk2-gk1;
}
int Q[maxn];
int main()
{
	int n,m;
	scanf("%d%d",&n,&m);
	for (int i=1;i<=n;i++)
		scanf("%d%d",&p[i].w,&p[i].h);
	sort(p+1,p+n+1);
	for (int i=1;i<=n;i++)
		prewh[i]=prewh[i-1]+1ll*p[i].w*p[i].h;
	for (int i=1;i<=n;i++)
		prew[i]=prew[i-1]+p[i].w;
	for (int i=1;i<=n;i++)
		dp[i][0]=1e18;
	for (int j=1;j<=m;j++)
	{
		int Left=1,Right=1;
		for(int i=1;i<=n;i++) 
		{
            while (Left<Right && 
getup(Q[Left],Q[Left+1],j)<=
prew[i]*getdown(Q[Left],Q[Left+1]))
				Left++;
            int Front=Q[Left];
            dp[i][j]=dp[Front][j-1]+prewh[i]-prewh[Front]
-1ll*p[Front+1].h*(prew[i]-prew[Front]);
            while (Left<Right && 
getup(Q[Right],i,j)*getdown(Q[Right-1],Q[Right])<=
getup(Q[Right-1],Q[Right],j)*getdown(Q[Right],i))
				Right--;
            Q[++Right]=i; 
        }
	}
	printf("%lld\n",dp[n][m]);
	return 0;
}
```
#### 特别注意：若$ \frac{f[k_2]-f[k_1]}{g[k_1]-g[k_2]}<h[i] $，其中$ h[i] $不是单调递增的函数，则需要二分。
题意：N个任务,每个任务有一个完成所需时间$ t_i $，试将这N个任务分成若干批，在每批任务开始前，机器需要启动时间S，这批任务完成所需时间为各个任务需要时间的总和(同一批任务将在同一时刻完成)。每个任务的费用是它的完成时刻乘上一个费用系数$ c_i $。求最小总费用。
分析：
设f[i]表示把前i批任务分成若干批的最小费用.运用"费用提前计算"的思想,有如下转移方程:
$ f[i]=min(f[j]+sumt[i] \times (sumc[i]-sumc[j]) + S \times (sumc[N]-sumc[j])) $。
设 $ k<j $，且k比j优，则有$ \frac{f[j]-f[k]}{sumc[j]-sumc[k]}<S+sumt[i] $。
但是t[i]不保证大于0，所以sumt[i]不是单调递增的，所以所以我们必须要维护整个队列，队头也不一定最优，每次转移时需要二分查找，找到一个位置j，j左侧线段的斜率比$ S+sumt[i] $小,j右侧线段的斜率比$ S+sumt[i] $大。
```C++
const int N=300005;
LL sumt[N],sumc[N],q[N],f[N];
inline int erfen(int i,int j,int l,int r)
{
    if (l==r)return q[l];
    while (l<r)
	{
        int mid=(l+r)>>1;
        if (f[q[mid+1]]-f[q[mid]]<=j*(sumc[q[mid+1]]-sumc[q[mid]]))
			l=mid+1;
        else r=mid;
    }
    return q[l];
}
int main()
{
    int n=read(),S=read();
    for(int i=1;i<=n;i++)
	{
        int t=read(),c=read();
        sumt[i]=sumt[i-1]+t;
        sumc[i]=sumc[i-1]+c;
    }
    int l=1,r=1;
    for(int i=1;i<=n;i++)
	{
        int j=erfen(i,S+sumt[i],l,r);
        f[i]=f[j]+(sumt[i]*sumc[i]+S*sumc[n])-(S+sumt[i])*sumc[j];
        while (l<r && 
(f[q[r]]-f[q[r-1]])*(sumc[i]-sumc[q[r]])>=(f[i]-f[q[r]])*
(sumc[q[r]]-sumc[q[r-1]]))
			r--;
        q[++r]=i;
    }
    printf("%lld\n",f[n]);
    return 0;
}
```
### 5、可逆背包问题
#### 题意：给定一个n个数的数列，每次以等概率抽取一个数（不放回），若抽出的数和大于a且小于等于b，则获胜；若大于b，则告负；否则继续抽牌。问获胜的概率。
分析：设dp[j][k]表示抽了j张牌，和为k的方案数，（$ a<b\leq500,n\leq500 $）。
枚举获胜前最后抽的一张牌，计算答案贡献复杂度为($ O(500^2) $)，但dp复杂度为($ O(500^3) $)，总复杂度为 ($ O(500^4) $)。但是可以利用“可逆背包”的做法，先预处理出dp[j][k]，然后枚举每张牌时，可以($ O(500^2) $) 撤销这张牌对dp[j][k]的贡献，再计算答案贡献。
```C++
#include <bits/stdc++.h>
using namespace std;
typedef long double ll;
ll dp[505][505];
ll f[505][505];
int p[505];
ll fac[505];
int main()
{
	int n,a,b;
	scanf("%d%d%d",&n,&a,&b);
	for (int i=0;i<n;i++)
		scanf("%d",&p[i]);
	fac[0]=1;
	for (int i=1;i<=n;i++)
		fac[i]=fac[i-1]*i;
	dp[0][0]=1;
	ll ans=0;
	for (int i=0;i<n;i++)
		for (int j=a-p[i];j>=0;j--)
			for (int k=n-1;k>=0;k--)
				dp[j+p[i]][k+1]+=dp[j][k];
	for (int i=n-1;i>=0;i--)
	{
		for (int j=0;j<=a;j++)
			for (int k=0;k<=n;k++)
			{
				if (j>=p[i] && k>=1) f[j][k]=dp[j][k]-f[j-p[i]][k-1];
				else f[j][k]=dp[j][k];
			}
		for (int j=0;j<=a;j++)
			for (int k=0;k<=n;k++)
			{
				if (j+p[i]>a && j+p[i]<=b) 
ans+=f[j][k]*fac[k]*fac[n-k-1];
			}
	}
	cout<<fixed<<setprecision(10)<<ans/fac[n]<<endl;
return 0;
} 
```
### 6、好题整理
#### 题意：钱在一个区间[a,b]内（未知），每次可以选择取x元钱： 
如果卡内余额不足x，花费b，卡内余额不动；
如果卡内余额大于x，花费a，卡内余额减少x。
问：如何操作可以花费最小的钱把钱全部取光。
分析：
可以发现，区间是没有用的，可以调整成[0,r−l]。 
但是又会发现，本来0开头的的和本来不是0开头的计算方式不同。 
于是可以分两种情况写出dp方程：
$ f[0,x]=min(max(f[0,i-1]+b,f[0,x-i]+a)|i\in[1,x]) $
$ f[1,x]=min(max(f[1,i-1]+b,f[0,x-i]+a)|i\in[1,x]) $
 
其中$ f[0,x] $表示钱的区间为[0,x]， $ f[1,x] $表示钱的区间为[a,a+x](a>0)。
注意到$ f[0,i-1]+b $单调增，$ f[0,x-i]+a $ 单调减，于是可以二分使得 $ f[0,i-1]+b $ 和$ f[0,x-i]+a $ 尽可能接近，以获取平衡的最小代价。得到一个近似最优解之后暴力在附近找。 
```C++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
ll f[200000+5],g[200000+5];
int main()
{
	int t,x,y;
	ll a,b;
	scanf("%d",&t);
	while (t--)
	{
		scanf("%d%d%lld%lld",&x,&y,&a,&b);
		int n=y-x;
		memset(f,0x3f,(n+5)*sizeof(f[0]));
		memset(g,0x3f,(n+5)*sizeof(g[0]));
		f[0]=0;
		for (int i=1;i<=n;i++)
		{
			int l=1,r=i;
			while (l<=r)
			{
				int mid=(l+r)>>1;
				if (f[mid-1]+b>f[i-mid]+a)  r=mid-1;
				else l=mid+1;
			}
			for (int k=max(l-5,1);k<=min(r+5,i);k++)
				f[i]=min(f[i],max(f[k-1]+b,f[i-k]+a));
		}
		g[0]=a;
		for (int i=1;i<=n;i++)
		{
			int l=1,r=i;
			while (l<=r)
			{
				int mid=(l+r)>>1;
				if (g[mid-1]+b>f[i-mid]+a)  r=mid-1;
				else l=mid+1;
			}
			for (int k=max(l-5,1);k<=min(r+5,i);k++)
				g[i]=min(g[i],max(g[k-1]+b,f[i-k]+a));
		}
		printf("%lld\n",x==0?f[n]:g[n]);
	}
	return 0;
}
```
