#include <bits/stdc++.h>
using namespace std;
int n,m;
int row[2][300000+5][4];
int col[300000+5][2][4];
char s[300000+5];
int p[4];
vector <int> a[300000+5]; 
vector <int> b[300000+5];
vector <int> c[300000+5];
int idx(char c)
{
	if (c=='A') return 0;
	if (c=='C') return 1;
	if (c=='G') return 2;
	if (c=='T') return 3;
}
char idxinv(int k)
{
	if (k==0) return 'A';
	if (k==1) return 'C';
	if (k==2) return 'G';
	if (k==3) return 'T';
}
void get_ready()
{
	for (int i=0;i<n;i++)
		for (int j=0;j<m;j++)
			for (int k=0;k<4;k++)
				if (a[i][j]!=k)
				{
					row[i%2][j][k]++;
					col[i][j%2][k]++;
				}					
} 
int best;
void solve()
{
	int ret=m*n;
	int ans=0;
	b[0][0]=p[0];
	b[0][1]=p[1];
	b[1][0]=p[2];
	b[1][1]=p[3];
	ans+=row[0][0][p[0]];
	ans+=row[1][0][p[2]];
	ans+=row[0][1][p[1]];
	ans+=row[1][1][p[3]];
	for (int i=2;i<n;i++)
		b[i][0]=b[i%2][0];
	for (int j=2;j<m;j++)
	{
		int k0=row[0][j][p[0+j%2]]+row[1][j][p[2+j%2]];
		int k1=row[0][j][p[2+j%2]]+row[1][j][p[0+j%2]];
		ans+=min(k0,k1);
		if (k0<k1)
			b[0][j]=p[0+j%2];
		else b[0][j]=p[2+j%2];
	}
	if (ans<best)
	{
		best=ans;
		for (int i=0;i<n;i++)
			c[i][0]=b[i][0];
		for (int j=0;j<m;j++)
			c[0][j]=b[0][j];
	}
	ans=0;
	ans+=col[0][0][p[0]];
	ans+=col[1][0][p[2]];
	ans+=col[0][1][p[1]];
	ans+=col[1][1][p[3]];
	for (int j=2;j<m;j++)
		b[0][j]=b[0][j%2];
	for (int i=2;i<n;i++)
	{
		int k0=col[i][0][p[0+2*(i%2)]]+col[i][1][p[1+2*(i%2)]];
		int k1=col[i][0][p[1+2*(i%2)]]+col[i][1][p[0+2*(i%2)]];
		ans+=min(k0,k1);
		if (k0<k1)
			b[i][0]=p[0+2*(i%2)];
		else b[i][0]=p[1+2*(i%2)];
	}
	if (ans<best)
	{
		best=ans;
		for (int i=0;i<n;i++)
			c[i][0]=b[i][0];
		for (int j=0;j<m;j++)
			c[0][j]=b[0][j];
	}
}
void get_ans()
{
	for (int i=1;i<n;i++)
		for (int j=1;j<m;j++)
		{
			int bit=0;
			bit|=1<<c[i-1][j-1];
			bit|=1<<c[i-1][j];
			bit|=1<<c[i][j-1];
			for (int k=0;k<4;k++)
				if ((bit|(1<<k))==15)
					c[i][j]=k;
		}
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<m;j++)
			putchar(idxinv(c[i][j]));
		putchar('\n');
	}
}
int main()
{
	scanf("%d%d",&n,&m);
	for (int i=0;i<n;i++)
	{
		a[i].resize(m);
		b[i].resize(m);
		c[i].resize(m);
	}
	for (int i=0;i<n;i++)
	{
		scanf("%s",s);
		for (int j=0;j<m;j++)
			a[i][j]=idx(s[j]);
	}
	get_ready();
	best=m*n;
	p[0]=0,p[1]=2,p[2]=1,p[3]=3;
	solve();
	for (int i=0;i<4;i++)
		p[i]=i;
	int ans=n*m;
	do
	{
		solve();
	}
	while (next_permutation(p,p+4));
	get_ans();
	return 0;
}
