#include <bits/stdc++.h>
using namespace std;
int g[1000][1000],xb[667],yb[667],cx[1000],cy[1000];
void update(int a,int b,int c)
{
	g[xb[a]][yb[a]]--;
	cx[xb[a]]--;
	cy[yb[a]]--;
	xb[a]=b;
	yb[a]=c;
	g[xb[a]][yb[a]]++;
	cx[xb[a]]++;
	cy[yb[a]]++;
}
bool canwin(int x,int y,int &nx,int &ny)
{
	for (int i=-1;i<=1;i++)
		for (int j=-1;j<=1;j++)
		{
			if (i==0 && j==0) continue;
			nx=x+i,ny=y+j;
			if (nx<1 || nx>999 || ny<1 || ny>999) continue;
			if (g[nx][ny]==0 && (cx[nx]>0 || cy[ny]>0))
				return true;
		}
	return false;
}
int main()
{
	int x,y;
	int cnt=0;
	scanf("%d%d",&x,&y);
	for (int i=1;i<=666;i++)
	{
		scanf("%d%d",&xb[i],&yb[i]);
		g[xb[i]][yb[i]]++;
		cx[xb[i]]++;
		cy[yb[i]]++;
	}
	bool finish=false;
	while (!finish)
	{
		int nx,ny;
		if (x==500 && y==500) break;
		bool win=canwin(x,y,nx,ny);
		if (!win)
		{
			int dirx=0,diry=0;
			if (x>500) dirx=-1;
			else if (x<500) dirx=1;
			if (y>500) diry=-1;
			else if (y<500) diry=1;
			nx=x+dirx,ny=y+diry;
		}
		cout<<nx<<' '<<ny<<endl;
		x=nx,y=ny;
		fflush(stdout);
		int a,b,c;
		cin>>a>>b>>c;
		if (a==-1 && b==-1 && c==-1) finish=true;
		else update(a,b,c);
	}
	int cnt1=0,cnt2=0,cnt3=0,cnt4=0;
	for (int i=1;i<=499;i++)
		for (int j=1;j<=499;j++)
			cnt1+=g[i][j];
	for (int i=1;i<=499;i++)
		for (int j=501;j<=999;j++)
			cnt2+=g[i][j];
	for (int i=501;i<=999;i++)
		for (int j=1;j<=499;j++)
			cnt3+=g[i][j];
	for (int i=501;i<=999;i++)
		for (int j=501;j<=999;j++)
			cnt4+=g[i][j];
	int minc=min(min(cnt1,cnt2),min(cnt3,cnt4));
	int dirx,diry;
	if (cnt1==minc) dirx=1,diry=1;
	if (cnt2==minc) dirx=1,diry=-1;
	if (cnt3==minc) dirx=-1,diry=1;
	if (cnt4==minc) dirx=-1,diry=-1;
	while (!finish)
	{
		if (++cnt>2000) break;
		int nx,ny;
		bool win=canwin(x,y,nx,ny);
		if (!win)
		{
			nx=x+dirx,ny=y+diry;
		}
		cout<<nx<<' '<<ny<<endl;
		x=nx,y=ny;
		fflush(stdout);
		int a,b,c;
		cin>>a>>b>>c;
		if (a==-1 && b==-1 && c==-1) finish=true;
		else update(a,b,c);
	}
	return 0;
}
