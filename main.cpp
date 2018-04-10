#include <iostream>
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <map>
#define sample 960
#define sample_test 412
#define attri 5
using namespace std;
float data[sample][attri],data2[sample_test][attri];
float pos_mean[attri-1]={0.0},neg_mean[attri-1]={0.0};
int pos_siz=0,neg_siz=0;
float S_w[attri-1][attri-1]={0.0};
pair<float,float> disc[sample];
void encode()
{
	std::ifstream file("train.txt");
	string line;
	int j=0,i;
	//cout<<"hi \n";
	for(j=0;j<sample;j++)
	{
	    int count=0;
	 	std::getline(file,line);
		std::stringstream linestream(line);
		string item;
		for(i=0;i<attri;i++)
        {

			std::getline(linestream,item,',');
			char *a = const_cast<char*>(item.c_str());
            if(i==attri-1)
                data[j][i]=float(atoi(a));
            else
                data[j][i]=float(atof(a));
        }
        if(data[j][attri-1]==1.0)
        {
           pos_siz+=1;
           for(i=0;i<attri-1;i++)
           {
               pos_mean[i]+=data[j][i];
           }
        }
        else
        {
            neg_siz+=1;
            for(i=0;i<attri-1;i++)
            {
                neg_mean[i]+=data[j][i];
            }
        }
	}
	for(i=0;i<attri-1;i++)
	{
	    pos_mean[i]=pos_mean[i]/float(pos_siz+0.0);
	    neg_mean[i]=neg_mean[i]/float(neg_siz+0.0);
	}
	std::ifstream file2("test.txt");
	for(j=0;j<sample_test;j++)
	{
	    int count=0;
	 	std::getline(file2,line);
		std::stringstream linestream(line);
		string item;
		for(i=0;i<attri;i++)
        {

			std::getline(linestream,item,',');
			char *a = const_cast<char*>(item.c_str());
            if(i==attri-1)
                data2[j][i]=float(atoi(a));
            else
                data2[j][i]=float(atof(a));
        }
	}
}
float covar(int i,int j)
{
    int m;
    float pos_cov=0.0,neg_cov=0.0;
    for(m=0;m<sample;m++)
    {
        if(data[m][attri-1]==1.0)
        {
            pos_cov=pos_cov+float((data[m][i]-pos_mean[i])*(data[m][j]-pos_mean[j]));
        }
        else
        {
            neg_cov=neg_cov+float((data[m][i]-neg_mean[i])*(data[m][j]-neg_mean[j]));
        }
    }
    return neg_cov+pos_cov;
}
void covar_matrix()
{
    int i,j;
    for(i=0;i<attri-1;i++)
    {
        for(j=0;j<attri-1;j++)
        {
            if(i<=j)
            {
                S_w[i][j]=covar(i,j);
            }
            else
            {
                S_w[i][j]=S_w[j][i];
            }
        }
    }
}
int inv_covar()
{
    int i,j,k,r,s;
    float I[attri-1][attri-1]={0.0},x;
    for(i=0;i<attri-1;i++)
        I[i][i]=1.0;
    cout<<"\n";
    for(i=0;i<attri-1;i++)
    {
        if(S_w[i][i]==0.0)
            return -1;
        x=S_w[i][i];
        for(j=0;j<attri-1;j++)
        {
            S_w[i][j]=S_w[i][j]/x;
            I[i][j]=I[i][j]/x;
        }
        for(j=0;j<attri-1;j++)
        {
            x=S_w[j][i];
            if(x!=0.0 && j!=i && x!=-0.0)
            {
                for(k=0;k<attri-1;k++)
                {
                    S_w[j][k]=S_w[j][k]+(-x*S_w[i][k]);
                    I[j][k]=I[j][k]+(-x*I[i][k]);
                }
            }
        }
    }
    for(i=0;i<attri-1;i++)
    {
        for(j=0;j<attri-1;j++)
        {
            S_w[i][j]=I[i][j];
        }
    }
    return 0;
}
void multiply(float A[attri-1][attri-1],float w_temp[attri-1],float w[attri-1])
{
    int i,j,k;
    for(i=0;i<attri-1;i++)
    {
        for(j=0;j<attri-1;j++)
        {
            w[i]+=A[i][j]*w_temp[j];
        }
    }

}
void fisher(float w[attri-1])
{
    int i,j;
    float val;
    for(i=0;i<sample;i++)
    {
        val=0.0;
        for(j=0;j<attri-1;j++)
        {
            val+=w[j]*data[i][j];
        }
        disc[i].first=val;
        disc[i].second=data[i][j];
    }
    sort(disc,disc+sample);
}
float log_two(float x)
{
    if(x==0.0 || x==1.0)
        return 0;
    else
        return log2(x);
}
float entropy(int p,int j)
{
    int b,e;
    if(p==1)
    {
        b=0;
        e=j;
    }
    else
    {
        b=j;
        e=sample-1;
    }
    float zcount=0.0,ocount=0.0,tot;
    for(;b<=e;b++)
    {
        if(disc[b].second==1.0)
            ocount+=1;
        else
            zcount+=1;
    }
    tot=zcount+ocount;
    return -(zcount/tot)*log_two(zcount/tot)-(ocount/tot)*log_two(ocount/tot);
}
int find_w0()
{
    int i,w0=0;
    float min_entro=1.0;
    for(i=0;disc[i].second!=1.0;i++);
    for(;i<sample;i++)
    {
        float temp=(((i+0.0)/sample)*entropy(1,i-1)+((sample-i+0.0)/sample)*entropy(0,i));
        if(temp<min_entro)
        {
            min_entro=temp;
            w0=i;
        }
    }
    cout<<w0<<" "<<min_entro<<" ";
    return w0;

}
void performance(float w0,float w[attri-1])
{
    int i,j,tp=0,tn=0,p=0,n=0;
    for(i=0;i<sample_test;i++)
    {
        float s=0.0;
        for(j=0;j<attri-1;j++)
        {
            s+=w[j]*data2[i][j];
        }
        if(data2[i][j]==1.0)
            p++;
        else
            n++;
        if(s>=w0 && data2[i][j]==1.0)
            tp++;
        if(s<w0 && data2[i][j]==0.0)
            tn++;
    }
    cout<<(tp+0.0)/(p+0.0)<<" "<<(tn+0.0)/(n+0.0)<<" ";
    cout<<p<<" "<<n;
}
int main()
{
    int i;
    encode();
    covar_matrix();
    for(i=0;i<attri-1;i++)
    {
        for(int j=0;j<attri-1;j++)
        {
            S_w[i][j]=S_w[i][j]/sample;
        }
    }
    int err=inv_covar();
    if(err==-1)
        cout<<"inverse does not exist";
    float w_temp[attri-1],w[attri-1]={0.0};
    for(i=0;i<attri-1;i++)
    {
        w_temp[i]=pos_mean[i]-neg_mean[i];
    }
    multiply(S_w,w_temp,w);
    fisher(w);
    int w0=find_w0();
    performance(disc[w0].first,w);
    cin>>i;
}
