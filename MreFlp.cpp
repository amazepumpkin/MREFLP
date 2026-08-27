#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <unistd.h>
#include <ctime>

using namespace std;

int n,total_n,dummy_n,c;    //n is the num of the real facilty,c is the column

struct Solution
{
    int* loc;
    int* fac;
    int value;
    void operator=(const Solution& s)
    {
        this->value=s.value;
        for(int i=0;i<total_n;++i){
            this->loc[i]=s.loc[i];
            this->fac[i]=s.fac[i];
        }
    }
};

//parameters
int m;                      //num of the rows
int** dis;             //distance between each facilty
int*  weight;               //weight of each facilty;in this case all weight equal to 1
string Infile,Outfile_result,Outfile_detail,InstanceName;
Solution curS,local_bestS,global_bestS,beginS;
double oneRunStartTime,timeLimit,find_best;
double** P_M ;             //Probility Matrix

//first local
int** quickm;               //help to quick caculate the value of a swap

//second local
int MaxIter;
int** calculatem;           //record the value of change

//InitSolutionP
double ran_c=0.20;           //the probility choose random

//UpdateP
double penalty_a=0.1,penalty_b=0.2,penalty_c=0.3; //in article b is not consistany

//SmoothP
double thre_s=0.95,smooth_s=0.30;

void Input();
void Input(int argc,char** argv);
void Output();
void AssinMemory();
void DestroryMemory();
void InitSolution();        //random choose for the first column
void InitSolution2();       //choose through the weight
void InitSolutionP();       //Init through the Probility Matrix
void FirstLocal(Solution& s); //选择一个设施一次判断
void NewFirstLocal(Solution& s); //每次判断全部都
void SecondLocal(Solution& s);
void NewSecondLocal(Solution& s); //tabu
void InitQucikm(const Solution& s);
void UpdateQucikm(int fac1,int fac2,int loc1,int loc2,const Solution& s);  //the facality exchanged and there location
void Perturbation(Solution& s,int num);
void GeneralAlgo();
void UpdateP(const Solution& s1,const Solution& s2,double** p);
void SmoothP(double** p);

void Input(int argc,char** argv)
{
    int opt=0;
    while((opt =getopt(argc,argv,"i:l:p:t:r:f:q:s:")))
    {
        switch(opt)
        {
            case 'i':
                InstanceName=string(optarg).c_str();
                break;
            case 'l':
                ran_c=atof(optarg);
                break;
            case 'p':
                MaxIter=atoi(optarg);
                break;
            case 't':
                penalty_a=atof(optarg);
                break;
            case 'r':
                penalty_b=atof(optarg);
                break;
            case 'f':
                penalty_c=atof(optarg);
                break;
            case 'q':
                thre_s=atof(optarg);
                break;
            case 's':
                smooth_s=atof(optarg);
                break;
        }
    }
}


void Input()
{
    double a;
    ifstream fin(Infile.c_str());
    if(fin.fail())
    {
        cout<<"error open"<<endl;   //可以再加入例子名称
    }
    fin>>n;
    //caculate the value of c needed
    if(m>n)
    {
        c=1;
    }else if(m<n&&n<(3*m+3)/2)
    {
        c=2;
    }else if(m==2&&n>=9)
    {
        c=ceil(double(2*n)/3)-1;
    }else if(m%2==1)
    {
        c=floor(double(2*n)/(m+1));
    }else if(m%2==0)
    {
        for(int i=1;i<100;++i)
        {
            a=m/2+1+(m+1)*i;
            if((a-m)<=n&&n<=a)
            {
                c=2*i+1;
                break;
            }
        }
    }
    total_n=c*m;
    dummy_n=total_n-n;
    AssinMemory();
    cout<<"vaildfacality:"<<n<<"  totalfacality:"<<total_n<<endl;
    cout<<"rows:"<<m<<"  column:"<<c<<endl;
    for(int i=0;i<n;++i)
        fin>>weight[i];
    
    for(int i=0;i<n;++i)
    {
        for(int j=0;j<n;++j)
        {
            fin>>dis[i][j];
            if(i>j)
                dis[i][j]=dis[j][i];
        }
    }
    
    for(int i=0;i<n;++i)
    {
        for(int j=n;j<total_n;++j)
            dis[i][j]=0;
    }
    
    for(int i=n;i<total_n;++i)
    {
        for(int j=0;j<total_n;++j)
            dis[i][j]=0;
    }
    
    fin.close();
}

void Output()
{
    FILE* fp;
    fp=fopen(Outfile_result.c_str(),"a+");
    //fprintf(fp, "%s\t%s\t%s\t\n","instance","cost","time");
    fprintf(fp, "%s\t%d\t%lf\t\n",InstanceName.c_str(),global_bestS.value,find_best);
    fclose(fp);
    fp=fopen(Outfile_detail.c_str(), "a+");
    fprintf(fp, "instance:%s n:%d total_n:%d\n",InstanceName.c_str(),n,total_n);
    for(int i=0;i<total_n;++i)
        fprintf(fp, " %d",global_bestS.loc[i]);
    fprintf(fp, "\n");
    fclose(fp);
}

void AssinMemory()
{
    dis=new int*[total_n];
    for(int i=0;i<total_n;++i)
        dis[i]=new int[total_n];
    
    weight=new int[n];
    
    curS.loc=new int[total_n];
    curS.fac=new int[total_n];
    local_bestS.loc=new int[total_n];
    local_bestS.fac=new int[total_n];
    global_bestS.loc=new int[total_n];
    global_bestS.fac=new int[total_n];
    beginS.loc=new int[total_n];
    beginS.fac=new int[total_n];
    
    quickm=new int*[total_n];
    for(int i=0;i<total_n;++i)
        quickm[i]=new int[c];
    
    calculatem=new int*[total_n];
    for(int i=0;i<total_n;++i)
        calculatem[i]=new int[c];
    
    P_M=new double*[total_n];
    for(int i=0;i<total_n;++i)
        P_M[i]=new double[c];
}

void DestroryMemory()
{
    for(int i=0;i<total_n;++i)
        delete []dis[i];
    delete []dis;
    
    delete []weight;
    
    delete []curS.loc;
    delete []curS.fac;
    delete []local_bestS.loc;
    delete []local_bestS.fac;
    delete []global_bestS.loc;
    delete []global_bestS.fac;
    delete []beginS.loc;
    delete []beginS.fac;
    
    for(int i=0;i<total_n;++i)
        delete []quickm[i];
    delete []quickm;
    
    for(int i=0;i<total_n;++i)
        delete []calculatem[i];
    delete []calculatem;
    
    for(int i=0;i<total_n;++i)
        delete []P_M[i];
    delete []P_M;
}

//greedy way
void InitSolution()
{
    int r,cl,objf=0,column,min,index;
    cl=n;
    int CL[n];
    int Obj[n];
    curS.value=0;
    for(int i=0;i<n;++i)
        CL[i]=i;
    srand((unsigned)time(NULL));
    //random choose m facality to fill the first column
    //cout<<"the 0 column: ";
    for(int i=0;i<m;++i)
    {
        r=rand()%cl;
        curS.fac[CL[r]]=i;
        curS.loc[i]=CL[r];
        cout<<CL[r]<<" ";
        cl-=1;
        swap(CL[r],CL[cl]);
    }
    cout<<endl;
    //greedy choose best m facality to the next column
    column=1;
    min=9999999;
    /*
    for(int i=0;i<n;++i)
        cout<<" "<<CL[i];
    cout<<" cl:"<<cl<<endl;  */
    while(cl>0)
    {
        //caculate the value
        if(cl>m)
        {
            for(int i=0;i<cl;++i)
            {
                for(int j=0;j<n-cl;++j)
                {
                    objf+=(column-(j/m))*dis[CL[i]][curS.loc[j]];
                    //cout<<dis[CL[i]][curS.loc[j]]<<endl;
                }
                Obj[i]=objf;
                //cout<<"facality:"<<CL[i]<<" value:"<<objf<<endl;
                objf=0;
            }
            for(int j=0;j<m;++j)
            {
                for(int i=0;i<cl;++i)
                {
                    if(Obj[i]<min)
                    {
                        min=Obj[i];
                        index=i;
                    }
                }
              //  cout<<"facality: "<<CL[index]<<"  value:"<<min<<"  ";
                curS.fac[CL[index]]=column*m+j;
                curS.loc[column*m+j]=CL[index];
                curS.value+=min;
                cl--;
                swap(CL[index],CL[cl]);
                swap(Obj[index],Obj[cl]);
                /*
                for(int i=0;i<n;++i)
                    cout<<" "<<CL[i];
                cout<<" cl:"<<cl<<endl; */
                min=9999999;
            }
            cout<<endl;
            column++;
        }else
        {
            for(int i=0;i<cl;++i)
            {
                for(int j=0;j<n-cl;++j)
                {
                    objf+=(column-(j/m))*dis[CL[i]][curS.loc[j]];
                    //cout<<dis[CL[i]][curS.loc[j]]<<endl;
                }
                curS.fac[CL[i]]=column*m+i;
                curS.loc[column*m+i]=CL[i];
                curS.value+=objf;
                objf=0;
            }
            cl=0;
        }
    }
   // cout<<"assign the dummy facality"<<endl;
    for(int i=n;i<total_n;++i)
    {
        curS.fac[i]=i;
        curS.loc[i]=i;
    }

}

//for test
void InitSolution2()
{
    int r,cl,objf=0,column,min,index;
    cl=n;
    int CL[n];
    int Obj[n];
    int weig[n];
    curS.value=0;
    for(int i=0;i<n;++i)
    {
        CL[i]=i;
        weig[i]=0;
    }
    srand((unsigned)time(NULL));
    // choose m facality with less weight to fill the first column
    for(int i=0;i<n;++i)
    {
        for(int j=0;j<n;++j)
        {
            weig[i]+=dis[i][j];
        }
    }
    min=999999;
    for(int i=0;i<m;++i)
    {
        if(weig[i]<min)
        {
            min=weig[i];
            index=i;
        }
        curS.fac[index]=i;
        curS.loc[i]=index;
        cl-=1;
        swap(CL[index],CL[cl]);
        swap(weig[index],weig[cl]);
    }
    //greedy choose best m facality to the next column
    column=1;
    min=9999999;
    /*
    for(int i=0;i<n;++i)
        cout<<" "<<CL[i];
    cout<<" cl:"<<cl<<endl;  */
    while(cl>0)
    {
        //caculate the value
        if(cl>m)
        {
            for(int i=0;i<cl;++i)
            {
                for(int j=0;j<n-cl;++j)
                {
                    objf+=(column-(j/m))*dis[CL[i]][curS.loc[j]];
                    //cout<<dis[CL[i]][curS.loc[j]]<<endl;
                }
                Obj[i]=objf;
                //cout<<"facality:"<<CL[i]<<" value:"<<objf<<endl;
                objf=0;
            }
            for(int j=0;j<m;++j)
            {
                for(int i=0;i<cl;++i)
                {
                    if(Obj[i]<min)
                    {
                        min=Obj[i];
                        index=i;
                    }
                }
              //  cout<<"facality: "<<CL[index]<<"  value:"<<min<<"  ";
                curS.fac[CL[index]]=column*m+j;
                curS.loc[column*m+j]=CL[index];
                curS.value+=min;
                cl--;
                swap(CL[index],CL[cl]);
                swap(Obj[index],Obj[cl]);
                /*
                for(int i=0;i<n;++i)
                    cout<<" "<<CL[i];
                cout<<" cl:"<<cl<<endl; */
                min=9999999;
            }
            cout<<endl;
            column++;
        }else
        {
            for(int i=0;i<cl;++i)
            {
                for(int j=0;j<n-cl;++j)
                {
                    objf+=(column-(j/m))*dis[CL[i]][curS.loc[j]];
                    //cout<<dis[CL[i]][curS.loc[j]]<<endl;
                }
                curS.fac[CL[i]]=column*m+i;
                curS.loc[column*m+i]=CL[i];
                curS.value+=objf;
                objf=0;
            }
            cl=0;
        }
    }
   // cout<<"assign the dummy facality"<<endl;
    for(int i=n;i<total_n;++i)
    {
        curS.fac[i]=i;
        curS.loc[i]=i;
    }
}

void InitSolutionP()
{
    double ran;
    int ran_l,index;
    bool isget=false;
    double best=0;
    double* pro=new double[c];
    int* num_c=new int[c];         //record the num of facility in a column
    for(int i=0;i<c;++i)
        num_c[i]=0;
    srand((unsigned)time(NULL));
    curS.value=0;
    for(int i=0;i<total_n;++i)
    {
        ran=rand()%(100)/(100.0);
        if(ran<ran_c)
        {
            //random choose a location
            ran_l=rand()%c;
            while(num_c[ran_l]==m)
                ran_l=rand()%c;
            curS.fac[i]=ran_l*m+num_c[ran_l];
            curS.loc[ran_l*m+num_c[ran_l]]=i;
            num_c[ran_l]++;
        }
        else
        {
            //choose through probility matrix
            //select the highest probility
            isget=false;
            for(int j=0;j<c;++j)
                pro[j]=P_M[i][j];
            while (isget==false) {
                best=-1;
                for(int j=0;j<c;++j)
                {
                    if(pro[j]>best)
                    {
                        best=pro[j];
                        index=j;
                    }
                }
                //add if the column isn't full
                if(num_c[index]<m)
                {
                    curS.fac[i]=index*m+num_c[index];
                    curS.loc[index*m+num_c[index]]=i;
                    num_c[index]++;
                    isget=true;
                }
                else
                {
                    pro[index]=-1;          //choose the next
                }
            }
        }
    }
    for(int i=0;i<total_n;++i)
    {
        for(int j=i+1;j<total_n;++j)
            curS.value+=abs(j/m-i/m)*dis[curS.loc[i]][curS.loc[j]];
    }
}

void InitQucikm(const Solution& s)
{
    int c_n,delta=0,c_del;
    for(int i=0;i<n;++i)
    {
        c_n=s.fac[i]/m;
        for(int j=0;j<c;++j)
        {
            if(j==c_n)
            {
                quickm[i][j]=0;
            }else
            {
                for(int k=i+1;k<n;++k)
                {
                    c_del=s.fac[k]/m;
                    delta+=(abs(j-c_del)-abs(c_n-c_del))*dis[i][k];
                }
                for(int k=0;k<i;++k)
                {
                    c_del=s.fac[k]/m;
                    delta+=(abs(j-c_del)-abs(c_n-c_del))*dis[i][k];
                }
                quickm[i][j]=delta;
                delta=0;
            }
        }
    }
    for(int i=n;i<total_n;++i)
    {
        for(int j=0;j<c;++j)
            quickm[i][j]=0;
    }
}


void UpdateQucikm(int fac1,int fac2,int loc1,int loc2,const Solution& s)
{
    int c_n,c_del,delta=0;
    for(int i=0;i<total_n;++i)
    {
        c_n=curS.fac[i]/m;
        if(i==fac1||i==fac2)
        {
            c_n=s.fac[i]/m;
            for(int j=0;j<c;++j)
            {
                if(j==c_n)
                {
                    quickm[i][j]=0;
                }else
                {
                    for(int k=i+1;k<n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs(j-c_del)-abs(c_n-c_del))*dis[i][k];
                    }
                    for(int k=0;k<i;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs(j-c_del)-abs(c_n-c_del))*dis[i][k];
                    }
                    quickm[i][j]=delta;
                    delta=0;
                }
            }
        }else
        {
            
            for(int j=0;j<c;++j)
            {
                if(j==c)
                {
                    
                }else
                {
                    quickm[i][j]=quickm[i][j]-(abs(j-loc1)-abs(c_n-loc1))*dis[i][fac1]+(abs(j-loc2)-abs(c_n-loc2))*dis[i][fac1]-(abs(j-loc2)-abs(c_n-loc2))*dis[i][fac2]+(abs(j-loc1)-abs(c_n-loc1))*dis[i][fac2];

                }
            }
        }
    }
}

//use only N1
void FirstLocal(Solution& s)
{
    bool isimprove=true;
    int minvalue,value,faca,loc,index,num=0;
    srand((unsigned)time(NULL));
    InitQucikm(s);
    minvalue=999999;
    while(isimprove)
    {
        isimprove=false;
        //random choose one facality
        faca=rand()%total_n;
        while(num<total_n&&isimprove==false)
        {
            loc=faca/m;
            //find the best exchange for that facality
            for(int i=0;i<loc*m;++i)
            {
                value=quickm[s.loc[faca]][i/m]+quickm[s.loc[i]][loc]+2*abs(i/m-loc)*dis[s.loc[faca]][s.loc[i]];
                if(value<minvalue)
                {
                    minvalue=value;
                    index=i;
                }
            }
            for(int i=(loc+1)*m;i<total_n;++i)
            {
                value=quickm[s.loc[faca]][i/m]+quickm[s.loc[i]][loc]+2*abs(i/m-loc)*dis[s.loc[faca]][s.loc[i]];
                if(value<minvalue)
                {
                    minvalue=value;
                    index=i;
                }
            }
            //if this exchange improve the solution,do it,else visit next facality
            if(minvalue<0)
            {
                curS.value+=minvalue;
                swap(s.loc[faca],s.loc[index]);
                swap(s.fac[curS.loc[faca]],s.fac[curS.loc[index]]);
                UpdateQucikm(s.loc[faca],s.loc[index],index/m,faca/m,s);
                isimprove=true;
            }else
            {
                num++;
                faca=(faca+1)%total_n;
                minvalue=9999999;
            }
        }
        num=0;
        minvalue=9999999;
    }
    
}

void NewFirstLocal(Solution& s)
{
    bool isimprove=true;
    int minvalue,value,faca,loc,index,num=0,index2;
    InitQucikm(s);
    minvalue=999999;
    while(isimprove)
    {
        isimprove=false;
        //random choose one facality
            for(int faca=0;faca<total_n;++faca)
            {
                loc=faca/m;
                for(int i=0;i<loc*m;++i)
                {
                    value=quickm[s.loc[faca]][i/m]+quickm[s.loc[i]][loc]+2*abs(i/m-loc)*dis[s.loc[faca]][s.loc[i]];
                    if(value<minvalue)
                    {
                        minvalue=value;
                        index=i;
                        index2=faca;
                    }
                }
                for(int i=(loc+1)*m;i<total_n;++i)
                {
                    value=quickm[s.loc[faca]][i/m]+quickm[s.loc[i]][loc]+2*abs(i/m-loc)*dis[s.loc[faca]][s.loc[i]];
                    if(value<minvalue)
                    {
                        minvalue=value;
                        index=i;
                        index2=faca;
                    }
                }
            }
            //find the best exchange for that facality
            //if this exchange improve the solution,do it,else visit next facality
            if(minvalue<0)
            {
                curS.value+=minvalue;
                swap(s.loc[index2],s.loc[index]);
                swap(s.fac[curS.loc[index2]],s.fac[curS.loc[index]]);
                UpdateQucikm(s.loc[index2],s.loc[index],index/m,index2/m,s);
                isimprove=true;
            }
        minvalue=9999999;
    }
}

//use both N1 and N2
void SecondLocal(Solution& s)
{
    int Noimprove=0,loc,fac_n,delta=0,loc_c,c_del;
    int loc_e,value;
    bool isexchange=true;
    int type,minvalue=999999;              //type 0:exchange,1:add
    //type=0,fac1,fac2 represent location of two facality will exchange,type=1,fac1 represent the facality ,column represent the column will add to
    int fac1,fac2,column;
    local_bestS=s;
    InitQucikm(s);
    while(Noimprove<MaxIter)
    {
        minvalue=999999;
        //caculate the value of each add(N2),find best one
        for(int i=0;i<total_n;++i)
        {
            loc=i/m;
            calculatem[s.loc[i]][loc]=0;
            
            for(int j=loc+1;j<c;++j)
            {
                loc_c=j*m;                  //the location of the facality below
                fac_n=s.loc[j*m];           //the facality that been changed
                if(fac_n<n)
                {
                    for(int k=fac_n+1;k<n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs((j-1)-c_del)-abs(j-c_del))*dis[fac_n][k];
                    }
                    for(int k=0;k<fac_n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs((j-1)-c_del)-abs(j-c_del))*dis[fac_n][k];
                }
                }
                if(s.loc[i]<n)
                {
                    for(int k=s.loc[i]+1;k<n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs(j-c_del)-abs((j-1)-c_del))*dis[s.loc[i]][k];
                    }
                    for(int k=0;k<s.loc[i];++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs(j-c_del)-abs((j-1)-c_del))*dis[s.loc[i]][k];
                    }
                }
                calculatem[s.loc[i]][j]=calculatem[s.loc[i]][j-1]+delta+2*dis[fac_n][s.loc[i]];
                if(calculatem[s.loc[i]][j]<minvalue)
                {
                    type=1;
                    minvalue=calculatem[s.loc[i]][j];
                    fac1=s.loc[i];
                    column=j;
                }
                delta=0;
            }
            
            
            for(int j=loc-1;j>0;--j)
            {
                loc_c=(j+1)*m-1;                  //the location of the facality below
                fac_n=s.loc[loc_c];           //the facality that been changed
                if(fac_n<n)
                {
                    for(int k=fac_n+1;k<n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs((j+1)-c_del)-abs(j-c_del))*dis[fac_n][k];
                    }
                    for(int k=0;k<fac_n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs((j+1)-c_del)-abs(j-c_del))*dis[fac_n][k];
                    }
                }
                
                if(s.loc[i]<n)
                {
                    for(int k=s.loc[i]+1;k<n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs(j-c_del)-abs((j+1)-c_del))*dis[s.loc[i]][k];
                    }
                    for(int k=0;k<s.loc[i];++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs(j-c_del)-abs((j+1)-c_del))*dis[s.loc[i]][k];
                    }
                }
                calculatem[s.loc[i]][j]=calculatem[s.loc[i]][j+1]+delta+2*dis[s.loc[i]][fac_n];
                if(calculatem[s.loc[i]][j]<minvalue)
                {
                    type=1;
                    minvalue=calculatem[s.loc[i]][j];
                    fac1=s.loc[i];
                    column=j;
                }
                delta=0;
            }
                calculatem[s.loc[i]][loc]=999999;
        }
        
        //caculate the value of exchange(N1)
        if(isexchange==false)
            InitQucikm(s);
        for(int i=0;i<total_n;++i)
        {
            loc_e=i/m;              //the column
            for(int j=0;j<loc_e*m;++j)
            {
                value=quickm[s.loc[i]][j/m]+quickm[s.loc[j]][loc_e]+2*abs(loc_e-j/m)*dis[s.loc[i]][s.loc[j]];
                if(value<minvalue)
                {
                    type=0;
                    minvalue=value;
                    fac1=i;
                    fac2=j;
                }
            }
            for(int j=(loc_e+1)*m;j<total_n;++j)
            {
                value=quickm[s.loc[i]][j/m]+quickm[s.loc[j]][loc_e]+2*abs(loc_e-j/m)*dis[s.loc[i]][s.loc[j]];
                if(value<minvalue)
                {
                    type=0;
                    minvalue=value;
                    fac1=i;
                    fac2=j;
                }
            }
        }
        
        //do the best
        if(type==0)
        {
            s.value+=minvalue;
            swap(s.loc[fac1],s.loc[fac2]);
            swap(s.fac[s.loc[fac1]],s.fac[s.loc[fac2]]);
            UpdateQucikm(s.loc[fac1], s.loc[fac2], fac2/m, fac1/m, s);
            isexchange=true;

            if(s.value<local_bestS.value)
            {
                Noimprove=0;
                local_bestS=s;
            }else{
                Noimprove++;
            }
        }else if(type==1)
        {
            s.value+=minvalue;
            isexchange=false;

            if(column>s.fac[fac1]/m)
            {
                for(int i=s.fac[fac1];i<column*m+1;++i)
                {
                    s.fac[s.loc[i]]--;
                    s.loc[i]=s.loc[i+1];
                }
                s.loc[column*m]=fac1;
                s.fac[fac1]=column*m;
            }else
            {
                for(int i=s.fac[fac1];i>(column+1)*m-2;--i)
                {
                    s.fac[s.loc[i]]++;
                    s.loc[i]=s.loc[i-1];
                }
                s.loc[(column+1)*m-1]=fac1;
                s.fac[fac1]=(column+1)*m-1;
            }
            
            if(s.value<local_bestS.value)
            {
                Noimprove=0;
                local_bestS=s;
            }else
            {
                Noimprove++;
            }
        }
        
    }
}

void NewSecondLocal(Solution& s)
{
    int Noimprove=0,loc,fac_n,delta=0,loc_c,c_del,iter=0,row;
    int loc_e,value;
    bool isexchange=true;
    int type,minvalue=999999;              //type 0:exchange,1:add
    //type=0,fac1,fac2 represent location of two facality will exchange,type=1,fac1 represent the facality ,column represent the column will add to
    int fac1,fac2,column;
    int tabu[total_n]; //禁忌表
    for(int i=0;i<total_n;++i)
    {
        tabu[i]=0;
    }
    local_bestS=s;
    InitQucikm(s);
    srand((unsigned)time(NULL));
    while(Noimprove<MaxIter)
    {
        iter++;
        minvalue=999999;
        //caculate the value of each add(N2),find best one
        for(int i=0;i<total_n;++i)
        {
            loc=i/m;
            calculatem[s.loc[i]][loc]=0;
            
            for(int j=loc+1;j<c;++j)
            {
                loc_c=j*m;                  //the location of the facality below
                fac_n=s.loc[j*m];           //the facality that been changed
                if(fac_n<n)
                {
                    for(int k=fac_n+1;k<n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs((j-1)-c_del)-abs(j-c_del))*dis[fac_n][k];
                    }
                    for(int k=0;k<fac_n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs((j-1)-c_del)-abs(j-c_del))*dis[fac_n][k];
                }
                }
                if(s.loc[i]<n)
                {
                    for(int k=s.loc[i]+1;k<n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs(j-c_del)-abs((j-1)-c_del))*dis[s.loc[i]][k];
                    }
                    for(int k=0;k<s.loc[i];++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs(j-c_del)-abs((j-1)-c_del))*dis[s.loc[i]][k];
                    }
                }
                calculatem[s.loc[i]][j]=calculatem[s.loc[i]][j-1]+delta+2*dis[fac_n][s.loc[i]];
#if 0
                if(calculatem[s.loc[i]][j]<minvalue)
                {
                    type=1;
                    minvalue=calculatem[s.loc[i]][j];
                    fac1=s.loc[i];
                    column=j;
                }
#endif
                if(s.value+calculatem[s.loc[i]][j]<local_bestS.value)
                {
                    type=1;
                    minvalue=calculatem[s.loc[i]][j];
                    fac1=s.loc[i];
                    column=j;
                    row=0;
                }else if(calculatem[s.loc[i]][j]<minvalue&& tabu[s.loc[i]]<iter)
                {
                    for(int aa=0;aa<m;++aa)
                    {
                        if(tabu[s.loc[j*m+aa]]<iter)
                        {
                            type=1;
                            minvalue=calculatem[s.loc[i]][j];
                            fac1=s.loc[i];
                            column=j;
                            row=aa;
                            break;
                        }
                    }
                }
                delta=0;
            }
            
            
            for(int j=loc-1;j>0;--j)
            {
                loc_c=(j+1)*m-1;                  //the location of the facality below
                fac_n=s.loc[loc_c];           //the facality that been changed
                if(fac_n<n)
                {
                    for(int k=fac_n+1;k<n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs((j+1)-c_del)-abs(j-c_del))*dis[fac_n][k];
                    }
                    for(int k=0;k<fac_n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs((j+1)-c_del)-abs(j-c_del))*dis[fac_n][k];
                    }
                }
                
                if(s.loc[i]<n)
                {
                    for(int k=s.loc[i]+1;k<n;++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs(j-c_del)-abs((j+1)-c_del))*dis[s.loc[i]][k];
                    }
                    for(int k=0;k<s.loc[i];++k)
                    {
                        c_del=s.fac[k]/m;
                        delta+=(abs(j-c_del)-abs((j+1)-c_del))*dis[s.loc[i]][k];
                    }
                }
                calculatem[s.loc[i]][j]=calculatem[s.loc[i]][j+1]+delta+2*dis[s.loc[i]][fac_n];
#if 0
                if(calculatem[s.loc[i]][j]<minvalue)
                {
                    type=1;
                    minvalue=calculatem[s.loc[i]][j];
                    fac1=s.loc[i];
                    column=j;
                }
#endif
                if(s.value+calculatem[s.loc[i]][j]<local_bestS.value)
                {
                    type=1;
                    minvalue=calculatem[s.loc[i]][j];
                    fac1=s.loc[i];
                    column=j;
                    row=0;
                }else if(calculatem[s.loc[i]][j]<minvalue&& tabu[s.loc[i]]<iter)
                {
                    for(int aa=0;aa<m;++aa)
                    {
                        if(tabu[s.loc[column*m+aa]]<iter)
                        {
                            type=1;
                            minvalue=calculatem[s.loc[i]][j];
                            fac1=s.loc[i];
                            column=j;
                            row=aa;
                            break;
                        }
                    }
                }
                delta=0;
            }
                calculatem[s.loc[i]][loc]=999999;
        }
        
        //caculate the value of exchange(N1)
        if(isexchange==false)
            InitQucikm(s);
        for(int i=0;i<total_n;++i)
        {
            loc_e=i/m;              //the column
            for(int j=0;j<loc_e*m;++j)
            {
                value=quickm[s.loc[i]][j/m]+quickm[s.loc[j]][loc_e]+2*abs(loc_e-j/m)*dis[s.loc[i]][s.loc[j]];
                if(value<minvalue)
                {
                    type=0;
                    minvalue=value;
                    fac1=i;
                    fac2=j;
                }
            }
            for(int j=(loc_e+1)*m;j<total_n;++j)
            {
                value=quickm[s.loc[i]][j/m]+quickm[s.loc[j]][loc_e]+2*abs(loc_e-j/m)*dis[s.loc[i]][s.loc[j]];
                if((value<minvalue &&tabu[s.loc[i]]<iter &&tabu[s.loc[j]]<iter)|| s.value+value<local_bestS.value)
                {
                    type=0;
                    minvalue=value;
                    fac1=i;
                    fac2=j;
                }
            }
        }
        
        //do the best
        if(type==0)
        {
            s.value+=minvalue;
            swap(s.loc[fac1],s.loc[fac2]);
            swap(s.fac[s.loc[fac1]],s.fac[s.loc[fac2]]);
            UpdateQucikm(s.loc[fac1], s.loc[fac2], fac2/m, fac1/m, s);
            isexchange=true;
            
            int tt=rand()%4;
            if(tt==0) tt=1;
            tabu[s.loc[fac1]]+=tt;
            tabu[s.loc[fac2]]+=tt;
            
            if(s.value<local_bestS.value)
            {
                Noimprove=0;
                local_bestS=s;
            }else{
                Noimprove++;
            }
        }else if(type==1)
        {
            s.value+=minvalue;
            isexchange=false;
            
            int tt=rand()%4;
            if(tt==0) tt=1;
            tabu[fac1]+=tt;
            tabu[s.loc[column*m+row]]+=tt;
            
            if(column>s.fac[fac1]/m)
            {
                for(int i=s.fac[fac1];i<column*m+1;++i)
                {
                    s.fac[s.loc[i]]--;
                    s.loc[i]=s.loc[i+1];
                }
                s.loc[column*m]=fac1;
                s.fac[fac1]=column*m;
            }else
            {
                for(int i=s.fac[fac1];i>(column+1)*m-2;--i)
                {
                    s.fac[s.loc[i]]++;
                    s.loc[i]=s.loc[i-1];
                }
                s.loc[(column+1)*m-1]=fac1;
                s.fac[fac1]=(column+1)*m-1;
            }
            
            if(s.value<local_bestS.value)
            {
                Noimprove=0;
                local_bestS=s;
            }else
            {
                Noimprove++;
            }
        }
        
    }
}

void Perturbation(Solution& s,int num)
{
    srand((unsigned)time(NULL));
    int fac1,fac2;
    InitQucikm(s);
    for(int i=0;i<num;++i)
    {
        fac1=rand()%n;
        do
        {
            fac2=rand()%n;
        }while(s.fac[fac2]>(s.fac[fac1]/m)*m-1&&s.fac[fac2]<(s.fac[fac1]/m+1)*m);
        s.value+=quickm[fac1][s.fac[fac2]/m]+quickm[fac2][s.fac[fac1]/m]+2*abs(s.fac[fac1]/m-s.fac[fac2]/m)*dis[fac1][fac2];
        swap(s.loc[s.fac[fac1]],s.loc[s.fac[fac2]]);
        swap(s.fac[fac1],s.fac[fac2]);
        UpdateQucikm(fac1, fac2, s.fac[fac2]/m, s.fac[fac1]/m, s);
    }
}

void UpdateP(const Solution& s1,const Solution& s2,double** p)
{
    for(int i=0;i<total_n;++i)
    {
        if(s1.fac[i]/m==s2.fac[i]/m)
        {
            //facility i stay in the same column
            for(int j=0;j<c;++j)
            {
                if(j==s1.fac[i]/m){
                    p[i][j]=penalty_a+(1-penalty_a)*p[i][j];
                }
                else
                {
                    p[i][j]=(1-penalty_a)*p[i][j];
                }
            }
        }
        else
        {
            for(int j=0;j<c;++j)
            {
                if(j==s1.fac[i]/m)
                {
                    p[i][j]=(1-penalty_c)*(1-penalty_b)*p[i][j];
                }
                else
                {
                   if(j==s2.fac[i]/m)
                   {
                       p[i][j]=penalty_c+(1-penalty_c)*penalty_b/(float)(c-1)+(1-penalty_c)*(1-penalty_b)*p[i][j];
                   }
                    else
                    {
                        p[i][j]=(1-penalty_c)*penalty_b/(float)(c-1)+(1-penalty_c)*(1-penalty_b)*p[i][j];
                    }
                }
            }
        }
    }
}

void SmoothP(double** p)
{
    for(int i=0;i<total_n;++i)
    {
        for(int j=0;j<c;++j)
        {
            if(p[i][j]>thre_s)
            {
                p[i][j]=(1-smooth_s)*p[i][j];
            }else if(p[i][j]<1-thre_s)
            {
                p[i][j]=smooth_s+(1-smooth_s)*p[i][j];
            }
        }
    }
}

void GeneralAlgo()
{
    for(int i=0;i<total_n;++i)
    {
        for(int j=0;j<c;++j)
            P_M[i][j]=1.0/c;
    }
    oneRunStartTime=clock();
    //InitSolution();
    //global_bestS=curS;
    global_bestS.value=999999;
    while ((clock() - oneRunStartTime) / CLOCKS_PER_SEC < timeLimit)
    {
        InitSolutionP();
        beginS=curS;
       // FirstLocal(curS);
       NewSecondLocal(curS);
        NewFirstLocal(curS);
        //first为swap decent
        //second为tabu
        
        cout<<"after local:"<<local_bestS.value<<" global:"<<global_bestS.value<<endl;
        if(local_bestS.value<global_bestS.value)
        {
            global_bestS=local_bestS;
            find_best=(clock()-oneRunStartTime)/CLOCKS_PER_SEC;
        }
        //curS=local_bestS;
        //Perturbation(curS, 15);
        //cout<<"after perturbation:"<<curS.value<<endl;
        UpdateP(beginS, local_bestS, P_M);
        SmoothP(P_M);
    }
}

int main(int argc,char ** argv) {
    
    if(argc==1)
    {
        Infile="./AnKeVa_2005_80dept_set5.txt";
        Outfile_result="./resultlearn.txt";
        Outfile_detail="./detaillearn.txt";
        InstanceName=Infile.substr(Infile.find_last_of("/")+1);
        InstanceName=InstanceName.substr(0,InstanceName.find_last_of("."));
        m=3;
    }else
    {
#if 0
        int opt=0;
        while((opt = getopt(argc,argv,"i:l:p:t:r:f:q:s:")))
        {
            switch(opt)
            {
                case 'i':
                    Infile=string(optarg).c_str();
                    break;
                case 'l':
                    ran_c=atof(optarg);
                    break;
                case 'p':
                    MaxIter=atoi(optarg);
                    break;
                case 't':
                    penalty_a=atof(optarg);
                    break;
                case 'r':
                    penalty_b=atof(optarg);
                    break;
                case 'f':
                    penalty_c=atof(optarg);
                    break;
                case 'q':
                    thre_s=atof(optarg);
                    break;
                case 's':
                    smooth_s=atof(optarg);
                    break;
            }
        }
#endif
        Infile=argv[1];
        Outfile_result=argv[2];
        Outfile_detail=argv[3];
        InstanceName=Infile.substr(Infile.find_last_of("/")+1);
        InstanceName=InstanceName.substr(0,InstanceName.find_last_of("."));
        m=atoi(argv[4]);
    }
    MaxIter=20;
    timeLimit=10;
    Input();
#if 0
    cout<<"input:"<<endl;
    for(int i=0;i<total_n;++i)
    {
        for(int j=0;j<total_n;++j)
            cout<<dis[i][j]<<" ";
        cout<<endl;
    }
#endif
#if 0
    curS.loc[0]=5;
    curS.loc[1]=3;
    curS.loc[2]=2;
    curS.loc[3]=0;
    curS.loc[4]=4;
    curS.loc[5]=7;
    curS.loc[6]=6;
    curS.loc[7]=1;
    curS.loc[8]=8;

    int vvv=0;
    for(int i=0;i<9;++i)
    {
        for(int j=i+1;j<9;++j)
        {
            vvv+=dis[curS.loc[i]][curS.loc[j]]*(j/3-i/3);
        }
    }
    cout<<"xixixix:"<<vvv<<endl;
#endif
#if 0
    cout<<endl;
    for(int i=0;i<total_n;++i)
        cout<<curS.loc[i]<<" ";
    cout<<endl;
    for(int i=0;i<total_n;++i)
        cout<<curS.fac[i]<<" "
    cout<<curS.value<<endl;
#endif

#if 0
    FirstLocal(curS);
    cout<<"after local:"<<curS.value<<endl;
    for(int i=0;i<total_n;++i)
        cout<<" "<<curS.loc[i];
    cout<<endl;
#endif

#if 0
    SecondLocal(curS);
    cout<<"after local:"<<local_bestS.value<<endl;
    for(int i=0;i<total_n;++i)
        cout<<" "<<local_bestS.loc[i];
    cout<<endl;
#endif


    GeneralAlgo();
    cout<<"value:"<<global_bestS.value<<"time:"<<find_best<<endl;
    for(int i=0;i<total_n;++i)
        cout<<" "<<global_bestS.loc[i];
    //Output();
   // Ouput();

#if 0
    InitSolution();
    cout<<"after way1:"<<curS.value<<endl;
    InitSolution2();
    cout<<"after way2:"<<curS.value<<endl;
#endif
    DestroryMemory();
    return 0;
}
