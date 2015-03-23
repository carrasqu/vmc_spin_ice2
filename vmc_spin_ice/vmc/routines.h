//
//  initialize.c
//  vmc_spin_ice
//
//  Created by Zhihao Hao on 2014-09-05.
//  Copyright (c) 2014 Zhihao Hao. All rights reserved.
//
using namespace std;
#include <iostream>     // std::cout
#include <cmath>        // std::abs

class charge_pair{
public:
    vector<int> pposeven;//position of a positive charge on the even sublattice
    vector<int> nposeven;//position of a negative charge on the even sublattice
    vector<int> pposodd;//position of a positive charge on the odd sublattice
    vector<int> nposodd;//position of a negative charge on the odd sublattice
};

void initialize(int *config,int &L,double &seedin)
{
    //give random configuration. First define a random number.
    MTRand *myrand=new MTRand();
    myrand->seed(seedin);
    int x,y,z,pyro,site,temp;
    int zpro=16*pow(L,2);
    int ypro=16*L;
    int xpro=16;
    for(z=0;z<L;z++)
    {
        for(y=0;y<L;y++)
        {
            for(x=0;x<L;x++)
            {
                for(pyro=0;pyro<4;pyro++)
                {
                    for(site=0;site<4;site++)
                    {
                        temp=myrand->randInt(1);
                        //project temp from 0,1 to -1 and 1
                        config[z*zpro+y*ypro+x*xpro+pyro*4+site]=(2*temp-1);
                    }
                }
            }
        }
    }
    return;
}

                                   
void initialize_spinice_q0(int *config,int spinon[],int &L,int &ntetra)
{
    //This subroutine we initialize a spin ice configuration at q=0
    int x,y,z,pyro,site,temp;
    int zpro=16*pow(L,2);
    int ypro=16*L;
    int xpro=16;
    for(z=0;z<L;z++)
    {
        for(y=0;y<L;y++)
        {
            for(x=0;x<L;x++)
            {
                for(pyro=0;pyro<4;pyro++)
                {
                    for(site=0;site<4;site++)
                    {
                        temp=site % 2;
                        config[z*zpro+y*ypro+x*xpro+pyro*4+site]=(2*temp-1);
                    }
                }
            }
        }
    }
    
    for(z=0;z<ntetra;z++)
    {
     spinon[z]=0;
    }
     
}

void initialize_spinice_X(int *config,int &L)
{
    //Warning!!!!!!!!!!
    //This initialize only works if L is even.
    int x,y,z,pyro,site,temp,qr;
    int zpro=16*pow(L,2);
    int ypro=16*L;
    int xpro=16;
    for(z=0;z<L;z++)
    {
        for(y=0;y<L;y++)
        {
            for(x=0;x<L;x++)
            {
                for(pyro=0;pyro<4;pyro++)
                {
                    for(site=0;site<4;site++)
                    {
                        temp=pyro*4+site;
                        qr=2*((x+y)%2)-1;
                        if((temp==0)||(temp==12)||(temp==4)||(temp==11)||(temp==5)||(temp==9)||(temp==1)||(temp==14))
                        {
                            config[z*zpro+y*ypro+x*xpro+temp]=qr;
                        }
                        else
                        {
                            config[z*zpro+y*ypro+x*xpro+temp]=-qr;
                        }
                    }
                }
            }
        }
    }
}
// Did it work?

/*obsolete code---------------------*/

//We use a routine to output the charge before and after we flip a spin.
//c1: uptetrahedron, before flipping. c2:down tetrahedron, before flipping. f1:uptetrahedron, after flipping. f2:downtetrahedron after flipping.
/* Based on Juan's tables and setup, we are reconstructing the single spin sweep codes
 ---------------------------------------*/
//we are including the new tables made by Juan! In this case, we don't need ivic...
inline int qcharge(int &t,int tetra[][4],int *config)
{
    int result=0,l;
    for(l=0;l<4;l++)
    {
        result+=config[tetra[t][l]];
    }
    //stagger it.
    if(t%2==0)
    {
        //cout<<"charge"<<result<<"\n";   
        return result;
    }
    else
    {   
        //cout<<"charge"<<-result<<"\n"; 
        return -result;
    }
}

void singlespin_update_new(int *config,int tetra[][4],int connect[][2],int &L,int &pos,double &prob,double &densitysquare, int &flag)
{
    //Again, we can operate in two modes. Thermalization mode flag=1 and updating mode, flag=0.
    int c1,c2,f1,f2,t;
    double endiff;
    //charge_new(config,tetra,connect,L,pos,c1,c2,f1,f2);
    t=connect[pos][0];
    c1=qcharge(t,tetra,config);
    //f1=c1-2*config[pos];
    t=connect[pos][1];
    c2=qcharge(t,tetra,config);
    //f2=c2+2*config[pos];
    config[pos]=-config[pos];
    t=connect[pos][0]; 
    f1=qcharge(t,tetra,config);
    t=connect[pos][1];
    f2=qcharge(t,tetra,config);
    config[pos]=-config[pos]; 
    //We note the stagered definition of charge.
    //now that we have charge before and after, we can divide into two different cases.
    if(flag)
    {
        //thermalize with proper 1/(4*2), 1/4 for spin length. 1/2 for charge square.
        endiff=(double)(pow(f1,2)+pow(f2,2)-pow(c1,2)-pow(c2,2))/8.0;
        if(endiff<=0)
        {//flip
            config[pos]=-config[pos];
        }
        else
        {
            double temp=exp(-endiff/densitysquare);
            if(temp>prob)
            {
                config[pos]=-config[pos];
            }
        }
    }
    else
    {
        //update
        //only two generic possibilities: creation/annilation a pair of monoples.
        //or the motion of a single monople.
        if((c1==0)&&(c2==0))
        {
            if(densitysquare>prob)
            {//flip if probability (here is just density) is larger than the random number
                config[pos]=-config[pos];
            }
        }
        if(((c1==2)&&(c2==-2))||((c1==-2)&&(c2==2)))
        {
            if((f1==0)&&(f2==0))
            {
                config[pos]=-config[pos];
            }
        }
        if(((c1==2)&&(c2==0))||((c1==-2)&&(c2==0)))
        {
            if(((f1==0)&&(f2==2))||((f1==0)&&(f2==-2)))
            {
                config[pos]=-config[pos];
            }
        }
        if(((c1==0)&&(c2==2))||((c1==0)&&(c2==-2)))
        {
            if(((f1==2)&&(f2==0))||((f1==-2)&&(f2==0)))
            {
                config[pos]=-config[pos];
            }
        }
    }
    return;
}

/*We are contributing two routines: 1. create a table of data using routine spinon_correlation so that we can
 model the positive and negative monopole correlator using non-interacting bosons. 
 2. create a pair of static monopoles along high symmetry directions
*/
inline double correl_compute(double &z,double &y,double&x,int&L,double &t_tilde,double &density)
{
    double pi=3.14159265359,result=0.0,rho;
    double kx,ky,kz;
    for(int nx=-L;nx<L;nx++){
        for(int ny=-L;ny<L;ny++){
            for(int nz=-L;nz<L;nz++){
                //compute momentum
                kx=2*pi*((double)nx)/((double)(L));
                ky=2*pi*((double)ny)/((double)(L));
                kz=2*pi*((double)nz)/((double)(L));
                //compute rho factor
                rho=4.0*(cos(kx/2.0)*cos(ky/2.0)+cos(kx/2.0)*cos(kz/2.0)+cos(ky/2.0)*cos(kz/2.0));
                //accumulate
                result+=cos(kx*x+ky*y+kz*z)*0.5*t_tilde*rho/pow(1.0-t_tilde*rho,0.5);
            }
        }
    }
    //divide by two factors of twos. The first one is from the reciprocal space normalization. the second one is from the formula to compute the correlation.
    //log the full corrlation function 
    return log(result/((double)pow(L,3)*2.0*2.0)+density);
}

void spinon_correlation(double correl[][16],double &t_tilde,double &eta,int &L,double &density)
{
    //correl is the holder for the table, t_tilde/4 is the hopping amplitude of spinons, eta is the overlap between the RK state and the interacting ground state. eta is usually approximated to be 1.
    //double pi=3.14159265359;
    int x,y,z;
    double ztemp,ytemp,xtemp;
    int zprod=pow(L,2);
    int mu_1,mu_2;
    //define an array storing all the displacement vectors
    double disvec[4][3]={{0.0,0.0,0.0},{0.0,0.5,0.5},{0.5,0.0,0.5},{0.5,0.5,0.0}};
    double temp;
    int indtemp;
    for(z=0;z<L;z++){
        for(y=0; y<L; y++) {
            for(x=0;x<L;x++){
                indtemp=z*zprod+y*L+x;
                //the four contributions: the first contribution: diagonal. 00,11,22,33. The "displacement is just x\hat{x}+y\hat{y}+\z\hat{z}
                for(mu_1=0;mu_1<4;mu_1++){
                    for(mu_2=0;mu_2<4;mu_2++){
                        //determine the real dispacement vector
                        ztemp=(double)z+disvec[mu_1][2]-disvec[mu_2][2];
                        ytemp=(double)y+disvec[mu_1][1]-disvec[mu_2][1];
                        xtemp=(double)x+disvec[mu_1][0]-disvec[mu_2][0];
                        //determine the correlator
                        temp=eta*correl_compute(ztemp,ytemp,xtemp,L,t_tilde,density);
                        //stock the table
                        correl[indtemp][mu_1*4+mu_2]=temp;
                    }
                }
            }
        }
    }
    return;
}
//Let us use one more simple routine to compute indtemp from x1,y1,z1,x2,y2,z2
int correl_index(int &t1,int &t2,int&L,int &flag)
{
    int zpro=pow(L,2);
    flag=0;
    int x1,y1,z1,x2,y2,z2;
    z1=t1/zpro;
    y1=(t1-z1*zpro)/L;
    x1=t1-zpro*z1-y1*L;
    z2=t2/zpro;
    y2=(t2-z2*zpro)/L;
    x2=t2-z2*zpro-y2*L;
    //8 different cases
    if((z1>=z2)&&(y1>=y2)&&(x1>=x2))
    {
        return (z1-z2)*zpro+(y1-y2)*L+(x1-x2);
    }
    else if((z1<z2)&&(y1>=y2)&&(x1>=x2))
    {
        return (z1-z2+L)*zpro+(y1-y2)*L+(x1-x2);
    }
    else if((z1>=z2)&&(y1<y2)&&(x1>=x2))
    {
        return (z1-z2)*zpro+(y1-y2+L)*L+(x1-x2);
    }
    else if((z1>=z2)&&(y1>=y2)&&(x1<x2))
    {
        return (z1-z2)*zpro+(y1-y2)*L+(x1-x2+L);
    }
    //the next four cases, the "2" is in front of 1.
    else if((z1<z2)&&(y1<y2)&&(x1>=x2))
    {
        flag=1;
        return (z2-z1)*zpro+(y2-y1)*L+(x2-x1+L)%L;
    }
    else if((z1>=z2)&&(y1<y2)&&(x1<x2))
    {
        flag=1;
        return ((z2-z1+L)%L)*zpro+(y2-y1)*L+(x2-x1);
    }
    else if((z1<z2)&&(y1>=y2)&&(x1<x2))
    {
        flag=1;
        return (z2-z1)*zpro+((y2-y1+L)%L)*L+(x2-x1);
    }
    else
    {
        flag=1;
        return (z2-z1)*zpro+(y2-y1)*L+(x2-x1);
    }
}


//we need to define more functions which take a special spinon configuration, charge_pair, to compute the would-be coefficient of the configuration
//amp is the computed amplitude


/* we also need code to attempt a pair spin flip. Is this better/needed?
 we generate the seed from a random process in the main program*/
void pair_flip(int *config,int ivic[][6],int tetra[][4],int connect[][2],int &L,double &densitysquare,int&pos,int &pos2,double &prob,charge_pair&chargepairs,double correl[][16],MTRand *myrand)
{
    //pos2 will be a random number from 0 to 5.
    //We choose a random position which opposite of pos.
  
    int flag=0;
//    
//    cout<<"config"<< "\n ";  
//    nhh=pow(L,3)*16;
//
//    for(l=0;l<nhh;l++)
//    {
//    cout<<config[l]<<"\n";
//      }
//    cout<<"-------------------"<<"\n";

    //charge_pair cpairtemp;
    pos2=ivic[pos][pos2];
    if(config[pos]*config[pos2]==1)
    {
       // cout<<"nospin flip"<<"\n";
        return;
    }
    int t1,t2;
    t1=connect[pos][0];
    t2=connect[pos][1];
    //repeated tetrahedra are deleted.
    if(t1==connect[pos2][0])
    {
        //In this case, the pair of spins connect two tetrahedron on "odd sublattice"
        t1=connect[pos2][1];
        flag=1;
    }
    else if(t2==connect[pos2][1])
    {
        //the pair of spins connect two tetrahedron on "even sublattice"
        t2=connect[pos2][0];
        flag=0;
    }
    //
    int c1,c2,f1,f2;
    c1=qcharge(t1,tetra,config);
    c2=qcharge(t2,tetra,config);
    
    config[pos]=-config[pos];
    config[pos2]=-config[pos2];
    f1=qcharge(t1,tetra,config);
    f2=qcharge(t2,tetra,config);
    config[pos]=-config[pos];
    config[pos2]=-config[pos2];
    

 
    
//    int temp=t1%2;
//    if(temp==0)
//    {
//       f1=c1-2*config[pos];
//       f2=c2-2*config[pos2];
//    }
//    else
//    {
//       f1=c1+2*config[pos];
//       f2=c2+2*config[pos2];
//    }

//    cout<<"which case"<<"\n";
//    cout<<"pos pos2 "<<pos<<" "<<pos2<<" \n";   
//    cout<<"tetras t1 t2 "<<t1<<" "<<t2<<" \n";
//    cout<<"c1,c2,f1,f2 " <<c1<<" "<<c2<<" "<<f1<<" "<<f2<<" \n";
//    cout<<"accepted rho^2 prob "<<densitysquare<<" "<<prob<<" \n";
//    cout <<"\n";  

    //Now we know the charges. Should we add some kind of thermal update here as well?
    //now we copy the "update" part from the single spin sweep case.
    if((c1==0)&&(c2==0))
    {
        if(densitysquare>prob)
        {//flip if probability (here is just density) is larger than the random number
            config[pos]=-config[pos];
            config[pos2]=-config[pos2];
            //cout << "accepted??? WTF"<<"\n";
            //The proposal is accepted! we now keep track the position of charges created!
            //Now we determine where are the charges and what values are them
            if(flag==0)
            {//even sublattice
                if(f1==2){
                    //now we push t1 into the vectors
                    chargepairs.pposeven.push_back(t1);
                    chargepairs.nposeven.push_back(t2);
                }
                else{
                    chargepairs.pposeven.push_back(t2);
                    chargepairs.nposeven.push_back(t1);
                }
            }
            else if(flag==1)
            {//odd sublattice
                if(f1==2){
                    //we push t1 and t2 to odd stacks
                    chargepairs.pposodd.push_back(t1);
                    chargepairs.nposodd.push_back(t2);
                }
                else{
                    chargepairs.pposodd.push_back(t2);
                    chargepairs.nposodd.push_back(t1);
                }
            }
        }
    }
    if(((c1==2)&&(c2==-2))||((c1==-2)&&(c2==2))) // S+S- annihilates oposite charges on two tetrahaedra on the same sublatttice
    {
        if((f1==0)&&(f2==0))
        {
            config[pos]=-config[pos];
            config[pos2]=-config[pos2];
            //we now need to delete the two charges from the respective stacks!
            if(flag==0)
            {
                //delete the two charges on the even sublattice
                if(c1==2)
                {//we need to search and delete. This generically cost order length of the vectors.
                    chargepairs.pposeven.erase(std::remove(chargepairs.pposeven.begin(),chargepairs.pposeven.end(),t1),chargepairs.pposeven.end());
                    chargepairs.nposeven.erase(std::remove(chargepairs.nposeven.begin(),chargepairs.nposeven.end(),t2),chargepairs.nposeven.end());
                }
                else{
                    chargepairs.pposeven.erase(std::remove(chargepairs.pposeven.begin(),chargepairs.pposeven.end(),t2),chargepairs.pposeven.end());
                    chargepairs.nposeven.erase(std::remove(chargepairs.nposeven.begin(),chargepairs.nposeven.end(),t1),chargepairs.nposeven.end());
                }
            }
            else{
                if(c1==2){
                    chargepairs.pposodd.erase(std::remove(chargepairs.pposodd.begin(),chargepairs.pposodd.end(),t1),chargepairs.pposodd.end());
                    chargepairs.nposodd.erase(std::remove(chargepairs.nposodd.begin(),chargepairs.nposodd.end(),t2),chargepairs.nposodd.end());
                }
                else{
                    chargepairs.pposodd.erase(std::remove(chargepairs.pposodd.begin(),chargepairs.pposodd.end(),t2),chargepairs.pposodd.end());
                    chargepairs.nposodd.erase(std::remove(chargepairs.nposodd.begin(),chargepairs.nposodd.end(),t1),chargepairs.nposodd.end());
                }
            }
        }
    }
    if(((c1==2)&&(c2==0))||((c1==-2)&&(c2==0))) // S+S- moves charges from tetrahedon 1 to 2 and viceversa
    {
        if(((f1==0)&&(f2==2))||((f1==0)&&(f2==-2)))
        {
            config[pos]=-config[pos];
            config[pos2]=-config[pos2];
        }
    }
    if(((c1==0)&&(c2==2))||((c1==0)&&(c2==-2))) // S+S- moves charges from tetrahedon 1 to 2 and viceversa
    {
        if(((f1==2)&&(f2==0))||((f1==-2)&&(f2==0)))
        {
            config[pos]=-config[pos];
            config[pos2]=-config[pos2];
        }
    }
    return;
}



void jastrow(int &ntetra, int spinonc[],double table[][2],double jast[])
{

int n,i;

for(n=0;n<ntetra;n++)
{
   table[n][0]=0;
   table[n][1]=0;
   
   if(n%2==0)
   { 
      for(i=0;i<ntetra;i+=2)
      {
        if((i!=n)&&(spinonc[i]==1))
        {
          //table[n][0]=table[n][0]+jast[i][n]*spinonc[i];
          //cout<<" n i "<<n<<" "<<i<<"\n"; 
          table[n][0]=table[n][0]+jast[i+ntetra*n]*spinonc[i];  
        } 
        if((i!=n)&&(spinonc[i]==-1))
        {
          //table[n][1]=table[n][1]+jast[i][n]*spinonc[i];
          // cout<<" n i "<<n<<" "<<i<<"\n"; 
          table[n][1]=table[n][1]+jast[i+ntetra*n]*spinonc[i];
        }  
         
      } 
   
   }
   else if(n%2==1)
   {
      for(i=1;i<ntetra;i+=2)
      {
        if((i!=n)&&(spinonc[i]==1))
        {
           //cout<<" n i "<<n<<" "<<i<<"\n"; 
          table[n][0]=table[n][0]+jast[i+ntetra*n]*spinonc[i];
        }
        if((i!=n)&&(spinonc[i]==-1))
        {
           //cout<<" n i "<<n<<" "<<i<<"\n";
          table[n][1]=table[n][1]+jast[i+ntetra*n]*spinonc[i];
        }
      } 
   } 
}

return;
}

void updatejas(double table[][2],double jast[],int &t1,int &c1,int &f1,int &ntetra)
{
int n;


if(c1/2==0&&f1/2==1)
{
  if(t1%2==0)
  {
     for(n=0;n<ntetra;n+=2)
      {
       if(n!=t1)
       {  
        table[n][0]=table[n][0]+jast[t1+ntetra*n];
       } 
      }
  }
  else if(t1%2==1)
  {
     for(n=1;n<ntetra;n+=2)
      {
       if(n!=t1)
       {
        table[n][0]=table[n][0]+jast[t1+ntetra*n];
       }
      }  
  } 
}

if(c1/2==0&&f1/2==-1)
{
  if(t1%2==0)
  {
     for(n=0;n<ntetra;n+=2)
      {
       if(n!=t1)
       {
        table[n][1]=table[n][1]-jast[t1+ntetra*n];
       }
      }
  }
  else if(t1%2==1)
  {
     for(n=1;n<ntetra;n+=2)
      {
       if(n!=t1)
       {
        table[n][1]=table[n][1]-jast[t1+ntetra*n];
       }
      }
  }
}

if(c1/2==1&&f1/2==0)
{
  if(t1%2==0)
  {
     for(n=0;n<ntetra;n+=2)
      {
       if(n!=t1)
       {
        table[n][0]=table[n][0]-jast[t1+ntetra*n];
       }
      }
  }
  else if(t1%2==1)
  {
     for(n=1;n<ntetra;n+=2)
      {
       if(n!=t1)
       {
        table[n][0]=table[n][0]-jast[t1+ntetra*n];
       }
      }
  }
}

if(c1/2==-1&&f1/2==0)
{
  if(t1%2==0)
  {
     for(n=0;n<ntetra;n+=2)
      {
       if(n!=t1)
       {
        table[n][1]=table[n][1]+jast[t1+ntetra*n];
       }
      }
  }
  else if(t1%2==1)
  {
     for(n=1;n<ntetra;n+=2)
      {
       if(n!=t1)
       {
        table[n][1]=table[n][1]+jast[t1+ntetra*n];
       }
      }
  }
}

return;
}


double jastrowfactor(double jast[],int spinonc[],int &ntetra)
{

int i,j;
double par,impar;

par=0;
impar=0;
for (i=0;i<ntetra;i+=2)
{
 for (j=i+2;j<ntetra;j+=2)
  if(spinonc[i]*spinonc[j]==-1)
  {
    par=par-jast[i+ntetra*j];
  }
}


for (i=1;i<ntetra;i+=2)
{
 for (j=i+2;j<ntetra;j+=2)
  if(spinonc[i]*spinonc[j]==-1)
  {
    impar=impar-jast[i+ntetra*j];
  }
}

par=exp(par+impar);

/*for (i=0;i<ntetra;i++)
{
 cout<<"tetra spinonc "<<i<<" "<<spinonc[i]<<"\n";  
}

cout<<"par "<<par<<" \n";

cout<<"manual 110 136 "<<exp(jast[110+ntetra*136]*spinonc[110]*spinonc[136])<<" \n";

cout<<"jastrow[110 136]"<<jast[110+ntetra*136]<<" "<<jast[136+ntetra*110]<<" \n";
*/
return par;

}
void pair_flip2(int *config,int spinonc[],int ivic[][6],int tetra[][4],int connect[][2],int &L,double &densitysquare,double &density,int&pos,int &pos2,double &prob,int &ntetra,double table[][2], double jast[])
{
   
    int l,nhh;
    double rat,jastr,jasi,jasf;
    double tablec[ntetra][2]; 
 
    pos2=ivic[pos][pos2];
    //cout<<"spins"<<"\n";
    //cout<<"positions "<<pos<<" "<<pos2<<"\n";
    //cout<<"conf "<<config[pos]<<" "<<config[pos2]<<"\n"; 
    if(config[pos]*config[pos2]==1)
    {
       // cout<<"nospin flip"<<"\n";
        return;
    }
    int t1,t2;
    t1=connect[pos][0];
    t2=connect[pos][1];
    //repeated tetrahedra are deleted.
    if(t1==connect[pos2][0])
    {
        //In this case, the pair of spins connect two tetrahedron on "odd sublattice"
        t1=connect[pos2][1];
    }
    else if(t2==connect[pos2][1])
    {
        //the pair of spins connect two tetrahedron on "even sublattice"
        t2=connect[pos2][0];
    }
    //
    int c1,c2,f1,f2;
    c1=qcharge(t1,tetra,config);
    c2=qcharge(t2,tetra,config);
    
    config[pos]=-config[pos];
    config[pos2]=-config[pos2];
    f1=qcharge(t1,tetra,config);
    f2=qcharge(t2,tetra,config);
    config[pos]=-config[pos];
    config[pos2]=-config[pos2];
    
    if(abs(f1)==4||abs(f2)==4)
    {
         // wave function for Q =\pm 2  charges is zero  	 
     	 return;

    }  
  
     //cout<<"configuration x \n";
     //for(l=0;l<ntetra;l++)
     //{
     //cout<<l<<" "<<spinonc[l]<<" table "<<table[l][0]<<" "<<table[l][1]<<"\n"; 
        
     //}  
     
    //cout<<"position changes t1="<<t1<<" t2="<<t2<<"\n";

    //cout<<"charges c1 c2 f1 f2 "<<c1/2<<" "<<c2/2<<" "<<f1/2<<" "<<f2/2<<" \n"; 
    //cout<<"prob "<<prob<<" \n"; 

    

    rat=pow(density,abs(f1/2)+abs(f2/2))/pow(density,abs(c1/2)+abs(c2/2));

    jastr=1.0; 
    if(-(f1/2+c1/2)==-1)
    { 
    jastr=jastr*exp(table[t1][1]*(f1/2-c1/2));
    } 
    else if(-(f1/2+c1/2)==1)
    {
    jastr=jastr*exp(table[t1][0]*(f1/2-c1/2));
    }
   
    if(-(f2/2+c2/2)==-1)
    {
    jastr=jastr*exp(table[t2][1]*(f2/2-c2/2));
    }  
    else if(-(f2/2+c2/2)==1)
    {
    jastr=jastr*exp(table[t2][0]*(f2/2-c2/2));
    }

    //n=1 m=2
      //cout<<"t1 t2 "<<t1<<" "<<t2<<"\n";  
      //cout<<"f1 f2 "<< f1 << " "<<f2 <<"\n";
     //  cout<< "f1*f2 " << jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))<<"\n";
      // cout<<"c1 c2 "<< c1 << " "<<c2 <<"\n";
     //  cout<< "c1*c2 " << jast[t1+ntetra*t2]*double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2))<<"\n";

     //cout<<"should be 1 "<<jastr<<" \n";  
     jastr=jastr*exp(   
                        +jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))
                       
                        +jast[t1+ntetra*t2]*double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2)));   
             //cout<<exp(   -jast[t1+ntetra*t2]*double(c1/2)*double(f2/2)*0.5*(1.0-double(c1/2)*double(f2/2))
             //           +jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))
             //           -jast[t1+ntetra*t2]*double(c2/2)*double(f1/2)*0.5*(1.0-double(c2/2)*double(f1/2))
             //           +jast[t1+ntetra*t2]*double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2)))<<"\n";
    
     //cout<<"jastrow= "<<jastr<<"manual"<< exp( jast[t1+t2*ntetra]*f1*f2/4 )<<"\n";
     //cout<<"initial rat is "<<rat<<"\n"; 
     rat=rat*pow(jastr,2.0); 
     //cout<<"rat is "<<rat<<" jastr is "<<jastr<<" prob is "<<prob<<"\n";
     //jasi=jastrowfactor(jast,spinonc,ntetra);
     //spinonc[t1]=f1/2;
     //spinonc[t2]=f2/2;
     //jasf=jastrowfactor(jast,spinonc,ntetra);  
     //spinonc[t1]=c1/2;
     //spinonc[t2]=c2/2; 
     //rat=rat*pow(jasf/jasi,2.0);
     //cout<<"jastrow[t1 t2]"<<jast[t1+ntetra*t2]<<" "<<jast[t2+ntetra*t1]<<" \n";  
     //cout<<"comparison jasf/jasi vs ratio"<< jasf/jasi<<" "<<jastr <<" \n";

     //if(abs(jasf/jasi-jastr)>0.00001) 
     //{
     // cout<<"eeeeeeeeeeeeeeeeeeeeeeeee"<< abs(jasf/jasi-jastr)<< " \n";
         
     // exit(0);
     //}  
    //cout<<"rat= "<<rat<<" prob "<<prob<<"\n";
     //cout<<"ratio "<<rat<<" \n"; 
    if(rat<1.0)
    {
      if(prob<rat)
      {
      	config[pos]=-config[pos];
        config[pos2]=-config[pos2];
        spinonc[t1]=f1/2;
        spinonc[t2]=f2/2;
        //cout<<"accepted"<<"\n";
        updatejas(table,jast,t1,c1,f1,ntetra);
        updatejas(table,jast,t2,c2,f2,ntetra);
         //jastrow(ntetra,spinonc, table, jast);



     //for(l=0;l<ntetra;l++)
    //{
     //cout<<l<<" "<<spinonc[l]<<" table "<<table[l][0]<<" "<<table[l][1]<<"\n";
     //cout<<"V04 "<<jast[0+ntetra*4]<<" V02 "<<jast[0+ntetra*2]<<"\n";

     //} 
           
        
      }	
    	
    }
    else
    {
     	config[pos]=-config[pos];
        config[pos2]=-config[pos2];
        spinonc[t1]=f1/2;
        spinonc[t2]=f2/2;
        updatejas(table,jast,t1,c1,f1,ntetra);
        updatejas(table,jast,t2,c2,f2,ntetra);
        //jastrow(ntetra,spinonc, table, jast); 
    } 
   
    return;
}
//include the new Jastrow factor and new way of implement the density
void pair_flip3(int *config,int spinonc[],int ivic[][6],int tetra[][4],int connect[][2],int &L,int&pos,int &pos2,double &prob,int &ntetra,double table[][2], double jast[])
{
   
    int l,nhh;
    double rat,jastr,jasi,jasf;
    double tablec[ntetra][2]; 
 
    pos2=ivic[pos][pos2];
    //cout<<"spins"<<"\n";
    //cout<<"positions "<<pos<<" "<<pos2<<"\n";
    //cout<<"conf "<<config[pos]<<" "<<config[pos2]<<"\n"; 
    if(config[pos]*config[pos2]==1)
    {
       // cout<<"nospin flip"<<"\n";
        return;
    }
    int t1,t2;
    t1=connect[pos][0];
    t2=connect[pos][1];
    //repeated tetrahedra are deleted.
    if(t1==connect[pos2][0])
    {
        //In this case, the pair of spins connect two tetrahedron on "odd sublattice"
        t1=connect[pos2][1];
    }
    else if(t2==connect[pos2][1])
    {
        //the pair of spins connect two tetrahedron on "even sublattice"
        t2=connect[pos2][0];
    }
    //
    int c1,c2,f1,f2;
    c1=qcharge(t1,tetra,config);
    c2=qcharge(t2,tetra,config);
    
    config[pos]=-config[pos];
    config[pos2]=-config[pos2];
    f1=qcharge(t1,tetra,config);
    f2=qcharge(t2,tetra,config);
    config[pos]=-config[pos];
    config[pos2]=-config[pos2];
    
    if(abs(f1)==4||abs(f2)==4)
    {
         // wave function for Q =\pm 2  charges is zero  	 
     	 return;

    }  
  
     //cout<<"configuration x \n";
     //for(l=0;l<ntetra;l++)
     //{
     //cout<<l<<" "<<spinonc[l]<<" table "<<table[l][0]<<" "<<table[l][1]<<"\n"; 
        
     //}  
     
    //cout<<"position changes t1="<<t1<<" t2="<<t2<<"\n";

    //cout<<"charges c1 c2 f1 f2 "<<c1/2<<" "<<c2/2<<" "<<f1/2<<" "<<f2/2<<" \n"; 
    //cout<<"prob "<<prob<<" \n"; 

    

    jastr=1.0; 
    if(-(f1/2+c1/2)==-1)
    { 
    jastr=jastr*exp(table[t1][1]*(f1/2-c1/2));
    } 
    else if(-(f1/2+c1/2)==1)
    {
    jastr=jastr*exp(table[t1][0]*(f1/2-c1/2));
    }
   
    if(-(f2/2+c2/2)==-1)
    {
    jastr=jastr*exp(table[t2][1]*(f2/2-c2/2));
    }  
    else if(-(f2/2+c2/2)==1)
    {
    jastr=jastr*exp(table[t2][0]*(f2/2-c2/2));
    }

    //n=1 m=2
      //cout<<"t1 t2 "<<t1<<" "<<t2<<"\n";  
      //cout<<"f1 f2 "<< f1 << " "<<f2 <<"\n";
     //  cout<< "f1*f2 " << jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))<<"\n";
      // cout<<"c1 c2 "<< c1 << " "<<c2 <<"\n";
     //  cout<< "c1*c2 " << jast[t1+ntetra*t2]*double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2))<<"\n";

     //cout<<"should be 1 "<<jastr<<" \n";  
     jastr=jastr*exp(   
                        +jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))
                       
                        +jast[t1+ntetra*t2]*double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2)));   
             //cout<<exp(   -jast[t1+ntetra*t2]*double(c1/2)*double(f2/2)*0.5*(1.0-double(c1/2)*double(f2/2))
             //           +jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))
             //           -jast[t1+ntetra*t2]*double(c2/2)*double(f1/2)*0.5*(1.0-double(c2/2)*double(f1/2))
             //           +jast[t1+ntetra*t2]*double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2)))<<"\n";
    
     //cout<<"jastrow= "<<jastr<<"manual"<< exp( jast[t1+t2*ntetra]*f1*f2/4 )<<"\n";
     //cout<<"initial rat is "<<rat<<"\n"; 
     rat=pow(jastr,2.0); 
     //cout<<"rat is "<<rat<<" jastr is "<<jastr<<" prob is "<<prob<<"\n";
     //jasi=jastrowfactor(jast,spinonc,ntetra);
     //spinonc[t1]=f1/2;
     //spinonc[t2]=f2/2;
     //jasf=jastrowfactor(jast,spinonc,ntetra);  
     //spinonc[t1]=c1/2;
     //spinonc[t2]=c2/2; 
     //rat=rat*pow(jasf/jasi,2.0);
     //cout<<"jastrow[t1 t2]"<<jast[t1+ntetra*t2]<<" "<<jast[t2+ntetra*t1]<<" \n";  
     //cout<<"comparison jasf/jasi vs ratio"<< jasf/jasi<<" "<<jastr <<" \n";

     //if(abs(jasf/jasi-jastr)>0.00001) 
     //{
     // cout<<"eeeeeeeeeeeeeeeeeeeeeeeee"<< abs(jasf/jasi-jastr)<< " \n";
         
     // exit(0);
     //}  
    //cout<<"rat= "<<rat<<" prob "<<prob<<"\n";
     //cout<<"ratio "<<rat<<" \n"; 
    if(rat<1.0)
    {
      if(prob<rat)
      {
      	config[pos]=-config[pos];
        config[pos2]=-config[pos2];
        spinonc[t1]=f1/2;
        spinonc[t2]=f2/2;
        //cout<<"accepted"<<"\n";
        updatejas(table,jast,t1,c1,f1,ntetra);
        updatejas(table,jast,t2,c2,f2,ntetra);
         //jastrow(ntetra,spinonc, table, jast);



     //for(l=0;l<ntetra;l++)
    //{
     //cout<<l<<" "<<spinonc[l]<<" table "<<table[l][0]<<" "<<table[l][1]<<"\n";
     //cout<<"V04 "<<jast[0+ntetra*4]<<" V02 "<<jast[0+ntetra*2]<<"\n";

     //} 
           
        
      }	
    	
    }
    else
    {
     	config[pos]=-config[pos];
        config[pos2]=-config[pos2];
        spinonc[t1]=f1/2;
        spinonc[t2]=f2/2;
        updatejas(table,jast,t1,c1,f1,ntetra);
        updatejas(table,jast,t2,c2,f2,ntetra);
        //jastrow(ntetra,spinonc, table, jast); 
    } 
   
    return;
}

void singlespin_sweep_new(int *config,int ivic[][6],int tetra[][4],int connect[][2],int &L,double &densitysquare,int&flag,MTRand *myrand)
{
    int count,total=16*(int)pow(L,3),pos;
    double prob;
    //update the system once.
    for(count=0;count<total;count++)
    {
        pos=myrand->randInt(total-1);
        prob=myrand->rand();
        //cout<< pos <<"  " << "  "<< prob <<"  \n";
        singlespin_update_new(config,tetra,connect,L,pos,prob,densitysquare,flag);
    }
    return;
}

//measuring the total energy of the system
void e0total(int *config,int tetra[][4],int &ntetra,int &L,double &estep)
{
    //we only go through the zero sublattice. 
    int x,y,z,t,pyro;
    int ztemp=8*pow(L,2);
    int ytemp=8*L;
    int temp;
    estep=0;

    for(z=0;z<ntetra;z++)
    {
     estep+=pow((double)qcharge(z,tetra,config),2.0)/8.0;
     //cout<<"z="<<z; 
    }
    return;
}  
    //cout<<"energylocal tetras"<<estep<<"\n";
//    estep=0;  
//    for(z=0;z<L;z++)
//    {
//        for(y=0;y<L;y++)
//        {
//            for(x=0;x<L;x++)
//            {
//                for(pyro=0;pyro<4;pyro++)
//                {
//                    ///the last part should multiply by 2. We only have 0, 1,2,3
//                    t=z*ztemp+y*ytemp+x*8+pyro*2;
//                    temp=qcharge(t,tetra,config);
//                     //cout<<"t="<<t;     
//                    estep+=((double)pow(temp,2.0))/8.0;
//                    t+=1;
//                    
//                    temp=qcharge(t,tetra,config);
//                    // cout<<"t="<<t;
//                    estep+=((double)pow(temp,2.0))/8.0;
//                    
//                }
//            }
//        }
//    }
//    cout<<"energylocal?"<<estep<<"\n";
//    return;
//}

//Now we need the routines to do measurement
//
//the step for ground state energy. We probably need more stuff to determine step.
double e0_step(int *config,int &c1,int&c2,int&f1,int&f2,double&density)
{
    double step=0;
    if((c1==0)&&(c2==0))
    {
        step=-density;
    }
    else if((c1==0)&&(c2==2)&&(f1==2)&&(f2==0))
    {
        step=-1;
    }
    else if((c1==0)&&(c2==-2)&&(f1==-2)&&(f2==0))
    {
        step=-1;
    }
    else if((c1==2)&&(c2==0)&&(f1==0)&&(f2==2))
    {
        step=-1;
    }
    else if((c1==-2)&&(c2==0)&&(f1==0)&&(f2==-2))
    {
        step=-1;
    }
    else if((c1==2)&&(c2==-2)&&(f1==0)&&(f2==0))
    {
        step=-1/density;
    }
    else if((c1==-2)&&(c2==2)&&(f1==0)&&(f2==0))
    {
        step=-1/density;
    }
    return step;
}




//auxilliary inline function

/*Juan's routines starts..... Here
 -----------------------------------------------------*/
int searchit (int site, int ntetra, int z, int tetra[][4])
{
    int tetraout;
    int zz,ii;
    
    if (z%2==0)
    {
        for (zz=1;zz<ntetra;zz=zz+2)
        {
            //cout<<"zz odd"<<zz<<"\n";
            for (ii=0;ii<4;ii++)
            {
                if ( tetra[zz][ii]==site )
                {
                    tetraout=zz;
                    //cout<<"success site"<<site<<"tetras "<<z<<" "<<tetraout<<"\n";
                    return tetraout;
                }
            }
            
        }
        
    }
    else if (z%2==1)
    {
        for (zz=0;zz<ntetra;zz=zz+2)
        {
            //cout<<"zz even"<<zz<<"\n";
            for (ii=0;ii<4;ii++)
            {
                if ( tetra[zz][ii]==site )
                {
                    tetraout=zz;
                    //cout<<"success site"<<site<<"tetras "<<z<<" "<<tetraout<<"\n";
                    return tetraout;
                }
            }
            
        }
        
    }
    return tetraout;
}


void boundary(int& xt, int& yt,int& zt, int& L,int& loc)
{
    int zp,yp,xp;
    zp=L*L*16;
    yp=L*16;
    xp=16;
    if (xt<0)xt=xt+L;
    if (xt>L-1)xt=xt-L;
    if (yt<0)yt=yt+L;
    if (yt>L-1)yt=yt-L;
    if (zt<0)zt=zt+L;
    if (zt>L-1)zt=zt-L;
    loc=zt*zp+yt*yp+xt*xp;
}
//ivic,tetra,connect,L,nh,ntetra
void latt(int ivic[][6],int tetra[][4],int connect[][2], int &L, int &nh, int &ntetra)
{
    int x,y,z,pyro,xt,yt,zt,zp,yp,xp,loc;
    int site,counter,tetracount,tetraout;
    
    zp=L*L*16;
    yp=L*16;
    xp=16;
    
    counter=0;
    tetracount=0;
    for (z=0;z<L;z++)
    {
        for (y=0;y<L;y++)
        {
            for (x=0;x<L;x++)
            {
                for (pyro=0;pyro<4;pyro++)
                {
                    for (site=0;site<4;site++)
                    {
                        
                        if (site==0)
                        {
                            // neighbors inside the same unit cell x,y,z
                            ivic[counter][0]=counter+1;
                            ivic[counter][1]=counter+2;
                            ivic[counter][2]=counter+3;
                            
                            if (pyro==0)
                            {
                                //x,y-1,z-1
                                xt=x;yt=y-1;zt=z-1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+5;  //1st pyro
                                //cout<<loc<<"\n";
                                
                                //x-1,y,z-1
                                xt=x-1;yt=y;zt=z-1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+10; // 2nd pyro
                                // cout<<loc<<"\n";
                                
                                //x-1,y-1,z
                                xt=x-1;yt=y-1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+15; //3th pyro
                                //cout<<loc<<"\n";
                            }
                            else if (pyro==1)
                            {
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+1; //0th pyro
                                
                                //x-1,y,z
                                xt=x-1;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+14; //3th pyro
                                
                                //x-1,y,z
                                xt=x-1;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+11; //2nd pyro
                                //cout<<loc<<"\n";
                                //return;
                            }
                            else if (pyro==2)
                            {
                                //x,y-1,z
                                xt=x;yt=y-1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+13; //3th pyro
                                
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+2; //0th pyro
                                
                                //x,y-1,z
                                xt=x;yt=y-1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+7; //1st pyro
                                
                            }
                            else if (pyro==3)
                            {
                                //x,y,z-1
                                xt=x;yt=y;zt=z-1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+9; //2th pyro
                                
                                //x,y,z-1
                                xt=x;yt=y;zt=z-1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+6; //1th pyro
                                
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+3; //0th pyro
                                
                            }
                            
                            
                        }
                        else if (site==1)
                        {
                            
                            // neighbors inside the same unit cell x,y,z
                            xt=x;yt=y;zt=z;
                            boundary(xt,yt,zt,L,loc);
                            ivic[counter][0]=counter-1;
                            ivic[counter][1]=counter+1;
                            ivic[counter][2]=counter+2;
                            
                            if (pyro==0)
                            {
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+4; //1st
                                
                                //x-1,y,z
                                xt=x-1;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+14; //3rd pyro
                                
                                //x-1,y,z
                                xt=x-1;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+11; //2nd pyro
                            }
                            else if (pyro==1)
                            {
                                //x,y+1,z+1
                                xt=x;yt=y+1;zt=z+1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+0; //0th
                                
                                //x-1,y+1,z
                                xt=x-1;yt=y+1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+10; // 2nd
                                
                                //x-1,y,z+1
                                xt=x-1;yt=y;zt=z+1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+15; //  3rd
                                
                            }
                            else if (pyro==2)
                            {
                                //x,y,z+1
                                xt=x;yt=y;zt=z+1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+12; //3rd
                                
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+6; //1st
                                
                                //x,y,z+1
                                xt=x;yt=y;zt=z+1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+3; //0th
                                
                            }
                            else if (pyro==3)
                            {
                                //x,y+1,z
                                xt=x;yt=y+1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+8; // 2nd
                                
                                //x,y+1,z
                                xt=x;yt=y+1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+2; //
                                
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+7; //1st
                                
                            }
                            
                        }
                        else if (site==2)
                        {
                            
                            // neighbors inside the same unit cell x,y,z
                            xt=x;yt=y;zt=z;
                            boundary(xt,yt,zt,L,loc);
                            ivic[counter][0]=counter-2;
                            ivic[counter][1]=counter-1;
                            ivic[counter][2]=counter+1;
                            
                            if (pyro==0)
                            {
                                //x,y-1,z
                                xt=x;yt=y-1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+7; //1st
                                
                                //x,y-1,z
                                xt=x;yt=y-1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+13; //3rd
                                
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+8; //2nd
                            }
                            else if (pyro==1)
                            {
                                //x,y,z+1
                                xt=x;yt=y;zt=z+1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+3; //0th
                                
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+9; //2nd
                                
                                //x,y,z+1
                                xt=x;yt=y;zt=z+1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+12; //3rd
                            }
                            else if (pyro==2)
                            {
                                //x,y-1,z+1
                                xt=x;yt=y-1;zt=z+1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+15; //3rd
                                
                                //x+1,y-1,z
                                xt=x+1;yt=y-1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+5; //1st
                                
                                //x+1,y,z+1
                                xt=x+1;yt=y;zt=z+1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+0; //0th
                                
                            }
                            else if (pyro==3)
                            {
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+11; //2nd
                                
                                //x+1,y,z
                                xt=x+1;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+1; //0th
                                
                                //x+1,y,z
                                xt=x+1;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+4; //1st
                                
                            }
                        }
                        else if (site==3)
                        {
                            // neighbors inside the same unit cell x,y,z
                            xt=x;yt=y;zt=z;
                            boundary(xt,yt,zt,L,loc);
                            ivic[counter][0]=counter-3;
                            ivic[counter][1]=counter-2;
                            ivic[counter][2]=counter-1;
                            
                            if (pyro==0)
                            {
                                //x,y,z-1
                                xt=x;yt=y;zt=z-1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+9; //2nd
                                
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+12; //3rd
                                
                                //x,y,z-1
                                xt=x;yt=y;zt=z-1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+6; //1st
                                
                            }
                            else if (pyro==1)
                            {
                                //13 8 2
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+13; //3rd
                                
                                //x,y+1,z
                                xt=x;yt=y+1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+8; // 2nd
                                
                                //x,y+1,z
                                xt=x;yt=y+1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+2; // 0th
                                
                            }
                            else if (pyro==2)
                            {  //1 4 14
                                //x+1,y,z
                                xt=x+1;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+1; //0th
                                
                                //x+1,y,z
                                xt=x+1;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+4; //1st
                                
                                //x,y,z
                                xt=x;yt=y;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+14; //3rd
                            }
                            else if (pyro==3)
                            { //5 0 10
                                //x+1,y,z-1
                                xt=x+1;yt=y;zt=z-1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][3]=loc+5; //1st
                                
                                //x+1,y+1,z
                                xt=x+1;yt=y+1;zt=z;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][4]=loc+0; //0th
                                
                                //x,y+1,z-1
                                xt=x;yt=y+1;zt=z-1;
                                boundary(xt,yt,zt,L,loc);
                                ivic[counter][5]=loc+10; // 2nd
                                
                            }
                        }
                        counter=counter+1;
                    }
                }
            }
        }
    }
    
    
    // list of tetrahedra
    tetracount=0;
    for (z=0;z<L;z++)
    {
        for (y=0;y<L;y++)
        {
            for (x=0;x<L;x++)
            {
                xt=x;yt=y;zt=z;
                boundary(xt,yt,zt,L,loc);
                
                for (pyro=0;pyro<4;pyro++)
                {
                    
                    tetra[tetracount][0]=loc+pyro*4+0;
                    tetra[tetracount][1]=loc+pyro*4+1;
                    tetra[tetracount][2]=loc+pyro*4+2;
                    tetra[tetracount][3]=loc+pyro*4+3;
                    tetracount=tetracount+1;
                    
                    tetra[tetracount][0]=loc+pyro*4+3;
                    tetra[tetracount][1]=ivic[loc+pyro*4+3][3];
                    tetra[tetracount][2]=ivic[loc+pyro*4+3][4];
                    tetra[tetracount][3]=ivic[loc+pyro*4+3][5];
                    tetracount=tetracount+1;
                    
                }
            }
        }
    }
    
    /*cout<<ntetra<<"\n"<<tetracount<<"\n";
    //cout<<"sites belongin to each tetrahedra z"<<"\n";
    for (z=0;z<ntetra;z++)
    {
        cout<<z<<"    "<<tetra[z][0]<<" "<<tetra[z][1]<<" "<<tetra[z][2]<<" "<<tetra[z][3]<<"\n";
    }
    cout<<"\n";*/
    
    
    //finding which tetrahedra each site connects to
    //initialize
    for (z=0;z<nh;z++)
    {
        for (x=0;x<2;x++)
        {
            connect[z][x]=-1;
        }
    }
    
    for (z=0;z<ntetra;z++)
    {
        for (x=0;x<4;x++)
        {
            site=tetra[z][x];
            tetraout=searchit(site,ntetra,z,tetra);
            //cout<<"success site"<<site<<"tetras "<<z<<" "<<tetraout<<" connsofar "<<connect[site][1]<<"\n";
            
            if (connect[site][0]==-1)
            {
                if (z%2==0)
                {
                    connect[site][0]=z;
                    connect[site][1]=tetraout;
                }
                else if (z%2==1)
                {
                    connect[site][0]=tetraout;
                    connect[site][1]=z;
                }
            }
            
        }
    }
    
    /*cout<<" site x connects which tetrahedra"<<"\n";
    for (z=0;z<nh;z++)
    {
        cout<<"x="<<z<<"   tetrahedra "<<connect[z][0]<<" "<<connect[z][1]<<"\n";
    }

    cout<<" table of nearest neighbors of site x "<<"\n";*/
    /*for (z=0;z<nh;z++)
    {
        cout<<"x="<<z<<" "<<ivic[z][0]<<" "<<ivic[z][1]<<" "<<ivic[z][2]<<" "<<ivic[z][3]<<" "<<ivic[z][4]<<" "<<ivic[z][5]<<"\n ";
    }*/

    
}

int defected(int & ctetra, int *config, int tetra[][4])
{
    int charge;
    int cc;
    
    charge=0;
    for (cc=0;cc<4;cc++)
    {
        charge=charge+config[tetra[ctetra][cc]];
    }
    
    //cout<<"charge"<<charge<<"tetra"<<ctetra;
    return charge;
}

void printv(int * config,int & nh)
{
    int i;
    cout<<"\n";
    for (i=0;i<nh;i++)
    {
        cout<<i<<"  "<<config[i]<<"\n";
    }
    cout<<"\n \n";
    
}
void inoutspin(int & ctetra,int *config, int tetra[][4],int spin[2])
{
    int coun,ar;
    
    spin[0]=0;
    spin[1]=0;
    
    if (ctetra%2==0)
    {
        coun=0;
        for (ar=0;ar<4;ar++)
        {
            if (config[tetra[ctetra][ar]]==-1)
            {
                spin[coun]=tetra[ctetra][ar];
                coun=coun+1;
            }
        }
    }
    else if (ctetra%2==1)
    {
        coun=0;
        for (ar=0;ar<4;ar++)
        {
            if (config[tetra[ctetra][ar]]==1)
            {
                spin[coun]=tetra[ctetra][ar];
                coun=coun+1;
            }
        }
        
    }
    return;
}

void collect(double&tempature, double&eclassical,double&esquare,double data[],double data2[],int&ndat ,int&nh, int&msteps, int&i)
{
    double err[ndat];
    double aver[ndat];
    int ii;
    data[0]=data[0]+eclassical/(msteps);
    data[1]=data[1]+esquare/msteps;
    data2[0]=data2[0]+pow(eclassical/(msteps),2.0);
    data2[1]=data2[1]+pow(esquare/msteps,2.0);
   
    for(ii=0;ii<ndat;ii++)
    {
     aver[ii]=data[ii]/((double)i+1.0);  
     err[ii]=sqrt( abs( pow(data[ii]/((double)i+1.0),2.0)-data2[ii]/((double)i+1.0)  )/( (double)i+1.0)); 
    }

    ofstream output;
    output.open("results.txt");
    output << "bins:          " << i+1  <<"\n";  
    output << "Energy:        " << aver[0]/(nh) <<" pm "<<err[0]/(nh)<< "\n";
    output << "Variance TotE: " << (aver[1]-pow(aver[0],2.0)) <<" pm "<<(err[1]+2.0*err[0]*aver[0])<< "\n"; 
    //output << "Specific: " << (aver[1]-pow(aver[0],2.0))/(pow(tempature,2.0))/(nh/16) <<" pm "<<(err[1]+2.0*err[0]*aver[0] )/(pow(tempature,2.0))/(nh/16)<< "\n"; 
    output.close();
    eclassical=0.0;
    esquare=0.0; 

}


void loopupdate(int *config, int ivic[][6],int tetra[][4],int connect[][2], int &L, int &nh, int &ntetra,int &visits,int &went, MTRand *myrand)
{
    int itetra,v,ctetra,ii,vstart,vfinal,tet,toflip;
    int go = 0,ar,outs,tcharge;
    int spin[2],coun,countervisits = 0;
    int visited[ntetra][3]; // 0: visited or not. 1: order in the line of visits. 2: which spin was visited.
    int ordered[ntetra];   //  which tetrahedra are visited in order
  
     
    //initializating variables tracking which tetrahedra are visited and in what order.
    tcharge=0;
    for (v=0;v<ntetra;v++)
    {
        tcharge=tcharge+defected(v,config,tetra);
        ordered[v]=-1;  
        for (ii=0;ii<3;ii++)
        {
            visited[v][ii]=-1;
        } 
    }
    
    //if(tcharge!=0) cout<<"total charge in configuration= "<<tcharge<<"\n";
 

    // select a random tetrahedra to start the loop   
    //itetra=rand() % ntetra; // old random number generator
     itetra=myrand->randInt(ntetra-1); 
   
    ctetra=itetra; // current tetrahedra ctetra
    //cout<<"random tetrahedron"<<itetra<<"\n";
    
    
    do {
        
        //cout<<"visited ctetra= "<<ctetra<<"visited="<< visited[ctetra][0]<<"  \n";
        if (visited[ctetra][0]==-1)
        {
            
            visited[ctetra][0]=1;
            ordered[countervisits]=ctetra; 
            visited[ctetra][1]=countervisits;
            countervisits=countervisits+1;
            
            //cout<<"defected?=    "<<defected(ctetra,config,tetra)<<"     "; 
            if (defected(ctetra,config,tetra)==0)
            {
                
                // if not defected then randomly choose outward spin and locate tetrahedron located through that spin   
                if (ctetra%2==0)
                {
                    // even tetrahedron: up (+1) means in; down (-1) means out.
                    inoutspin(ctetra,config,tetra,spin);
                    //outs=rand()%2; // old random number
                    outs=myrand->randInt(1);               
                    //cout<<"out spins"<<spin[0]<<spin[1]<<"chosen"<<spin[outs]<<"\n";
                    visited[ctetra][2]=spin[outs];
                    if (connect[spin[outs]][0]==ctetra)
                    { 
                        ctetra=connect[spin[outs]][1];
                    }
                    else if (connect[spin[outs]][1]==ctetra)
                    { 
                        ctetra=connect[spin[outs]][0];      
                    }
                    //cout << "next chosen tetrahedron "  <<ctetra<<"\n";
                    
                }
                else if (ctetra%2==1)
                {
                    // odd tetrahedron: up (+1) means out; down (-1) means in.
                    inoutspin(ctetra,config,tetra,spin);
                    //outs=rand()%2; 
                    outs=myrand->randInt(1); 
                    //cout<<"out spins"<<spin[0]<<spin[1]<<"chosen"<<spin[outs]<<"\n";   
                    visited[ctetra][2]=spin[outs];
                    if (connect[spin[outs]][0]==ctetra)
                    { 
                        ctetra=connect[spin[outs]][1];
                    } 
                    else if (connect[spin[outs]][1]==ctetra)
                    {
                        ctetra=connect[spin[outs]][0];   
                    }
                    //cout << "next chosen tetrahedron "  <<ctetra<<"\n"; 
                }  
                
            }
            else if (defected(ctetra,config,tetra)!=0)
            {
                // abort the loop
                go=-1; 
                break; 
            }  
            
        }
        else if (visited[ctetra][0]==1)
        {
            // loop found
            vstart=visited[ctetra][1]; 
            vfinal=countervisits-1; 
            go=1; 
        }
        
        //cout<<"  counting visits  "<<countervisits <<"  \n";
        
    } while (go==0);
    
    
    
    
    // flipping the spins along the loop
    if (go==1)
    {  
        //cout<<"THERE WAS A LOOP vstart"<<vstart<<"vfinal="<<vfinal<<"\n";
        for (ii=vstart;ii<=vfinal;ii++)  
        {
            
            tet=ordered[ii];
            toflip=visited[tet][2];
            config[toflip]=-config[toflip];
            
        }
    }
    
    tcharge=0;
    for (v=0;v<ntetra;v++)
    {
        tcharge=tcharge+defected(v,config,tetra);
    }
    
    //if(tcharge!=0)cout<<"total charge in final configuration= "<<tcharge<<"\n";

    went=go;
    visits=countervisits;
    
    return;
}

void collect2(double&tempature, double&eclassical,double&esquare,double data[],double data2[],int&ndat ,int&nh, int&msteps, int&i,std::ofstream &output_data,double &density,double &t_tilde)
{
    double err[ndat];
    double aver[ndat];
    int ii;
    data[0]=data[0]+eclassical/(msteps);
    data[1]=data[1]+esquare/msteps;
    data2[0]=data2[0]+pow(eclassical/(msteps),2.0);
    data2[1]=data2[1]+pow(esquare/msteps,2.0);
    
    for(ii=0;ii<ndat;ii++)
    {
        aver[ii]=data[ii]/((double)i+1.0);
        err[ii]=sqrt( abs( pow(data[ii]/((double)i+1.0),2.0)-data2[ii]/((double)i+1.0)  )/( (double)i+1.0));
    }
    
    //ofstream output;
    //output.open("results.txt");
    //output << "bins:          " << i+1  <<"\n";
    //output << "Energy:        " << aver[0]/(nh) <<" pm "<<err[0]/(nh)<< "\n";
    //output << "Variance TotE: " << (aver[1]-pow(aver[0],2.0)) <<" pm "<<(err[1]+2.0*err[0]*aver[0])<< "\n";
    //output << "Specific: " << (aver[1]-pow(aver[0],2.0))/(pow(tempature,2.0))/(nh/16) <<" pm "<<(err[1]+2.0*err[0]*aver[0] )/(pow(tempature,2.0))/(nh/16)<< "\n";
    
    //output.close();
    output_data<< density << "\t" << t_tilde << "\t";
    output_data<< i+1 << "\t" << aver[0]/(double)nh <<"\t"<<err[0]/(double)(nh)<< "\t";
    output_data<< (aver[1]-pow(aver[0],2.0))<<"\t"<<(err[1]+2.0*err[0]*aver[0])<< "\n";
    eclassical=0.0;
    esquare=0.0;
}
