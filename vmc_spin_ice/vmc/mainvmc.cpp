//
//  main.cpp
//  vmc_spin_ice
//
//  Created by Zhihao Hao on 2014-08-21.
//  Copyright (c) 2014 Zhihao Hao. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "MersenneTwister.h"
#include "routines.h"
#include "estimator.h"


int main()
{

    // insert code here...
    // system size L*L*L*4*4
    int L=2;
    cout << "give me L"<<"\n";
    cin >> L; 
    int nh;
    int ntetra;
    //density
    double rho;
    double seedin=20.0;
    //spin configuration and
    MTRand *myrand=new MTRand();
    myrand->seed(seedin);
    //seedin=myrand->rand();
    int *config;
    double *alpha;
    double cp;
    config=new int[L*L*L*4*4];
    //initialization!!!!!!!!!!!!!!!!!!!! TEST
    //Notice X only works with even L. 
    //initialize_spinice_X(config,L);
    initialize_spinice_q0(config,L);
    //
    nh=16*(int)pow(L,3);
    ntetra=8*(int)pow(L,3);
    int ivic[nh][6];  // each site its 6 neighbours
    int tetra[ntetra][4]; // each tetrahedron and their 4 sites
    int connect[nh][2]; // each site connects two tetrahedra
    //construction of tables.
    latt(ivic,tetra,connect,L,nh,ntetra);
    int pos,pos2;
    int c1,c2,f1,f2,t,flag;
    //Let's try the single spin update.
    int x,y,z,pyro,site,temp,count,x0,total=10;
    int stat=1000000;
    int zpro=16*pow(L,2);
    int ypro=16*L;
    int xpro=16;
    for(z=0;z<L;z++)
    {
        for(y=0;y<L;y++)
        {
            for(x=0;x<L;x++)
            {
                std::cout<<"(x,y,z) is:"<<x<<'\t'<<y<<'\t'<<z<<'\n';
                for(pyro=0;pyro<4;pyro++)
                {
                    for(site=0;site<4;site++)
                    {
                        std::cout<<config[z*zpro+y*ypro+x*xpro+pyro*4+site]<<'\n';
                    }
                }
            }
        }
    }
    e0total(config,tetra,ntetra,L,cp);
    cout <<"energy should be zero for proper initiali"<<cp<<" tetra"<<ntetra<<"\n";
 
// if(L==1)
// {   
//  config[0]=1;
//  config[1]=-1;
//  config[2]=-1;
//  config[3]=1;
//  config[4]=1;
//  config[5]=-1;
//  config[6]=-1;
//  config[7]=1;
//  config[8]=1;
//  config[9]=-1;
//  config[10]=-1;
//  config[11]=1;
//  config[12]=1;
//  config[13]=-1;
//  config[14]=-1;
// config[15]=1;
// }
    //some simple test of pair update.
    //pos=3;
    //t=connect[pos][0];
    //c1=qcharge(t,tetra,config);
    //f1=c1-2*config[pos];
    //t=connect[pos][1];
    //c2=qcharge(t,tetra,config);
    //f2=c2+2*config[pos];
    double prob,density,densitysquare,tempature,jp;
    //pos=3;
    //prob=0.1;
    jp=0.1; 
    density=0.25;
    cout << "give me density"<<"\n";
    cin >> density;
  
    densitysquare=pow(density,2.0);
    //flag=0;
    //test
    //int pos2=myrand->randInt(5);
    //pair_flip(config,ivic,tetra,connect,L,densitysquare,pos,pos2,prob);


    // VMC data
    int thbins; // number of bins during thermalization
    int nbins; // number of bins production run
    int msteps; // bin  length  
    int nloops; // number of loops per monte carlo step
    int k,i,lo; // counters
    int cloop,went,avvisit,visits; //
    double ave;
    int ndat=2;
    double data[ndat],data2[ndat];
    double esquare=0,eclassical=0,estep;
    
    thbins=20;
    nbins=500;
    msteps=2000;
    nloops=100;
    
    
    tempature=0.0;     
    flag=0;
    e0total(config,tetra,ntetra,L,estep); 
    cout <<"charge^2"<<estep<<"\n"; 
      
    // Thermalization
    seedin=myrand->randInt();
    for(i=0;i<thbins;i++)
    {
      cout<<i<<" "<<nloops <<"\n";
      for(k=0;k<msteps;k++)
      {
         // pair spin flip
         for(lo=0;lo<nh;lo++) 
         {
         pos=myrand->randInt(nh-1); // random spin  to flip     
         pos2=myrand->randInt(5);   // random neighbor  to flip
         prob=myrand->rand();
         //cout<<"pos pos2 nh prob "<<pos<<" "<<pos2<<" "<<nh<<" "<<prob<<" \n"; 
         //pair_flip(config,ivic,tetra,connect,L,densitysquare,pos,pos2,prob);    
         pair_flip2(config,ivic,tetra,connect,L,densitysquare,density,pos,pos2,prob);

         //e0total(config,tetra,ntetra,L,estep);
         //cout<<"diagonal energy"<<estep<<"\n";
         //cout << "                   "<<" \n"; 
         }
         
        // singlespin_sweep_new(config,ivic,tetra,connect,L,tempature,flag,myrand);
  
         //loop update
        
         avvisit=0; 
         cloop=0;
         for(lo=0;lo<nloops;lo++)
         {
           loopupdate(config,ivic,tetra,connect,L,nh,ntetra,visits,went,myrand);
           if(went==1)
           {        
             avvisit=avvisit+visits;
             cloop=cloop+1;   
           }
         }
        //e0total(config,tetra,ntetra,L,estep);
        //cout<<"diagonal energy"<<estep<<"\n";
        //cout << "                   "<<" \n";
         
         //adjusting the required number of loops 
         if(cloop>0)
         { 
           ave=(double)avvisit/(double)cloop;
           //cout<<i<<" "<<"avvisit "<<ave <<" nloops "<<nloops << "crit "<<nh/(2*(int)ave) <<"\n";
           if(nloops<nh/(2*(int)ave))
           {
              nloops=nloops+1; //*abs(nloops-nh/(2*(int)ave));
           }   
           else
           {
              nloops=nloops-1; //*abs(nloops-nh/(2*(int)ave)); 
              if(nloops<1)nloops=1;
           } 
         } 
          
      }  

    }
    cout<<i<<" "<<"avvisit "<<ave <<" nloops "<<nloops << "crit "<<nh/(2*(int)ave) <<"\n";

    // initialize measurements
    eclassical=0.0;
    estep=0.0;
    for(i=0;i<ndat;i++)
    {
     data[i]=0.0;
     data2[i]=0.0;
    }   
 
    for(i=0;i<nbins;i++) 
    {
      for(k=0;k<msteps;k++)
      {
         
         // pair spin flip
         for(lo=0;lo<nh;lo++)
         {
         pos=myrand->randInt(nh-1); // random spin  to flip
         pos2=myrand->randInt(5);   // random neighbor  to flip
         prob=myrand->rand();
         //pair_flip(config,ivic,tetra,connect,L,densitysquare,pos,pos2,prob);
         pair_flip2(config,ivic,tetra,connect,L,densitysquare,density,pos,pos2,prob);
         }
         //singlespin_sweep_new(config,ivic,tetra,connect,L,tempature,flag,myrand); 
         // loopupdate
         for(lo=0;lo<nloops;lo++)
         {
          loopupdate(config,ivic,tetra,connect,L,nh,ntetra,visits,went,myrand);
         }
        
         // measurements
         //e0total(config,tetra,ntetra,L,estep);
//         energy_est2(int *config,int tetra[][4],int &ntetra,int ivic[][6],int connect[][2],int &L,double &density,double &jp,double &estep)  
         //energy_est(config,tetra,connect,L,density,jp,estep);
         energy_est2(config,tetra,ntetra,ivic,connect,L,nh,density,jp,estep);
         eclassical+=estep;
         estep=pow(estep,2.0);
         esquare+=estep;   
      }
      collect(tempature,eclassical,esquare,data,data2,ndat,nh,msteps,i);
       //e0total(config,tetra,ntetra,L,estep);
       //cout<<"diagonal energy"<<estep<<"\n";
      //cout <<"eclassical zero?"<< eclassical <<"\n"; 
    }
 
    return 0;
}

