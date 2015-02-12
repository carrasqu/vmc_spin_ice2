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
#include <stdlib.h>
#include <algorithm>
#include "MersenneTwister.h"
#include "routines.h"
#include "estimator.h"


int main()
{

    // insert code here...
    // system size L*L*L*4*4
    int L=4;
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
    initialize_spinice_X(config,L);
    //
    nh=16*(int)pow(L,3);
    ntetra=8*(int)pow(L,3);
    int ivic[nh][6];  // each site its 6 neighbours
    int tetra[ntetra][4]; // each tetrahedron and their 4 sites
    int connect[nh][2]; // each site connects two tetrahedra
    double disvec[4][3]={{0.0,0.0,0.0},{0.0,0.5,0.5},{0.5,0.0,0.5},{0.5,0.5,0.0}};
    //location to determine the pair of charges.
    charge_pair chargepairs;
    //construction of tables.
    latt(ivic,tetra,connect,L,nh,ntetra);
    int pos;
    int c1,c2,f1,f2,t,flag;
    //Let's try the single spin update.
    int x,y,z,pyro,site,temp,count,x0,total=10;
    int stat=1000000;
    int zpro=16*pow(L,2);
    int ypro=16*L;
    //int xpro=16;
    //Now we are trying to test the energy estimator.
    //some simple test of pair update.
    /* we now test the new routine we set up*/
    int llgg=pow(L,3);
    double t_tilde=0.01,eta=1.0;
    double correl[llgg][16];
    double d2=0.016;  
    //now we call the function
    spinon_correlation(correl, t_tilde, eta, L,d2);
 
    zpro=pow(L,2);
    ypro=L;
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
                        std::cout<<x<<"\t"<<y<<"\t"<<z<<" true (x,y,z) is:"<<(double)x+disvec[pyro][0]-disvec[site][0]<<'\t'<<(double)y+disvec[pyro][1]-disvec[site][1]<<'\t'<<(double)z+disvec[pyro][2]-disvec[site][2]<<'\n';
                        std::cout<<"mu_1: "<<pyro<<", "<<"mu_2: "<<site<<"\t"<<correl[z*zpro+y*ypro+x][pyro*4+site]<<'\n';
                        
                    }
                }
            }
        }
    }
    

    
  int x1,y1,z1,x2,y2,z2,d,loc1,loc2,tlab1,tlab2,tindex1,tindex2,flag_here;
  double jast[ntetra][ntetra]; 
 
    for (z1=0;z1<L;z1++)
    {
        for (y1=0;y1<L;y1++)
        {
            for (x1=0;x1<L;x1++)
            {
              loc1=x1+y1*L+z1*L*L; 

                for (z2=0;z2<L;z2++)
                {
                     for (y2=0;y2<L;y2++)
                     {
                         for (x2=0;x2<L;x2++)
                         {
                            
                             loc2=x2+y2*L+z2*L*L;
                             d=correl_index(loc1,loc2,L,flag_here);
                             for(tlab1=0;tlab1<4;tlab1++)
                             {
                                 for(tlab2=0;tlab2<4;tlab2++)
                                 {
                                   tindex1=z1*L*L*8+y1*L*8+x1*8+2*tlab1;
                                   tindex2=z2*L*L*8+y2*L*8+x2*8+2*tlab2;
                                   // std::cout<<"tetras considered "<<tindex1<<" "<<tindex2<<"\n" ;
                                   if(flag_here==0)
                                   {
                                   jast[tindex1][tindex2]=correl[d][tlab1*4+tlab2];
                                   tindex1=tindex1+1;
                                   tindex2=tindex2+1;
                                   jast[tindex1][tindex2]=correl[d][tlab1*4+tlab2];
                                   tindex1=tindex1-1;
                                   jast[tindex1][tindex2]=0;
                                   tindex1=tindex1+1;
                                   tindex2=tindex2-1;  
                                   jast[tindex1][tindex2]=0;
                                   }
                                    else if(flag_here==1)
                                    {
                                        jast[tindex1][tindex2]=correl[d][tlab2*4+tlab1];
                                        tindex1=tindex1+1;
                                        tindex2=tindex2+1;
                                        jast[tindex1][tindex2]=correl[d][tlab2*4+tlab1];
                                        tindex1=tindex1-1;
                                        jast[tindex1][tindex2]=0;
                                        tindex1=tindex1+1;
                                        tindex2=tindex2-1;
                                        jast[tindex1][tindex2]=0;
                                    }
                                 }
                             }   
                         }
                     }
                }
            } 
        }
    }   

  
    for(tlab1=0;tlab1<ntetra;tlab1++)
    {
        for(tlab2=tlab1;tlab2<ntetra;tlab2++)
        {
            if(abs(jast[tlab1][tlab2]-jast[tlab2][tlab1])>0.000000001){
            std::cout<<"tetra1 "<<tlab1<<" tetra2 "<<tlab2<<" jas "<<jast[tlab1][tlab2]<<" jas2 "<<jast[tlab2][tlab1]<<" diff "<<jast[tlab1][tlab2]-jast[tlab2][tlab1]<<"\n";
            }
        }
    } 
    return 0;
    pos=3;
    t=connect[pos][0];
    c1=qcharge(t,tetra,config);
    f1=c1-2*config[pos];
    t=connect[pos][1];
    c2=qcharge(t,tetra,config);
    f2=c2+2*config[pos];
    double prob,density,densitysquare,tempature,jp,estep;
    pos=3;
    prob=0.1;
    density=0.4;
    jp=0.08;
    densitysquare=pow(density,2.0);
    //
    energy_est(config,tetra,connect,L,density,jp,estep);
    flag=0;
    //test
    int pos2=myrand->randInt(5);
    pair_flip(config,ivic,tetra,connect,L,densitysquare,pos,pos2,prob,chargepairs,correl,myrand);
    /*Specific heat: single spin flip----------------------------------------------------*/
    //define a file stream.
    double esquare=0,eclassical=0;
    ofstream spec_heat;
    spec_heat.open("/Users/zhao/Documents/XCode/vmc_spin_ice/vmc_spin_ice/spec.txt");
    //thermalization
    flag=1;
    for(count=0;count<total*zpro;count++)
    {
        seedin=myrand->rand();
        //seedin=30.0;
        //singlespin_sweep(config,L,tempature,flag,seedin);
        singlespin_sweep_new(config,ivic,tetra,connect,L,tempature,flag,myrand);
    }
    for(tempature=0.1;tempature<3.0;tempature+=0.05)
    {
        //We are trying to measure the specific heat
        esquare=0;
        eclassical=0;
    for(count=0;count<stat;count++)
    {
        //update
        for(x0=0;x0<total;x0++)
        {
            seedin=myrand->randInt();
            singlespin_sweep_new(config,ivic,tetra,connect,L,tempature,flag,myrand);
        }
        //now we measure
        //e0total(config,tetra,L,estep);
        eclassical+=estep;
        estep=pow(estep,2.0);
        esquare+=estep;
    }
    //we now compute specific heat per cubic unit cell.
        eclassical=eclassical/((double)stat);
        esquare=esquare/((double)stat);
        cp=(esquare-pow(eclassical,2.0))/tempature/(double)(pow(L,3));
        eclassical=eclassical/(double)(pow(L,3));
        spec_heat<<tempature<<'\t'<<eclassical<<'\t'<<cp<<'\n';
    }
    spec_heat.close();
    return 0;
}

