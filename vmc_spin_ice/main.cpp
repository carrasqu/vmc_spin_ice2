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
#include "MersenneTwister.h"
#include "routines.h"
#include "estimator.h"


int main()
{

    // insert code here...
    // system size L*L*L*4*4
    int L=1;
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
    //construction of tables.
    latt(ivic,tetra,connect,L,nh,ntetra);
    int pos;
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
    //Now we are trying to test the energy estimator.
    //some simple test of pair update.
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
    pair_flip(config,ivic,tetra,connect,L,densitysquare,pos,pos2,prob);
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
        e0total(config,tetra,L,estep);
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

