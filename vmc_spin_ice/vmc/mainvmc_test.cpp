
//  main.cpp
//  vmc_spin_ice
//
//  Created by Zhihao Hao on 2014-08-21.
//  Copyright (c) 2014 Zhihao Hao. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include "MersenneTwister.h"
#include "routines.h"
#include "estimator.h"


int main()
{

    // insert code here...
    // system size L*L*L*4*4
    int L=2;
    double prob,density,densitysquare,tempature,jp;
    double t_tilde=0.01,eta=1.0;
    double density_low,density_high,t_tilde_low,t_tilde_high,density_step,t_tilde_step;
    double seedin=20.0;
    double disvec[4][3]={{0.0,0.0,0.0},{0.0,0.5,0.5},{0.5,0.0,0.5},{0.5,0.5,0.0}};
    int thbins; // number of bins during thermalization
    int nbins; // number of bins production run
    int msteps; // bin  length
    int nloops; // number of loops per monte carlo step
    int k,i,lo,l1; // counters
    int cloop,went,avvisit,visits; //
    double ave;
    int ndat=2;
    double data[ndat],data2[ndat];
    double esquare=0,eclassical=0,estep;

    //we are trying to use a grid to get all the data in one simulation
    std::ifstream para;
    para.open("parameters.txt");
    para >> L >> jp >> seedin;
    para >> density_low >> density_high >> density_step >> t_tilde_low >> t_tilde_high >> t_tilde_step;
    para >> eta;
    para >> thbins >> nbins >> msteps >>nloops;
    //cout << "give me L"<<"\n";
    para.close();
    cout << "L is " <<L<<"\n";
    std::ofstream output_data;
    output_data.open("results.txt");
    output_data << "density\t t_tilde\t bins\t Energy\t EnVariance\t VarianceofTotE\t VarianceOfV \n";
    int nh;
    int ntetra,flag_here;
    //density
    double rho;
    //spin configuration and
    MTRand *myrand=new MTRand();
    myrand->seed(seedin);
    //seedin=myrand->rand();
    int *config;
    double *alpha;
    double cp;
    config=new int[L*L*L*4*4];
    int *neighbour;
    vector<double> dist;
    vector<int> coord;
    //initialization!!!!!!!!!!!!!!!!!!!! TEST
    //Notice X only works with even L. 
    //initialize_spinice_X(config,L);
    //
    nh=16*(int)pow(L,3);
    ntetra=8*(int)pow(L,3);
    neighbour=new int[ntetra*ntetra/4];
    int ivic[nh][6];  // each site its 6 neighbours
    int tetra[ntetra][4]; // each tetrahedron and their 4 sites
    int connect[nh][2]; // each site connects two tetrahedra
    int spinonc[ntetra];//spinon configuration 
    // initialization of config and spinonc
    int pos,pos2;
    int c1,c2,f1,f2,t,flag;
    initialize_spinice_q0(config,spinonc,L,ntetra);
    latt(ivic,tetra,connect,L,nh,ntetra);
    //test of our neighbour figuring out routine. 
    count_neighbour(L,disvec,neighbour,dist,coord); 
    //out put some stuff;
    //Let's try the single spin update.
    int x,y,z,pyro,site,temp,count,x0,total=10;
    double jast[ntetra*ntetra]; // jastrow factor
    double *mu_n;
    //mu_n stores the nth neighbour jastrow factor
    mu_n=new double[dist.size()];
    //this is just for testing
    for(x=0;x<dist.size();x++)
    {
	mu_n[x]=0.1*(double)(x);
    }
    build_jast(L,ntetra,mu_n,neighbour,jast); 
    //print stuff
    for(x=0;x<4*L*L*L;x++)
    {
      for(y=0;y<4*L*L*L;y++)
      {
        cout<<"tetra 1 is: "<<x<<" tetra2 is: "<<y<<" they are "<<neighbour[x*4*L*L*L+y]<<" th neighbour\n";
	if((abs(jast[2*x*ntetra+2*y]-0.1*neighbour[4*x*L*L*L+y])>0.000000001)||(abs(jast[(2*x+1)*ntetra+2*y+1]-0.1*neighbour[4*x*L*L*L+y])>0.00000001)){cout<<"something is wrong at "<<x<<" and "<<y<<" \n"<<"\n";}
	//cout<<"even sublattice: "<<jast[2*x*ntetra+2*y]<<" odd sublattice: "<<jast[(2*x+1)*ntetra+2*y+1]<<"\n";
      }
    }
    for(x=0;x<dist.size();x++)
    {
      cout<<x<<"th neighbour distance is "<<dist[x]<<"\n";
      cout<<"the coordination number is "<<coord[x]<<"\n";
    }
    return 0;
    int stat=1000000;
    int zpro=16*pow(L,2);
    int ypro=16*L;
    int xpro=16;
    double table[ntetra][2];
    int llgg=pow(L,3);
    double correl[llgg][16];
    int x1,y1,z1,x2,y2,z2,d,loc1,loc2,tlab1,tlab2,tindex1,tindex2;
    // VMC data
    
    

    //looping
    for(density=density_low;density<density_high;density+=density_step){
        for(t_tilde=t_tilde_low;t_tilde<t_tilde_high;t_tilde+=t_tilde_step){
    //Now we for each t_tilde and density, we calculate energy.
    
    initialize_spinice_q0(config,spinonc,L,ntetra);
    //construction of tables.
    latt(ivic,tetra,connect,L,nh,ntetra);
    densitysquare=pow(density,2.0);


    // CONSTRUCTING THE JASTROW PROPOSED BY ZHIHAO
    
    
    //now we call the function
    spinon_correlation(correl, t_tilde, eta, L,density);
    zpro=pow(L,2);
    ypro=L;
    /*for(z=0;z<L;z++)
    {
        for(y=0;y<L;y++)
        {
            for(x=0;x<L;x++)
            {

                
                for(pyro=0;pyro<4;pyro++)
                {
                    for(site=0;site<4;site++)
                    {
                        std::cout<<"(x,y,z) is:"<<(double)x+disvec[pyro][0]-disvec[site][0]<<'\t'<<(double)y+disvec[pyro][1]-disvec[site][1]<<'\t'<<(double)z+disvec[pyro][2]-disvec[site][2]<<'\n';
                        std::cout<<"mu_1: "<<pyro<<", "<<"mu_2: "<<site<<"\t"<<correl[z*zpro+y*ypro+x][pyro*4+site]<<'\n';
                        
                        
                    }
                }
            }
        }
    }*/
    
   	
    //return 0;
    
  
 /*
    for (z1=0;z1<L;z1++)
    {
        for (y1=0;y1<L;y1++)
        {
            for (x1=0;x1<L;x1++)
            {
              loc1=x1+y1*L+z1*zpro; 

                for (z2=0;z2<L;z2++)
                {
                     for (y2=0;y2<L;y2++)
                     {
                         for (x2=0;x2<L;x2++)
                         {
                            
                             loc2=x2+y2*L+z2*zpro;
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
                                   jast[tindex1*ntetra+tindex2]=correl[d][tlab1*4+tlab2];
                                   tindex1=tindex1+1;
                                   tindex2=tindex2+1;
                                   jast[tindex1*ntetra+tindex2]=correl[d][tlab1*4+tlab2];
                                   tindex1=tindex1-1;
                                   jast[tindex1*ntetra+tindex2]=0;
                                   tindex1=tindex1+1;
                                   tindex2=tindex2-1; 
                                   jast[tindex1*ntetra+tindex2]=0;
                                   }
                                    else if(flag_here==1)
                                    {
                                        jast[tindex1*ntetra+tindex2]=correl[d][tlab2*4+tlab1];
                                        tindex1=tindex1+1;
                                        tindex2=tindex2+1;
                                        jast[tindex1*ntetra+tindex2]=correl[d][tlab2*4+tlab1];
                                        tindex1=tindex1-1;
                                        jast[tindex1*ntetra+tindex2]=0;
                                        tindex1=tindex1+1;
                                        tindex2=tindex2-1;
                                        jast[tindex1*ntetra+tindex2]=0;
                                    }
                                 }
                             }  
                         }
                     }
                }
            } 
        }
    }   
*/
    for (z1=0;z1<L;z1++)
    {
        for (y1=0;y1<L;y1++)
        {
            for (x1=0;x1<L;x1++)
            {
              loc1=x1+y1*L+z1*zpro; 

                for (z2=0;z2<L;z2++)
                {
                     for (y2=0;y2<L;y2++)
                     {
                         for (x2=0;x2<L;x2++)
                         {
                            
                             loc2=x2+y2*L+z2*zpro;
			     d=correl_index(loc1,loc2,L,flag_here);
                             for(tlab1=0;tlab1<4;tlab1++)
                             {
                                 for(tlab2=0;tlab2<4;tlab2++)
                                 {
                                   tindex1=z1*L*L*8+y1*L*8+x1*8+2*tlab1;
                                   tindex2=z2*L*L*8+y2*L*8+x2*8+2*tlab2;
                                   // std::cout<<"tetras considered "<<tindex1<<" "<<tindex2<<"\n" ;
                                   jast[tindex1*ntetra+tindex2]=eta*log(density);
                                   tindex1=tindex1+1;
                                   tindex2=tindex2+1;
                                   jast[tindex1*ntetra+tindex2]=eta*log(density);
                                   tindex1=tindex1-1;
                                   jast[tindex1*ntetra+tindex2]=0;
                                   tindex1=tindex1+1;
                                   tindex2=tindex2-1; 
                                   jast[tindex1*ntetra+tindex2]=0;
				 }
                             }  
                         }
                     }
                }
            } 
        }
    }   
    /*for(tlab1=0;tlab1<ntetra;tlab1++)
    {
        for(tlab2=tlab1;tlab2<ntetra;tlab2++)
        {
	   if(abs(jast[tlab1+ntetra*tlab2]-jast[tlab2+ntetra*tlab1])>pow(10.0,-8.0)){
           std::cout<<"tetra1 "<<tlab1<<" tetra2 "<<tlab2<<" jas "<<jast[tlab1+ntetra*tlab2]<<" and "<<jast[tlab2+ntetra*tlab1]<<"\n";
           }  
        }
    } */
    
    //return 0;
    
    tempature=0.0;     
    flag=0;
    e0total(config,tetra,ntetra,L,estep); 
    //cout <<"charge^2"<<estep<<"\n";
      
    // Thermalization
    seedin=myrand->randInt();
    for(i=0;i<thbins;i++)
    {
      //compute jastrow from scratch
       jastrow(ntetra,spinonc, table, jast); 
      //cout<<i<<" "<<nloops <<"\n";
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
         pair_flip3(config,spinonc,ivic,tetra,connect,L,pos,pos2,prob,ntetra,table,jast);
         }
         
//    e0total(config,tetra,ntetra,L,estep); 
    //cout <<"charge^2: "<<estep<<"\n";
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
    //cout<<i<<" "<<"avvisit "<<ave <<" nloops "<<nloops << "crit "<<nh/(2*(int)ave) <<"\n";

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

      //compute jastrow from scratch
       jastrow(ntetra,spinonc, table, jast);

      for(k=0;k<msteps;k++)
      {
         
         // pair spin flip
         for(lo=0;lo<nh;lo++)
         {
         pos=myrand->randInt(nh-1); // random spin  to flip
         pos2=myrand->randInt(5);   // random neighbor  to flip
         prob=myrand->rand();
         pair_flip3(config,spinonc,ivic,tetra,connect,L,pos,pos2,prob,ntetra,table,jast);
         }


         //singlespin_sweep_new(config,ivic,tetra,connect,L,tempature,flag,myrand); 
         // loopupdate
         for(lo=0;lo<nloops;lo++)
         {
          loopupdate(config,ivic,tetra,connect,L,nh,ntetra,visits,went,myrand);
         }
         energy_est3(config,tetra,ntetra,ivic,connect,L,nh,jp,estep,table,jast,spinonc);
         eclassical+=estep;
         //std::cout<<"the energy is "<<estep<<"\n";
         estep=pow(estep,2.0);
         esquare+=estep;   
      }
      collect2(tempature,eclassical,esquare,data,data2,ndat,nh,msteps,i,output_data,density,t_tilde);
    }
    }
    }
    output_data.close();
    return 0;
}

