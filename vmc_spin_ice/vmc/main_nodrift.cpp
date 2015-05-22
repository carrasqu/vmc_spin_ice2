/*We are implementing the "statistical minimization" developed in Sorella's paper for our VMC*/
/*We are taken a very simple wave function now. */
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <Accelerate/Accelerate.h>
#include "MersenneTwister.h"
#include "routines.h"
#include "estimator.h"

int main()
{
    //define variables
   int L=2;
   double prob,pair_density,tempature,jp;
   double t_tilde=0.01,eta=1.0;
   double density_low,density_high,t_tilde_low,t_tilde_high,density_step,t_tilde_step;
   double seedin=20.0;
   double disvec[4][3]={{0.0,0.0,0.0},{0.0,0.5,0.5},{0.5,0.0,0.5},{0.5,0.5,0.0}};
   int thbins; // number of bins during thermalization
   int nbins; // number of bins production run
   int msteps; // bin  length
   int nloops,error; // number of loops per monte carlo step
   int k,i,lo,l1; // counters
   int cloop,went,avvisit,visits; //
   int pos,pos2,max_iteration;
   double ave,delta_t,scale;
   int ndat=2;
   double data[ndat],data2[ndat];
   double esquare=0,eclassical=0,estep,estep2;
   int x,y,z,temp,dtemp;
    //we are trying to use a grid to get all the data in one simulation
   std::ifstream para;
   para.open("parameters.txt");
   para >> L >> jp >> seedin;
   //para >> density_low >> density_high >> density_step >> t_tilde_low >> t_tilde_high >> t_tilde_step;
   //para >> eta;
   para >> thbins >> nbins >> msteps >>nloops;
   para >> pair_density>>max_iteration;
    //cout << "give me L"<<"\n";
   para.close();
   cout << "L is " <<L<<"\n";
   MTRand *myrand=new MTRand();
   myrand->seed(seedin);
   int ntetra=8*L*L*L,nh=2*ntetra;
    //Tables
   int *config,*neighbour,*spinonc,*p_n,*pn_temp;
   config=new int[nh];
   neighbour=new int[ntetra*ntetra/4];
   spinonc=new int[ntetra];
   int ivic[nh][6],tetra[ntetra][4],connect[nh][2],first_neighbour[ntetra][12],neighbour_connect[ntetra*12][2];
   double table[ntetra][2];
   vector<int> coord;
   vector<double> dist;
   initialize_spinice_q0(config,spinonc,L,ntetra);
   //fill all tables. 
   latt(ivic,tetra,connect,L,nh,ntetra);
   count_neighbour(L,disvec,neighbour,dist,coord,first_neighbour);
   first_neighbour_table(ntetra,first_neighbour,neighbour_connect,connect,tetra);
   cout<<"I am here!\n";
   /*for(x=0;x<ntetra;x++){
       cout<<"tetra "<<x<<" 's first neighbours are ";
     for(y=0;y<12;y++)
     {
       cout<<first_neighbour[x][y]<<" ";
     }
     cout<<"\n";
   }
   for(x=0;x<ntetra/2;x++){
     for(y=0;y<ntetra/2;y++){
       cout << "first tetra is "<<x <<" second tetra is "<<y <<" they are "<<neighbour[x*ntetra/2+y]<<"th neighbours\n";
     }
   }*/
   //return 0;
   /*-----------------------------------------
     We fix this to be first neighbour
     -----------------------------------------*/
   //temp=2;
   //also initialize the smatrix
   /*start of the core of the minimization program.*/
     //now we start our "minimization process"
     //1. bulid jastrow CHECKED
     //2. thermailization
     for(i=0;i<thbins;i++){
       //compute jastrow from scratch;
       for(k=0;k<msteps;k++){
         //pair spin flip
         for(lo=0;lo<nh;lo++){
           pos=myrand->randInt(nh-1);//random spin
           pos2=myrand->randInt(5);//random neighbour
           prob=myrand->rand();
           pair_flip_nodrift(neighbour,neighbour_connect,first_neighbour,ntetra,config,spinonc,ivic,tetra,connect,L,prob,pos,pos2,pair_density);
         }
     /*look in
     for(x=0;x<ntetra;x++){cout<<"tetra "<<x<<" has charge "<<spinonc[x]<<"\n";}
     //we finished thermailization
     cout<<dist.size()<<"\n";
     pn_count(ntetra,spinonc,neighbour,pn_temp);
     for(x=1;x<dist.size();x++){
       cout<<"# of "<<x<<"th neighbour pairs is "<<pn_temp[x]<<"\n";
       pn_temp[x]=0;
     }*/
     //return 0;
         //loop update
         avvisit=0;
         cloop=0;
         for(lo=0;lo<nloops;lo++){
           loopupdate(config,ivic,tetra,connect,L,nh,ntetra,visits,went,myrand);
           if(went==1)
           {
             avvisit=avvisit+visits;
             cloop=cloop+1;  
           }
         }
         //adjust the required number of loops. 
         if(cloop>0){
           ave=(double)avvisit/(double)cloop;
           if(nloops<nh/(2*(int)ave)){
             nloops++;
           }
           else{
             nloops--;
             if(nloops<1)nloops=1;
           }
         }
       }
     }
     //initialize measurements
     eclassical=0.0;
     estep=0.0;
     for(i=0;i<ndat;i++){
       //holder for energy and energy square and their uncertainties. 
       data[i]=0.0;
       data2[i]=0.0;
     }
     for(i=0;i<nbins;i++){
       //again compute jastrow
       //reset energy and esquare as well
       eclassical=0.0;
       esquare=0.0;
       //a "Monte Carlo Step".
       for(k=0;k<msteps;k++){
         //reset pn_temp;
         //pair flip
         for(lo=0;lo<nh;lo++){
           pos=myrand->randInt(nh-1);
           pos2=myrand->randInt(5);
           prob=myrand->rand();
           pair_flip_nodrift(neighbour,neighbour_connect,first_neighbour,ntetra,config,spinonc,ivic,tetra,connect,L,prob,pos,pos2,pair_density);
         }
         //loop flip
         for(lo=0;lo<nloops;lo++){
           loopupdate(config,ivic,tetra,connect,L,nh,ntetra,visits,went,myrand);
         }
         energy_est_nodrift(neighbour,neighbour_connect,first_neighbour,ntetra,config,tetra,connect,L,jp,pair_density,estep);
         eclassical+=estep;
         //cout<<estep<<"\t";
         estep2=pow(estep,2.0);
         esquare+=estep2;
         //look in
         //measure p_ntemp
       }
       //return 0;
       //Here we need to somehow collect data. How should we do it?
       //we need to compute average of E and p_n. Then, we want to calculate f_n, which is generaled force. 
       eclassical=eclassical/((double)msteps);
       cout<<"bin "<<i<<" "<<"average of e is"<<eclassical<<"\n";
       esquare/=(double)msteps;
       data[0]+=eclassical;
       data2[0]+=pow(eclassical,2.0);
       if(i==36)
       {
         //we get to know this configuratoin!
         for(x=0;x<ntetra;x++){
           cout<<"tetra hedra "<<x<<" the spins are "<<config[tetra[x][0]]<<" "<<config[tetra[x][1]]<<" "<<config[tetra[x][2]]<<" "<<config[tetra[x][3]]<<" charge is";
           int cccc=config[tetra[x][0]]+config[tetra[x][1]]+config[tetra[x][2]]+config[tetra[x][3]];
           cout<< cccc<<"\n";
         }
         cout<<neighbour[22*(ntetra/2)+25]<<"\n";
         return 0;
       }
     }
     data[0]/=(double)nbins;
     data2[0]/=(double)nbins;
     data2[0]-=pow(data[0],2.0);
     data2[0]=pow(data2[0],0.5);
     cout<<"energy is: "<<data[0]<<" pm "<<data2[0]<<"\n";
     //compute the generalized force and change to mu
     //smatrix
     //cout<<"smatrix:"<<"\n";
     //cout<<"\n";
     //inverse the smatrix.
   return 0;
}
