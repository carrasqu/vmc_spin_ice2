/*We are implementing the "statistical minimization" developed in Sorella's paper for our VMC*/

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
   double prob,density,densitysquare,tempature,jp;
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
   para >> delta_t>>max_iteration;
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
   int ivic[nh][6],tetra[ntetra][4],connect[nh][2];
   double table[ntetra][2];
   double *jast,*mu_n,*f_n,*data_f,*smatrix,*smatrix_temp;//jastrow table,jastrow vector and generalized force holder and the "stochasitic" matrix.
   jast=new double[ntetra*ntetra];
   //we use a "pre-initialization" to determine the # of different neighbours in here. 
   vector<double> dist;
   vector<int> coord;
   initialize_spinice_q0(config,spinonc,L,ntetra);
   //fill all tables. 
   latt(ivic,tetra,connect,L,nh,ntetra);
   count_neighbour(L,disvec,neighbour,dist,coord);
   /*look in
   for(x=0;x<ntetra/2;x++){
     for(y=0;y<ntetra/2;y++){
       cout<<"first tetra: "<< x<<" second tetra: "<<y<<" .they are "<<neighbour[x*ntetra/2+y]<<" th neighbour.\n";
     }
   }*/
   //allocate memory for mu_n,f_n,p_n,smatrix. 
   temp=dist.size();
//   cout<<temp<<"\n";
   mu_n=new double[temp];
   f_n=new double[temp];
   p_n=new int[temp];
   pn_temp=new int[temp];
   data_f=new double[temp];
   dtemp=temp-1;
   smatrix=new double[(temp-1)*(temp-1)];
   smatrix_temp=new double[(temp-1)*(temp-1)];
   scale=0.2;
   //we need an initial guess for mu_n
   for(x=0;x<temp;x++){
	   //
	   p_n[x]=0;
	   //mu_n[x]=myrand->rand();
	   f_n[x]=1.0;
           data_f[x]=0.2;
           data_f[x]=1.0;
	   if((x==1)||(x==2)){
		   mu_n[x]=3.48;
           }
           else{
             //set large values of mu for long-range pairs to isolate the effect of the nn pairs. 
             mu_n[x]=20;
           }
   }
   /*some other mu_n initialization*/
   mu_n[1]=3.443;
   mu_n[2]=4.557;
   mu_n[3]=5.342;
   mu_n[4]=6.1;
   mu_n[5]=7.0;
   /*-----------------------------------------
     We fix this to be first neighbour
     -----------------------------------------*/
   //temp=2;
   //also initialize the smatrix
   for(x=1;x<temp;x++){
     for(y=1;y<temp;y++){
       //
       smatrix[(x-1)*(temp-1)+(y-1)]=0;
     }
   }
   int flag;
   double fnmag=fn_mag(data_f,temp),tol=5*pow(10.0,-5);


   /*start of the core of the minimization program.*/
   while(fnmag>tol){
     //now we start our "minimization process"
     //1. bulid jastrow CHECKED
     build_jast(L,ntetra,mu_n,neighbour,jast);
     //2. thermailization
     for(i=0;i<thbins;i++){
       //compute jastrow from scratch;
      jastrow(ntetra,spinonc,table,jast);
       for(k=0;k<msteps;k++){
         //pair spin flip
         for(lo=0;lo<nh;lo++){
           pos=myrand->randInt(nh-1);//random spin
           pos2=myrand->randInt(5);//random neighbour
           prob=myrand->rand();
           pair_flip3(config,spinonc,ivic,tetra,connect,L,pos,pos2,prob,ntetra,table,jast);
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
     for(i=0;i<temp;i++){
     //holder for generalized force
       data_f[i]=0.0;
     }
     for(x=1;x<temp;x++){
       for(y=1;y<temp;y++){
         smatrix[(x-1)*(temp-1)+(y-1)]=0.0;
       }
     }
     /*now we start to measure.*/
     for(i=0;i<nbins;i++){
       //again compute jastrow
       for(x=0;x<temp;x++){
         //holder for pairs and generalized force for each bin. 
         p_n[x]=0;
         f_n[x]=0.0;
       }
       for(x=1;x<temp;x++){
         for(y=1;y<temp;y++){
           //we set smatrix_temp to be zero
           smatrix_temp[(x-1)*(temp-1)+(y-1)]=0.0;
         }
       }
       //reset energy and esquare as well
       eclassical=0.0;
       esquare=0.0;
       jastrow(ntetra,spinonc,table,jast);
       //a "Monte Carlo Step".
       for(k=0;k<msteps;k++){
         //reset pn_temp;
         for(x=0;x<temp;x++){pn_temp[x]=0;}
         //pair flip
         for(lo=0;lo<nh;lo++){
           pos=myrand->randInt(nh-1);
           pos2=myrand->randInt(5);
           prob=myrand->rand();
           pair_flip3(config,spinonc,ivic,tetra,connect,L,pos,pos2,prob,ntetra,table,jast);
         }
         //loop flip
         for(lo=0;lo<nloops;lo++){
           loopupdate(config,ivic,tetra,connect,L,nh,ntetra,visits,went,myrand);
         }
         energy_est3(config,tetra,ntetra,ivic,connect,L,nh,jp,estep,table,jast,spinonc);
         eclassical+=estep;
         //cout<<estep<<"\t";
         estep2=pow(estep,2.0);
         esquare+=estep2;
         //look in
         //measure p_ntemp
         pn_count(ntetra,spinonc,neighbour,pn_temp);
         for(x=1;x<temp;x++){
           p_n[x]+=pn_temp[x];
           f_n[x]+=-2.0*estep*(double)pn_temp[x];
           //cout<<pn_temp[x]<<"\t"<<-2.0*estep*(double)pn_temp[x];
           for(y=1;y<temp;y++){
             //accumulate smatrix_temp
             smatrix_temp[(x-1)*(temp-1)+(y-1)]+=pn_temp[x]*pn_temp[y];
         //    cout<<smatrix_temp[(x-1)*(temp-1)+(y-1)]<<" ";
           }
       //    cout<<"\n";
         }
         for(x=1;x<temp;x++){
           pn_temp[x]=0;
         }
       }
       //return 0;
       //Here we need to somehow collect data. How should we do it?
       //we need to compute average of E and p_n. Then, we want to calculate f_n, which is generaled force. 
       eclassical=eclassical/((double)msteps);
      // cout<<"average of e is"<<eclassical<<"\n";
       esquare/=(double)msteps;
       data[0]+=eclassical;
       data2[0]+=pow(eclassical,2.0);
       for(x=1;x<temp;x++){
         double pntemp=(double)p_n[x]/((double)msteps);
       //  cout<<x<<"th neighbour pairs is expected to be "<<pntemp<<"\n";
         //collect general force for each bin. 
         data_f[x]+=f_n[x]/((double)msteps)+2*eclassical*pntemp;
         for(y=1;y<temp;y++){
           //accumulate smatrix
           smatrix[(x-1)*(temp-1)+(y-1)]+=smatrix_temp[(x-1)*(temp-1)+(y-1)]/((double)msteps)-pntemp*(double)p_n[y]/((double)msteps);
           //cout<<smatrix[(x-1)*(temp-1)+(y-1)]<<" ";
         }
         //cout<<"\n";
       }
     }
     data[0]/=(double)nbins;
     data2[0]/=(double)nbins;
     data2[0]-=pow(data[0],2.0);
     data2[0]=pow(data2[0],0.5);
     cout<<"iteration: "<<flag<<" .energy is: "<<data[0]<<" pm "<<data2[0]<<"\n";
     //compute the generalized force and change to mu
     //smatrix
     //cout<<"smatrix:"<<"\n";
     for(x=1;x<temp;x++){
       for(y=1;y<temp;y++){
         smatrix[(x-1)*(temp-1)+(y-1)]/=(double)nbins;
     //    cout<<smatrix[(x-1)*(temp-1)+(y-1)]<<" ";
       }
       //cout<<"\n";
     }
     //cout<<"\n";
     //inverse the smatrix.
     if(temp>2){
     get_inverse(dtemp,smatrix,error);
     }//both matrices are N by N. get_inverse is our portal to compute inverse
     /*cout<<"inverse of smatrix is:"<<"\n";
     for(x=1;x<temp;x++){
       for(y=1;y<temp;y++){
         cout<<smatrix[(x-1)*(temp-1)+(y-1)]<<" ";
       }
       cout<<"\n";
     }
     cout<<"\n";*/
     for(x=1;x<temp;x++){
       data_f[x]/=(double)nbins;
     }
     for(x=1;x<temp;x++){
       if(temp==2){
         mu_n[x]+=data_f[x]*delta_t/(smatrix[0]);
       }
       else{//now the matrix dimension is larger than 1. We need to do a matrix inversion. 
         /*Naive*/ 
         //mu_n[x]+=data_f[x]*delta_t;
         /*Sorella*/
         //scale=myrand->rand();
         for(y=1;y<temp;y++){
           mu_n[x]+=smatrix[(x-1)*(temp-1)+(y-1)]*data_f[y]*delta_t;
         }
       }
       cout<<"generalized force is: "<<data_f[x]<<" the new jastrow is "<<mu_n[x]<<"the `statistical factor' is "<<smatrix[0]<<"\n";
     }
     cout<<"\n";
     flag++;
     fnmag=fn_mag(data_f,temp);
     if(flag>max_iteration)
     {
       cout<<"the algorithm has not converged after "<<max_iteration<<" number of tries!\n";
       break;
     }
   }
   return 0;
}
