// pointers as arguments:
#include <iostream>
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "MersenneTwister.h"
using namespace std;

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

    cout<<ntetra<<"\n"<<tetracount<<"\n";
    cout<<"sites belongin to each tetrahedra z"<<"\n";
    for (z=0;z<ntetra;z++)
    {
      cout<<z<<"    "<<tetra[z][0]<<" "<<tetra[z][1]<<" "<<tetra[z][2]<<" "<<tetra[z][3]<<"\n"; 
    }    
    cout<<"\n"; 


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
    
   cout<<" site x connects which tetrahedra"<<"\n"; 
   for (z=0;z<nh;z++)
   {
     cout<<"x="<<z<<"   tetrahedra "<<connect[z][0]<<" "<<connect[z][1]<<"\n";
   } 

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
               
void loopupdate(int *config, int ivic[][6],int tetra[][4],int connect[][2], int &L, int &nh, int &ntetra)
{
 int itetra,v,ctetra,ii,vstart,vfinal,tet,toflip;
 int go = 0,ar,outs,tcharge;
 int spin[2],coun,countervisits = 0;
 int visited[ntetra][3]; // 0: visited or not. 1: order in the line of visits. 2: which spin was visited.
 int ordered[ntetra];   //  which tetrahedra are visited in order

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

 cout<<"total charge in configuration= "<<tcharge<<"\n";

 itetra=rand() % ntetra;
 ctetra=itetra;
 cout<<"random tetrahedron"<<itetra<<"\n";


do {

    cout<<"visited ctetra= "<<ctetra<<"visited="<< visited[ctetra][0]<<"  \n";
    if (visited[ctetra][0]==-1)
    {
       
       visited[ctetra][0]=1;
       ordered[countervisits]=ctetra; 
       visited[ctetra][1]=countervisits;
       countervisits=countervisits+1;

       cout<<"defected?=    "<<defected(ctetra,config,tetra)<<"     "; 
       if (defected(ctetra,config,tetra)==0)
       {
        
         // if not defected then randomly choose outward spin and locate tetrahedron located through that spin   
         if (ctetra%2==0)
         {
         // even tetrahedron: up (+1) means in; down (-1) means out.
          inoutspin(ctetra,config,tetra,spin);
          outs=rand()%2;
          cout<<"out spins"<<spin[0]<<spin[1]<<"chosen"<<spin[outs]<<"\n";
          visited[ctetra][2]=spin[outs];
          if (connect[spin[outs]][0]==ctetra)
          { 
             ctetra=connect[spin[outs]][1];
          }
          else if (connect[spin[outs]][1]==ctetra)
          { 
             ctetra=connect[spin[outs]][0];      
          }
          cout << "next chosen tetrahedron "  <<ctetra<<"\n";

         }
         else if (ctetra%2==1)
         {
         // odd tetrahedron: up (+1) means out; down (-1) means in.
          inoutspin(ctetra,config,tetra,spin);
          outs=rand()%2; 
          cout<<"out spins"<<spin[0]<<spin[1]<<"chosen"<<spin[outs]<<"\n";   
          visited[ctetra][2]=spin[outs];
          if (connect[spin[outs]][0]==ctetra)
          { 
           ctetra=connect[spin[outs]][1];
          } 
          else if (connect[spin[outs]][1]==ctetra)
          {
            ctetra=connect[spin[outs]][0];   
          }
           cout << "next chosen tetrahedron "  <<ctetra<<"\n"; 
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
 
   cout<<"  counting visits  "<<countervisits <<"  \n";

   } while (go==0);
    



   // flipping the spins along the loop
   if (go==1)
   {  
     cout<<"THERE WAS A LOOP vstart"<<vstart<<"vfinal="<<vfinal<<"\n";
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

 cout<<"total charge in final configuration= "<<tcharge<<"\n";


    return;
}


int main ()
{
  int L;
  int nh; 
  int ntetra;
  int nloops;
  int i;  
 
  cin>>L; 
  cout<<L<< "\n";
  nh=pow(L,3)*16;
  ntetra=pow(L,3)*2*4;
 
  int *config;
  config=new int[nh];  

  config[0]=-1;
  config[1]=-1;
  config[2]=-1;
  config[3]=1;
  config[4]=1;
  config[5]=-1;
  config[6]=-1;
  config[7]=1;
  config[8]=1;
  config[9]=-1;
  config[10]=-1;
  config[11]=1;
  config[12]=1;
  config[13]=-1;
  config[14]=-1;
  config[15]=1;

 
  int ivic[nh][6];  // each site its 6 neighbours 
  int tetra[ntetra][4]; // each tetrahedron and their 4 sites
  int connect[nh][2]; // each site connects two tetrahedra 

  /* initialize random seed: */
  srand (time(NULL));
  
  //construct the lattice tables 
  latt(ivic,tetra,connect,L,nh,ntetra);

  cout<<"neigbors"<<"\n"; 
  for (int n=0; n<nh; n++) {
      //for (int vic=0;vic<6;vic++){

      //cout<<n<<ivic[n][vic]<<"\n";
      cout<<n+1<<"    "<<ivic[n][0]+1<<" "<<ivic[n][1]+1<<" "<<ivic[n][2]+1<<" "<<ivic[n][3]+1<<" "<<ivic[n][4]+1<<" "<<ivic[n][5]+1<<"\n";    
      //}   

  }

  cout <<"before loop";
  nloops=1000;
  for (i=0;i<nloops;i++)  
  { 
   printv(config,nh); 
   loopupdate(config,ivic,tetra,connect,L,nh,ntetra);
   printv(config,nh); 
   cout<<"\n"<<"\n";
    
    
  }
  return 0;
}
