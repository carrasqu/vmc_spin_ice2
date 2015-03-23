//
//  estimator.h
//  vmc_spin_ice
//
//  Created by Zhihao Hao on 2014-12-12.
//  Copyright (c) 2014 Zhihao Hao. All rights reserved.
//
//-----------------------------
//We are creating a separate header file for estimators.
//-----------------------------

//energy estimator: give a spin configuration, we produce a energy
//density is the holder for the density of monopoles.
void energy_est(int *config,int tetra[][4],int connect[][2],int &L,double &density,double &jp,double &estep)
{
    //We go through each tetrahedron and compute both the diagonal and off-diagonal energies
    int x,y,z,ztetra,ytetra;
    ztetra=8*pow(L,2);
    ytetra=8*L;
    //0<=cc<8 count tetrahedra number
    int cc,temp,charge;
    //c1,c2: two adjacent spins in a tetrahedron.
    int c1,c2,t1,t2,f1,f2,spin1,spin2;
    //start accumulating
    estep=0;
    for(z=0;z<L;z++)
    {
        for(y=0;y<L;y++)
        {
            for(x=0;x<L;x++)
            {
                for(cc=0;cc<8;cc++)
                {
                    //the diagonal part of energy
                    //temp: tetrahedron number
                    temp=z*ztetra+y*ytetra+8*x+cc;
                    //diagonal energy
                    charge=qcharge(temp,tetra,config);
                    estep+=pow((double)charge,2.0)/8.0;
                    
                    //of diagonal ones. we go through every pair of nearest neighbor spins. 
                    for(c1=0;c1<4;c1++)
                    {
                        //c2>c1.
                        for(c2=c1+1;c2<4;c2++)
                        {
                            //the location of the two spins.
                            spin1=tetra[temp][c1];
                            spin2=tetra[temp][c2];
                            //We eliminate the case where the two spins are of the same state:
                            if(config[spin1]==config[spin2])
                            {
                                continue;
                            }
                            else
                            {
                                //We figure out all the charge states.
                                if(connect[spin1][0]!=temp)
                                {
                                    t1=connect[spin1][0];
                                }
                                else
                                {
                                    t1=connect[spin1][1];
                                }
                                if(connect[spin2][0]!=temp)
                                {
                                    t2=connect[spin2][0];
                                }
                                else
                                {
                                    t2=connect[spin2][1];
                                }
                                //Now we compute the charges.
                                f1=qcharge(t1,tetra,config);
                                f2=qcharge(t2,tetra,config);
                                //we don't involve charge 2 (4) states.
                                if((abs(f1)==4)&&(abs(f2)==4))
                                {
                                    continue;
                                }
                                //now t1 and t2 store charges.
                                t1=f1;
                                t2=f2;
                                if(temp%2==1)//In this case, temp is a "down tetrahedra",which means t1 and t2 are both up tetrahedron
                                {
                                    int test=config[spin1];
                                    f1+=-2*config[spin1];
                                    test=config[spin2];
                                    f2+=(-2)*test;
                                }
                                else
                                {
                            
                                    //In this case, temp is "up" tetrahedron and both t1 and t2 are down tetrahedron
                                    f1+=2*config[spin1];
                                    f2+=2*config[spin2];
                                }
                                //eliminate the state with charge 2 (or 4) monopoles.
                                if((abs(f1)==4)||(abs(f2)==4))
                                {
                                    continue;
                                }
                                //We should stop caring about initial state with charge 4 also.
                                //Now we consider all different allowed cases
                                if((t1==0)&&(t2==0))
                                {
                                    //this process increase the monopole number by 2.
                                    estep+=-0.5*jp*density;
                                }
                                else if((abs(t1)==2)&&(t2==0))
                                {
                                    //move a monopole.
                                    estep+=-0.5*jp;
                                }
                                else if((t1==0)&&(abs(t2)==2))
                                {
                                    //move a monopole
                                    estep+=-0.5*jp;
                                }
                                else if((abs(t1)==2)&&(abs(t2)==2))
                                {
                                    //decrease the monopole number by 2.
                                    estep+=-0.5*jp/density;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //we don't do any normalization here.
    return;
}

void energy_est3(int *config,int tetra[][4],int &ntetra,int ivic[][6],int connect[][2],int &L,int &nh,double &jp,double &estep,double table[][2],double jast[],int spinonc[])
{
    //We go through each tetrahedron and compute both the diagonal and off-diagonal energies
    int x,y,z,ztetra,ytetra,i,j,pos,pos2;
    ztetra=8*pow(L,2);
    ytetra=8*L;
    //0<=cc<8 count tetrahedra number
    int cc,temp,charge;
    //c1,c2: two adjacent spins in a tetrahedron.
    int c1,c2,t1,t2,f1,f2,spin1,spin2;
    double jastr,jasi,jasf;
    
    
    //potential energy
    estep=0;
    for(z=0;z<ntetra;z++)
    {
     estep+=pow((double)qcharge(z,tetra,config),2.0)/8.0;
     //cout<<"spinonc "<<z<<" "<<qcharge(z,tetra,config)/2<<" " <<spinonc[z]<<"\n";
    }
    //cout<<"table "<<"\n";
    //for(z=0;z<ntetra;z++)
    //{
    
     //cout<<"table 0 - 1  "<<z<<" "<<table[z][0]<<" " <<table[z][1]<<"\n";
    //} 
     
    
    //cout<<"diagonal energy "<< estep<< "\n";
    //kinetic energy
    
    for(i=0;i<nh;i++)
    {
      pos=i;	
      for(j=0;j<6;j++)
      {
        
        pos2=ivic[i][j];
        if(config[pos]*config[pos2]==1)
        {
         // cout<<"nospin flip"<<"\n";
          continue;
        }
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

        //cout<<"spins"<<"\n";
        //cout<<"positions "<<pos<<" "<<pos2<<"\n";
        //cout<<"conf "<<config[pos]<<" "<<config[pos2]<<"\n";

        c1=qcharge(t1,tetra,config);
        c2=qcharge(t2,tetra,config);

        config[pos]=-config[pos];
        config[pos2]=-config[pos2];
        f1=qcharge(t1,tetra,config);
        f2=qcharge(t2,tetra,config);
        config[pos]=-config[pos];
        config[pos2]=-config[pos2];

         //cout<<"position changes t1="<<t1<<" t2="<<t2<<"\n";

         //cout<<"charges c1 c2 f1 f2 "<<c1/2<<" "<<c2/2<<" "<<f1/2<<" "<<f2/2<<" \n";

        if(abs(f1)==4||abs(f2)==4)
        {
         // wave function for Q =\pm 2  charges is zero
          continue;
        }

       jastr=1.0;
       if(-(f1/2+c1/2)==-1)
       {
         jastr=jastr*exp(table[t1][1]*double(f1/2-c1/2));
       }
       else if(-(f1/2+c1/2)==1)
       {
         jastr=jastr*exp(table[t1][0]*double(f1/2-c1/2));
       }

       if(-(f2/2+c2/2)==-1)
       {
         jastr=jastr*exp(table[t2][1]*double(f2/2-c2/2));
       }
       else if(-(f2/2+c2/2)==1)
       {
         jastr=jastr*exp(table[t2][0]*double(f2/2-c2/2));
       }
  
       //cout<<"f1 f2 "<< f1 << " "<<f2 <<"\n";
       //cout<< "f1*f2 " << double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))<<"\n";
       //cout<<"c1 c2 "<< c1 << " "<<c2 <<"\n";
       //cout<< "c1*c2 " << double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2))<<"\n"; 

       jastr=jastr*exp( 
                        +jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))
                       
                        +jast[t1+ntetra*t2]*double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2)));
         
       //cout <<double(c1/2)*double(f2/2)*0.5*(1.0-double(c1/2)*double(f2/2))<<"\n";
       //cout <<double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))<<"\n";    
       //cout <<double(c2/2)*double(f1/2)*0.5*(1.0-double(c2/2)*double(f1/2))<<"\n";
       //cout <<double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2))<<"\n";  
       // cout<<exp(-jast[t1+ntetra*t2]*double(c1/2)*double(f2/2)*0.5*(1.0-double(c1/2)*double(f2/2))
       //                 +jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))
       //                 -jast[t1+ntetra*t2]*double(c2/2)*double(f1/2)*0.5*(1.0-double(c2/2)*double(f1/2))
       //                 +jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2)))<<"\n";   
        //cout<<"charges energy c1 c2 f1 f2 "<<c1<<" "<<c2<<" "<<f1<<" "<<f2<<" \n";

        //cout<<"-table[5,0] , table[7,0] "<<-table[5][0]<<" "<<table[7][0]<<"\n";

        //cout<<"jastrow rat "<<  jastr<< " manual "<< exp( -jast[3+ntetra*5]+ jast[3+ntetra*7]   )<<"\n" ;
       // jasi=jastrowfactor( jast,spinonc,ntetra);
       // spinonc[t1]=f1/2;
       // spinonc[t2]=f2/2;
       // jasf=jastrowfactor( jast,spinonc,ntetra);
       // spinonc[t1]=c1/2;
       // spinonc[t2]=c2/2;
        estep+=-0.25*jp*jastr;
       //  cout<<"KKKKKKKKKKKKKKcomparison jasi/jasf vs ratio"<< jasi/jasf<<" "<<jastr <<" \n";
       //estep+=-0.25*jp*pow(sqrt(density),abs(f1/2)+abs(f2/2))/(pow(sqrt(density),abs(c1/2)+abs(c2/2)))*jasf/jasi;
        

           
        //exit(0);

	
      }       
      
    } 
    
    //cout<<"total energy "<< estep<< "\n";
   return;
}
void energy_est2(int *config,int tetra[][4],int &ntetra,int ivic[][6],int connect[][2],int &L,int &nh,double &density,double &jp,double &estep,double table[][2],double jast[],int spinonc[])
{
    //We go through each tetrahedron and compute both the diagonal and off-diagonal energies
    int x,y,z,ztetra,ytetra,i,j,pos,pos2;
    ztetra=8*pow(L,2);
    ytetra=8*L;
    //0<=cc<8 count tetrahedra number
    int cc,temp,charge;
    //c1,c2: two adjacent spins in a tetrahedron.
    int c1,c2,t1,t2,f1,f2,spin1,spin2;
    double jastr,jasi,jasf;
    
    
    //potential energy
    estep=0;
    for(z=0;z<ntetra;z++)
    {
     estep+=pow((double)qcharge(z,tetra,config),2.0)/8.0;
     //cout<<"spinonc "<<z<<" "<<qcharge(z,tetra,config)/2<<" " <<spinonc[z]<<"\n";
    }
    //cout<<"table "<<"\n";
    //for(z=0;z<ntetra;z++)
    //{
    
     //cout<<"table 0 - 1  "<<z<<" "<<table[z][0]<<" " <<table[z][1]<<"\n";
    //} 
     
    
    //cout<<"diagonal energy "<< estep<< "\n";
    //kinetic energy
    
    for(i=0;i<nh;i++)
    {
      pos=i;	
      for(j=0;j<6;j++)
      {
        
        pos2=ivic[i][j];
        if(config[pos]*config[pos2]==1)
        {
         // cout<<"nospin flip"<<"\n";
          continue;
        }
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

        //cout<<"spins"<<"\n";
        //cout<<"positions "<<pos<<" "<<pos2<<"\n";
        //cout<<"conf "<<config[pos]<<" "<<config[pos2]<<"\n";

        c1=qcharge(t1,tetra,config);
        c2=qcharge(t2,tetra,config);

        config[pos]=-config[pos];
        config[pos2]=-config[pos2];
        f1=qcharge(t1,tetra,config);
        f2=qcharge(t2,tetra,config);
        config[pos]=-config[pos];
        config[pos2]=-config[pos2];

         //cout<<"position changes t1="<<t1<<" t2="<<t2<<"\n";

         //cout<<"charges c1 c2 f1 f2 "<<c1/2<<" "<<c2/2<<" "<<f1/2<<" "<<f2/2<<" \n";

        if(abs(f1)==4||abs(f2)==4)
        {
         // wave function for Q =\pm 2  charges is zero
          continue;
        }

       jastr=1.0;
       if(-(f1/2+c1/2)==-1)
       {
         jastr=jastr*exp(table[t1][1]*double(f1/2-c1/2));
       }
       else if(-(f1/2+c1/2)==1)
       {
         jastr=jastr*exp(table[t1][0]*double(f1/2-c1/2));
       }

       if(-(f2/2+c2/2)==-1)
       {
         jastr=jastr*exp(table[t2][1]*double(f2/2-c2/2));
       }
       else if(-(f2/2+c2/2)==1)
       {
         jastr=jastr*exp(table[t2][0]*double(f2/2-c2/2));
       }
  
       //cout<<"f1 f2 "<< f1 << " "<<f2 <<"\n";
       //cout<< "f1*f2 " << double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))<<"\n";
       //cout<<"c1 c2 "<< c1 << " "<<c2 <<"\n";
       //cout<< "c1*c2 " << double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2))<<"\n"; 

       jastr=jastr*exp( 
                        +jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))
                       
                        +jast[t1+ntetra*t2]*double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2)));
         
       //cout <<double(c1/2)*double(f2/2)*0.5*(1.0-double(c1/2)*double(f2/2))<<"\n";
       //cout <<double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))<<"\n";    
       //cout <<double(c2/2)*double(f1/2)*0.5*(1.0-double(c2/2)*double(f1/2))<<"\n";
       //cout <<double(c1/2)*double(c2/2)*0.5*(1.0-double(c1/2)*double(c2/2))<<"\n";  
       // cout<<exp(-jast[t1+ntetra*t2]*double(c1/2)*double(f2/2)*0.5*(1.0-double(c1/2)*double(f2/2))
       //                 +jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2))
       //                 -jast[t1+ntetra*t2]*double(c2/2)*double(f1/2)*0.5*(1.0-double(c2/2)*double(f1/2))
       //                 +jast[t1+ntetra*t2]*double(f1/2)*double(f2/2)*0.5*(1.0-double(f1/2)*double(f2/2)))<<"\n";   
        //cout<<"charges energy c1 c2 f1 f2 "<<c1<<" "<<c2<<" "<<f1<<" "<<f2<<" \n";

        //cout<<"-table[5,0] , table[7,0] "<<-table[5][0]<<" "<<table[7][0]<<"\n";

        //cout<<"jastrow rat "<<  jastr<< " manual "<< exp( -jast[3+ntetra*5]+ jast[3+ntetra*7]   )<<"\n" ;
       // jasi=jastrowfactor( jast,spinonc,ntetra);
       // spinonc[t1]=f1/2;
       // spinonc[t2]=f2/2;
       // jasf=jastrowfactor( jast,spinonc,ntetra);
       // spinonc[t1]=c1/2;
       // spinonc[t2]=c2/2;
        estep+=-0.25*jp*pow(sqrt(density),abs(f1/2)+abs(f2/2))/(pow(sqrt(density),abs(c1/2)+abs(c2/2)))*jastr;
       //  cout<<"KKKKKKKKKKKKKKcomparison jasi/jasf vs ratio"<< jasi/jasf<<" "<<jastr <<" \n";
       //estep+=-0.25*jp*pow(sqrt(density),abs(f1/2)+abs(f2/2))/(pow(sqrt(density),abs(c1/2)+abs(c2/2)))*jasf/jasi;
        

           
        //exit(0);

	
      }       
      
    } 
    
    //cout<<"total energy "<< estep<< "\n";
   return;
}
