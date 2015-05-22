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

//no drift
void energy_est_nodrift(int *neighbour,int neighbour_connect[][2],int first_neighbour[][12],int &ntetra,int *config,int tetra[][4],int connect[][2],int &L,double &jp,double &pair_density,double &estep){
  //simplest form of vmc
  //cout<<"in est\n";
  estep=0.0;
  int x,y,z,pos,pos2,total=ntetra/2;
  int ztetra=8*pow(L,2),ytetra=8*L;
  int xx,yy;
  int cter,cc,temp,c1,c2,t1,t2,f1,f2,spin1,spin2,charge,cc1,cc2,l1,l2,lt1,lt2;
  int flag1,flag2,flag,temp1,temp2;
  int range=1;
  for(z=0;z<L;z++)
  {
    for(y=0;y<L;y++)
    {
      for(x=0;x<L;x++)
      {
        //diagonal part;
        for(cc=0;cc<8;cc++)
        {
          temp=z*ztetra+y*ytetra+x*8+cc;
          charge=qcharge(temp,tetra,config);
          estep+=pow((double)charge,2.0)/8.0;
          //off-diagonal part;
          for(cc1=0;cc1<4;cc1++)
          {
            for(cc2=cc1+1;cc2<4;cc2++)
            {//c1 and c2 are the two spins on the particular tetrahedron.
              spin1=tetra[temp][cc1];
              spin2=tetra[temp][cc2];
              //if the bond is not flippable, next pair
              if(config[spin1]==config[spin2]){continue;}
              else{
                //we have a pair of flippalbe spins. What are locations of two tetrahedrons?
                if(connect[spin1][0]!=temp)
                {
                  t1=connect[spin1][0];
                }
                else
                {
                  t1=connect[spin1][1];
                }
                //similar for spin2
                if(connect[spin2][0]!=temp)
                {
                  t2=connect[spin2][0];
                }
                else
                {
                  t2=connect[spin2][1];
                }
                //now we only consider case
                c1=qcharge(t1,tetra,config);
                c2=qcharge(t2,tetra,config);
                config[spin1]=-config[spin1];
                config[spin2]=-config[spin2];
                f1=qcharge(t1,tetra,config);
                f2=qcharge(t2,tetra,config);
                config[spin1]=-config[spin1];
                config[spin2]=-config[spin2];
                flag1=0;
                flag2=0;
                flag=0;
                //only consider two cases
                if((c1==0)&&(c2==0))
                {//the final state have two charges.
                  estep+=-jp*0.5*pair_density;
                }
                else if((c1==-2)&&(c2==2)&&(f1==0)&&(f2==0))
                {
                  for(xx=0;xx<12;xx++)
                  {
                    if(first_neighbour[t1][xx]==t2){continue;}
                    l1=first_neighbour[t1][xx];
                    temp1=qcharge(l1,tetra,config);
                    if(temp1==2){flag1=1;}
                    for(yy=0;yy<12;yy++){
                      if(first_neighbour[t2][yy]==t1){continue;}
                      l2=first_neighbour[t2][yy];
                      temp2=qcharge(l2,tetra,config);
                      if(temp2==-2){flag2=1;}
                      if((flag1==1)&&(flag2==1))
                      {
                        lt1=l1/2;
                        lt2=l2/2;
                        if(neighbour[lt1*total+lt2]>range){flag=1;}
                      }
else if(neighbour[lt1*total+lt2]==1){
              //in this case, we need to see the two sites. 
             for(cter=0;cter<12;cter++){
                if(l2==first_neighbour[l1][cter]){break;}
             }
             if(config[neighbour_connect[l1*12+cter][0]]==config[neighbour_connect[l1*12+cter][1]]){flag=1;}
            }
                      flag2=0;
                    }
                    flag1=0;
                  }
                  if(flag==0)
                  {
                  estep+=-jp*0.5/pair_density;
                  }
                }
                else if((c1==2)&&(c2==-2)&&(f1==0)&&(f2==0))
                { 
                  for(xx=0;xx<12;xx++)
                  {
                    if(first_neighbour[t1][xx]==t2){continue;}
                    l1=first_neighbour[t1][xx];
                    temp1=qcharge(l1,tetra,config);
                    if(temp1==-2){flag1=1;}
                    for(yy=0;yy<12;yy++){
                      if(first_neighbour[t2][yy]==t1){continue;}
                      l2=first_neighbour[t2][yy];
                      temp2=qcharge(l2,tetra,config);
                      if(temp2==2){flag2=1;}
                      if((flag1==1)&&(flag2==1))
                      {
                        lt1=l1/2;
                        lt2=l2/2;
                        if(neighbour[lt1*total+lt2]>range){flag=1;}
                      }
else if(neighbour[lt1*total+lt2]==1){
              //in this case, we need to see the two sites. 
             for(cter=0;cter<12;cter++){
                if(l2==first_neighbour[l1][cter]){break;}
             }
             if(config[neighbour_connect[l1*12+cter][0]]==config[neighbour_connect[l1*12+cter][1]]){flag=1;}

            }
                      flag2=0;
                    }
                    flag1=0;
                  }
                  if(flag==0){
                  estep+=-jp*0.5/pair_density;
                  }
                }
                else{
                  continue;
                }
              }
            }
          }
        }
      }
    }
  }
  //cout<<"out est\n";
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

//A new function to compute the # of nth neighbour pairs times the product of the charges on this two pairs. 

void pn_count(int &ntetra,int *spinonc,int *neighbour,int *pn)
{//pn is a vector of size dist.size() and each element is the number of + - charge pairs at nth neighbour. 
  int x,y,temp=ntetra/2,xtemp,ytemp;
  //even sublattice
  for(x=0;x<ntetra;x+=2){
     for(y=x;y<ntetra;y+=2){
       	if(x==y){continue;}
       	if((spinonc[x]==1)&&(spinonc[y]==-1)){
	  pn[neighbour[(x*temp+y)/2]]+=-1; 
       	}
	else if((spinonc[x]==-1)&&(spinonc[y]==1)){
          pn[neighbour[(x*temp+y)/2]]+=-1;
	}
     }
  }
  //odd sublattice
  for(x=1;x<ntetra;x+=2){
     for(y=x;y<ntetra;y+=2){
	if(x==y){continue;}
	xtemp=(x-1)/2;
	ytemp=(y-1)/2;
       	if((spinonc[x]==1)&&(spinonc[y]==-1)){
	  pn[neighbour[xtemp*temp+ytemp]]+=-1; 
       	}
	else if((spinonc[x]==-1)&&(spinonc[y]==1)){
          pn[neighbour[xtemp*temp+ytemp]]+=-1;
	}	
     }
  }
  return;
}

/*We are trying to write the center piece of the code. From a set of mu_n, we are trying to determine fn
void fn_evaluate(int *config,int *mu_n,int *neighbour,int *spinonc,int *tetra,int *connect,int &L,int &ntetra,MTRand,*myrand,int *fn)
{
  //variable definition. 
  int x,y,z;
  //use mu_n,neighbour,we construct jastrow
  return;
}*/
