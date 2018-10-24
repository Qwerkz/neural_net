#include <math.h>
#include <stdio.h>
/*   This creates a big network of Ne, Ni cells
     with distance-dependent coupling that is sparse*/
     
double lambdaee=.1;
double lambdaie=.05;
double lambdaei=.05;
double lambdaii=.05;
double pee=.45,pei=5,pie=1,pii=5;
int c[1000][500];
double w[1000][500];
double drand48();
int nc[500];
FILE *fp;
FILE *wt;
main()
{

  fp=fopen("cxee3.tab","w");
  wt=fopen("wxee3.tab","w");
  makecon(400,400,lambdaee,pee,1);
  fp=fopen("cxei3.tab","w");
  wt=fopen("wxei3.tab","w");
  makecon(80,400,lambdaei,pei,0);
  fp=fopen("cxie3.tab","w");
  wt=fopen("wxie3.tab","w");
  makecon(400,80,lambdaie,pie,1);
  fp=fopen("cxii3.tab","w");
  wt=fopen("wxii3.tab","w");
  makecon(80,80,lambdaii,pii,0);



}

makecon(int npre,int npost,double lambda,double pmax,int type)
{
  int i,j,k,l,kmax=0;
  double x,y,p,z;
  for(i=0;i<1000;i++)
    nc[i]=0;

  for(i=0;i<npost;i++){
    k=0;
    x=(double) i/(double)npost;
    for(j=0;j<npre;j++){
      y=(double) j/(double) npre;
      z=fabs(x-y)/lambda;
      switch(type){
      case 1:
	if(drand48()<pmax*exp(-z)){
	  c[i][k]=j;
	  w[i][k]=1.0;
	  k++;
	}
	break;
      case 0:
	if(z<=1){
	  c[i][k]=j;
	  w[i][k]=1.0;
	  k++;
	  }
	break;
      }

    }
    nc[i]=k;
    if(k>kmax)kmax=k;
  }
  printf("kmax=%d \n",kmax);
  for(i=0;i<npost;i++){
    k=nc[i];
    for(l=k;l<kmax;l++){
      c[i][l]=i;
      w[i][l]=0;
    }
  }
  fprintf(fp,"%d\n0\n%d\n",npost*kmax,npost*kmax-1);
  fprintf(wt,"%d\n0\n%d\n",npost*kmax,npost*kmax-1);
  for(i=0;i<npost;i++)
    {
      for(k=0;k<kmax;k++)
	{
	  fprintf(fp,"%d\n",c[i][k]);
	  fprintf(wt,"%g\n",w[i][k]);
	}
    }

  fclose(wt);
  fclose(fp);
}
