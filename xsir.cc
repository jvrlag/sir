// 180318-200318
// SIR model

#include"Matrix.h"
#include"Graph.h"
#include"EX.h"
#include"delaunay.h"

void Draw(const Matrix &X, const Graph &G, const List &V)
{
     long xsize=EX_CW->width;
     EX_Clear();
     EX_Color(0.5,0.5,0.5);
     for (long il=1;il<=G.Nl;il++)
     {
	  long i,j;
	  G.Link_Sites(i,j,il);
	  EX_Line((long)round(xsize*X(i,1)),
		  (long)round(xsize*X(i,2)),
		  (long)round(xsize*X(j,1)),
		  (long)round(xsize*X(j,2)));
     }

     for (long i=1;i<=X.N1;i++)
     {
	  if (V(i)==0) EX_Color(1,1,1);
	  if (V(i)==1) EX_Color(1,0,0);
	  if (V(i)==2) EX_Color(0,0.5,0.5);
	  EX_Fill_Circle( (long)round(xsize*X(i,1)),
			  (long)round(xsize*X(i,2)), 6);
     }
     
     EX_Flush();
}

// 0: susceptible
// 1: infected
// 2: recovered or vaccinated
void Evolve(const Graph &G, List &V, double inf, double rec)
{
     long N=V.N;
     long i=Rand_I(1,N);
     if (V(i)==2) return; // no reinfection
     if (V(i)==1)
     {
	  if (Rand()<rec) V(i)=2;
	  return; // already sick or vaccined
     }
     // Now, V(i)=0
     List Neigh=G.Neighbours(i);
     double sum=0;
     for (long k=1;k<=Neigh.N;k++)
	  sum+=(V(Neigh(k))==1);
     if (Rand()<sum*inf) V(i)=1;
     
}

int main(int argc, char *argv[])
{

     Rand_Open(time(0));
     long N=40;
     long Nv=30;
     double inf=0.1; // infection
     double rec=0.1; // recovery
     long graphics=1;
     long Nt=100;
     Input(N,"-N",argc,argv);
     Input(Nv,"-Nv",argc,argv);
     Input(inf,"-inf",argc,argv);
     Input(rec,"-rec",argc,argv);
     Input(graphics,"-graphics",argc,argv);
     Input(Nt,"-Nt",argc,argv);

     

     Matrix X(N,2);
     for (long i=1;i<=N;i++)
     {
	  X(i,1)=Rand(0.05,0.95);
	  X(i,2)=Rand(0.05,0.95);
     }
     Graph G=Delaunay_Graph(X);
     List V(N);
     V.Set(0);
     // launch a bunch of vaccinated people
     for (long i=1;i<=Nv;i++)
	  V(i+1)=2;
     V(1)=1;

     if (graphics)
     {
	  long xsize=700;
	  EX_Start(400,1,xsize,xsize);
	  EX_Enable_Buffer();
     }
     
     printf("# N: %ld, Nv: %ld, inf: %g, rec: %g\n",N,Nv,inf,rec);
     printf("T S I R\n");
     for (long in=1;in<=Nt*N;in++)
     {
	  if (graphics)
	  {
	       Draw(X,G,V);
	       if (EX_Key_Pressed()) do{}while(EX_Read_Key()!=' ');
	  }
	  Evolve(G,V,inf,rec);
	  if (in%N==0)
	       printf("%ld %ld %ld %ld\n",in/N,
		      V.Count(0),V.Count(1),V.Count(2));
     }

     if (graphics)
     {
	  EX_Read_Key();
	  EX_Close();
     }
}
