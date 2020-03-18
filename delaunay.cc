// 170715-180319
// Delaunay lattice creation
#include"delaunay.h"

void Draw_Structure(const Matrix &D, const Table &Tri)
{
     long xsize=EX_CW->width;
     for (long i=1;i<=Tri.N1;i++)
     {
	  long p1=Tri(i,1);
	  long p2=Tri(i,2);
	  long p3=Tri(i,3);
	  EX_Line( D(p1,1)*xsize, D(p1,2)*xsize, D(p2,1)*xsize, D(p2,2)*xsize);
	  EX_Line( D(p2,1)*xsize, D(p2,2)*xsize, D(p3,1)*xsize, D(p3,2)*xsize);
	  EX_Line( D(p3,1)*xsize, D(p3,2)*xsize, D(p1,1)*xsize, D(p1,2)*xsize);
     }
}


// returns true if the triangle is traversed clockwise
// triangle with indices T, coordinates taken from Matrix D
bool Is_Clockwise(const Matrix &D, const List &T)
{
     Vector V1=D.Row(T(2))-D.Row(T(1));
     Vector V2=D.Row(T(3))-D.Row(T(2));
     double cross=V1(1)*V2(2)-V1(2)*V2(1);
     return (cross>0.0);
}

// Find if point P is inside the CC of triangle Tri,
// taking the coordinates from rows of matrix D
bool Is_Inside_CC(const Matrix &D, const List &T, const Vector &P)
{
     bool clockwise=Is_Clockwise(D,T);
     Matrix W(3);
     W(1,1)=D(T(1),1)-P(1);
     W(2,1)=D(T(2),1)-P(1);
     W(3,1)=D(T(3),1)-P(1);
     W(1,2)=D(T(1),2)-P(2);
     W(2,2)=D(T(2),2)-P(2);
     W(3,2)=D(T(3),2)-P(2);
     W(1,3)=Sqr(D(T(1),1))+Sqr(D(T(1),2))-Sqr(P(1))-Sqr(P(2));
     W(2,3)=Sqr(D(T(2),1))+Sqr(D(T(2),2))-Sqr(P(1))-Sqr(P(2));
     W(3,3)=Sqr(D(T(3),1))+Sqr(D(T(3),2))-Sqr(P(1))-Sqr(P(2));
     bool pos=(Det(W)>0.0);
     return (clockwise ? pos : !pos );
}

// insert point inew from P into the triangle structure, change table T
void Bowyer_Watson_Step(Table &Tri, const Matrix &P, long inew)
{
//     printf("New point: "); P.Row(inew).Write();
     long Nt=Tri.N1;
     List Bad; // Find bad triangles
     for (long it=1;it<=Nt;it++)
     {
	  if (Is_Inside_CC(P,Tri.Row(it),P.Row(inew)))
	       Bad.Append(it);
     }
//     printf("Bad triangles: ");
//     Bad.Write();
     // Now, make new triangles
     // Check all the edges in the triangles, if they do not appear in other
     // triangles, we create a new triangle
     long Ntbad=Bad.N;
     for (long itb=1;itb<=Ntbad;itb++) // for all bad triangles
     {
	  for (long k=1;k<=3;k++) // for all edges of the triangle
	  {
	       long l1=Tri( Bad(itb), k );
	       long l2=Tri( Bad(itb), k<3 ? k+1 : 1 );
	       bool new_triangle=true;
	       // is it present in other bad triangle?
	       for (long itb2=1;itb2<=Ntbad;itb2++)
		    if (itb2!=itb) // not in our triangle!
		    {
			 List Itb2=Tri.Row( Bad(itb2) );
			 if (Itb2.Find(l1) && Itb2.Find(l2))
			 {
			      new_triangle=false;
			      break;
			 }
		    }
	       if (new_triangle) // add new triangles!
	       {
		    List New_T(3); New_T(1)=l1; New_T(2)=l2; New_T(3)=inew;
		    Tri.Append_Row(New_T);
	       }
	  }
     }
     // now, remove bad triangles
     for (long itb=Ntbad;itb>=1;itb--)
	  Tri.Remove_Row(Bad(itb));
}

Table Delaunay_Triangulation(Matrix &P)
{
     long N=P.N1;
     P.Resize(N+3,2);
     P(N+1,1)=0.0; P(N+1,2)=2.0;
     P(N+2,1)=0.0; P(N+2,2)=0.0;
     P(N+3,1)=2.0; P(N+3,2)=0.0; // supertriangle

     Table Tri(1,3);
     // the big external triangle
     Tri(1,1)=N+1; Tri(1,2)=N+2; Tri(1,3)=N+3;
//     long Nt=1; // number of triangles
     
     for (long inew=1;inew<=N;inew++)
	  Bowyer_Watson_Step(Tri,P,inew);

     // Now, final structure
     // Remove all triangles containing points 1, 2 and 3
     for (long i=Tri.N1;i>=1;i--)
     {
	  List T=Tri.Row(i);
	  if ( T.Find(N+1) || T.Find(N+2) || T.Find(N+3) ) Tri.Remove_Row(i);
     }
   
     return Tri;
}

void Triangulation_2_Graph(Graph &G, const Table &Tri, long N)
{
//     Graph G(N);
     G.Create(N);
     for (long k=1;k<=Tri.N1;k++)
     {
          long s1=Tri(k,1), s2=Tri(k,2), s3=Tri(k,3);
          G.Add_Link(s1,s2);
          G.Add_Link(s2,s3);
          G.Add_Link(s3,s1);
     }
     G.Update_Index();
}

Graph Delaunay_Graph(const Matrix &P)
{
     Matrix P2(P);
     Table Tri=Delaunay_Triangulation(P2);
     Graph G;
     Triangulation_2_Graph(G,Tri,P.N1);
     return G;
}
