// 170715
// Delaunay lattice creation
#include"Matrix.h"
#include"EX.h"
#include"Graph.h"

void Draw_Structure(const Matrix &D, const Table &Tri);

// returns true if the triangle is traversed clockwise
// triangle with indices T, coordinates taken from Matrix D
bool Is_Clockwise(const Matrix &D, const List &T);

// Find if point P is inside the CC of triangle Tri,
// taking the coordinates from rows of matrix D
bool Is_Inside_CC(const Matrix &D, const List &T, const Vector &P);

// insert point inew from P into the triangle structure, change table T
void Bowyer_Watson_Step(Table &Tri, const Matrix &P, long inew);

Table Delaunay_Triangulation(Matrix &P);

void Triangulation_2_Graph(Graph &G, const Table &Tri, long N);

Graph Delaunay_Graph(const Matrix &P);


