#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
  int nbVert;
  int *ListNbAdj;
  int **ListAdjacencies;
  int *BipartitionVector;
  double *VectA;
  double *VectB;
  double *VectC;
  double *TheDefect;
  double *TheDefectBis;
  int *SetPlus;
  double *TheRadiusA_Step1;
  double *TheRadiusA_Step2;
  double *MyAttempt;
  double *TheResult;
  double *RadiusVector;
  double *u2;
  double *u3;
  double *prov1;
  double *prov2;
  double *DirectionVect;
  double *SingleCoord;
  double *TheDiff;
  double **ListCoord;
  /* provisionnary variables */
  int *ListDone;
  int *ListStatus;

} PlanGraph;

void ReadGraph(char* TheFile, PlanGraph *TheGraph)
{
  FILE *DATAGRAPH=NULL;
  int nbVert, iVert, i, nbAdj;
  int ThePart;
  int TheVal;
  int eRet;
  DATAGRAPH=fopen(TheFile, "r");
  if (DATAGRAPH == NULL)
    {
      fprintf(stderr,"The file %s does not exist\n",TheFile);
      exit(EXIT_FAILURE);
    }
  eRet=fscanf(DATAGRAPH, "%d", &nbVert);
  TheGraph->nbVert=nbVert;
  if ((TheGraph->BipartitionVector = (int*)malloc(nbVert*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->SetPlus = (int*)malloc(nbVert*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->ListNbAdj = (int*)malloc(nbVert*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->ListDone = (int*)malloc(nbVert*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->ListStatus = (int*)malloc(nbVert*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->VectA = (double*)malloc(nbVert*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->VectB = (double*)malloc(nbVert*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->VectC = (double*)malloc(nbVert*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->TheDefect = (double*)malloc(nbVert*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->TheDefectBis = (double*)malloc(nbVert*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->TheRadiusA_Step1 = (double*)malloc(nbVert*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->TheRadiusA_Step2 = (double*)malloc(nbVert*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->MyAttempt = (double*)malloc(nbVert*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->TheResult = (double*)malloc(nbVert*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->RadiusVector = (double*)malloc(nbVert*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->u2 = (double*)malloc(3*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->u3 = (double*)malloc(3*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->prov1 = (double*)malloc(3*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->prov2 = (double*)malloc(3*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->DirectionVect = (double*)malloc(3*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->SingleCoord = (double*)malloc(3*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->TheDiff = (double*)malloc(3*sizeof(double))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->ListAdjacencies = (int**)malloc(nbVert*sizeof(int*))) == 0)
    exit (EXIT_FAILURE);
  if ((TheGraph->ListCoord = (double**)malloc(nbVert*sizeof(double*))) == 0)
    exit (EXIT_FAILURE);
  for (iVert=1; iVert<=nbVert; iVert++)
    {
      TheGraph->ListDone[iVert-1]=0;
      TheGraph->ListStatus[iVert-1]=0;
      if ((TheGraph->ListCoord[iVert-1] = (double*)malloc(3*sizeof(double))) == 0)
	exit (EXIT_FAILURE);
      for (i=1; i<=3; i++)
	TheGraph->ListCoord[iVert-1][i-1]=0;
    }
  for (iVert=1; iVert<=nbVert; iVert++)
    {
      eRet=fscanf(DATAGRAPH, "%d", &ThePart);
      TheGraph->BipartitionVector[iVert-1]=ThePart;
      eRet=fscanf(DATAGRAPH, "%d", &nbAdj);
      TheGraph->ListNbAdj[iVert-1]=nbAdj;
      if ((TheGraph->ListAdjacencies[iVert-1] = (int*)malloc(nbAdj*sizeof(int))) == 0)
	exit (EXIT_FAILURE);
      for (i=1; i<=nbAdj; i++)
	{
	  eRet=fscanf(DATAGRAPH, "%d", &TheVal);
	  TheGraph->ListAdjacencies[iVert-1][i-1]=TheVal;
	}
    }
  fclose(DATAGRAPH);
}




double Alpha(double Ru, double Rv)
{
  //  fprintf(stderr, "Ru=%f  Rv=%f\n", Ru, Rv);
  return atan(tan(Ru)/sin(Rv));
}


void difAlpha(double Ru, double Rv, double difRu, double difRv)
{
  double W, s1, s2;
  W=tan(Ru)/sin(Rv);
  s1=cos(Ru);
  s2=sin(Rv);
  difRu=1/(s2*s1*s1*(1+W*W));
  difRv=-tan(Ru)*cos(Rv)/(s2*s2*(1+W*W));
}



void DefectEvaluation(PlanGraph *TheGraph, double **TheRadius, double **TheDefect)
{
  double defect, pi, angle;
  long uVert, vVert;
  int i;
  pi=3.141592653589792;
  for (vVert=1; vVert<=TheGraph->nbVert; vVert++)
    {
      defect=pi;
      for (i=1; i<=TheGraph->ListNbAdj[vVert-1]; i++)
	{
	  uVert=TheGraph->ListAdjacencies[vVert-1][i-1];
	  angle=Alpha((*TheRadius)[uVert-1], (*TheRadius)[vVert-1]);
	  //	  fprintf(stderr, "angle=%f\n", angle);
	  defect=defect-angle;
	}
      (*TheDefect)[vVert-1]=defect;
    }
}


double SquareFunctional(PlanGraph *TheGraph, double **TheRadius)
{
  double defect, pi, angle;
  long uVert, vVert;
  double sqr;
  int i;
  pi=3.141592653589792;
  sqr=0;
  for (vVert=1; vVert<=TheGraph->nbVert; vVert++)
    {
      defect=pi;
      for (i=1; i<=TheGraph->ListNbAdj[vVert-1]; i++)
	{
	  uVert=TheGraph->ListAdjacencies[vVert-1][i-1];
	  angle=Alpha((*TheRadius)[uVert-1], (*TheRadius)[vVert-1]);
	  defect=defect-angle;
	}
      sqr=sqr+defect*defect;
    }
  return sqr;
}




void ScalarMultVector(double **InputVect, double **w, double val)
{
  int i;
  for (i=1; i<=3; i++)
    (*w)[i-1]=val*(*InputVect)[i-1];
}

void AdditionVector(double **V1, double **V2, double **w)
{
  int i;
  for (i=1; i<=3; i++)
    (*w)[i-1]=(*V1)[i-1]+(*V2)[i-1];
}

void VectorProduct(double **a, double **b, double **c)
{
  double c1, c2, c3;
  c1=(*a)[1]*(*b)[2]-(*a)[2]*(*b)[1];
  c2=(*a)[2]*(*b)[0]-(*a)[0]*(*b)[2];
  c3=(*a)[0]*(*b)[1]-(*a)[1]*(*b)[0];
  (*c)[0]=c1;
  (*c)[1]=c2;
  (*c)[2]=c3;
}

void RenormalizeVector(double **TheVect, double **RetVect)
{
  double TheNorm;
  int iCol;
  TheNorm=0;
  for (iCol=1; iCol<=3; iCol++)
    TheNorm=TheNorm+(*TheVect)[iCol-1]*(*TheVect)[iCol-1];
  for (iCol=1; iCol<=3; iCol++)
    (*RetVect)[iCol-1]=(*TheVect)[iCol-1]/sqrt(TheNorm);
}




double VectorScalarProduct(double **V1, double **V2)
{
  double scal;
  int i;
  scal=0;
  for (i=1; i<=3; i++)
    scal=scal+(*V1)[i-1]*(*V2)[i-1];
  return scal;
}






double FuncDistance(double Ru, double Rv)
{
  return acos(cos(Ru)*cos(Rv));
}




int NextIdx(int TheDeg, int idx)
{
  if (idx < TheDeg)
    return idx+1;
  return 1;
}



void VectorsU2U3(double **u, double **v, double **u2, double **u3)
{
  double w[3];
  double scal, norm;
  int i;
  scal=VectorScalarProduct(u, v);
  norm=0;
  for (i=1; i<=3; i++)
    {
      w[i-1]=(*v)[i-1]-scal*(*u)[i-1];
      norm=norm+w[i-1]*w[i-1];
    }
  for (i=1; i<=3; i++)
    (*u2)[i-1]=w[i-1]/sqrt(norm);
  VectorProduct(u, u2, u3);
}



void ComputeCoordinates(PlanGraph *TheGraph, double *MaximalDiscrepancy)
{
  double angle, TheDist, MinusOne, val1, val2, dist;
  int jVert, iAdj, PrevAdj, TheIAdj, i;
  int PrevVert, TheVert;
  int IsAllDone;
  int iVert;
  iVert=1;
  jVert=TheGraph->ListAdjacencies[0][0];
  TheDist=FuncDistance(TheGraph->RadiusVector[iVert-1], TheGraph->RadiusVector[jVert-1]);
  TheGraph->ListCoord[iVert-1][0]=1;
  TheGraph->ListCoord[jVert-1][0]=cos(TheDist);
  TheGraph->ListCoord[jVert-1][1]=sin(TheDist);
  TheGraph->ListStatus[iVert-1]=1;
  TheGraph->ListStatus[jVert-1]=1;
  *MaximalDiscrepancy=0;
  PrevAdj=-445;
  while(1)
    {
      IsAllDone=1;
      for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
	{
	  if (TheGraph->ListStatus[iVert-1] == 1 && TheGraph->ListDone[iVert-1] == 0)
	    {
	      TheGraph->ListDone[iVert-1]=1;
	      for (iAdj=1; iAdj<=TheGraph->ListNbAdj[iVert-1]; iAdj++)
		{
		  jVert=TheGraph->ListAdjacencies[iVert-1][iAdj-1];
		  if (TheGraph->ListStatus[jVert-1] == 1)
		    PrevAdj=iAdj;
		}
	      VectorsU2U3(&(TheGraph->ListCoord[iVert-1]), &(TheGraph->ListCoord[TheGraph->ListAdjacencies[iVert-1][PrevAdj-1]-1]), &(TheGraph->u2), &(TheGraph->u3));
	      angle=0;
	      for (iAdj=1; iAdj<=TheGraph->ListNbAdj[iVert-1]; iAdj++)
		{
		  TheIAdj=NextIdx(TheGraph->ListNbAdj[iVert-1], PrevAdj);
		  PrevVert=TheGraph->ListAdjacencies[iVert-1][PrevAdj-1];
		  TheVert=TheGraph->ListAdjacencies[iVert-1][TheIAdj-1];
		  angle=angle+Alpha(TheGraph->RadiusVector[PrevVert-1], TheGraph->RadiusVector[iVert-1])+Alpha(TheGraph->RadiusVector[TheVert-1], TheGraph->RadiusVector[iVert-1]);
		  val1=cos(angle);
		  ScalarMultVector(&(TheGraph->u2), &(TheGraph->prov1), val1);
		  val2=sin(angle);
		  ScalarMultVector(&(TheGraph->u3), &(TheGraph->prov2), val2);
		  AdditionVector(&(TheGraph->prov1), &(TheGraph->prov2), &(TheGraph->DirectionVect));
		  TheDist=FuncDistance(TheGraph->RadiusVector[iVert-1], TheGraph->RadiusVector[TheVert-1]);
		  val1=cos(TheDist);
		  ScalarMultVector(&(TheGraph->ListCoord[iVert-1]), &(TheGraph->prov1), val1);
		  val2=sin(TheDist);
		  ScalarMultVector(&(TheGraph->DirectionVect), &(TheGraph->prov2), val2);
		  AdditionVector(&(TheGraph->prov1), &(TheGraph->prov2), &(TheGraph->SingleCoord));
		  if (TheGraph->ListStatus[TheVert-1] == 0)
		    {
		      for (i=1; i<=3; i++)
			TheGraph->ListCoord[TheVert-1][i-1]=TheGraph->SingleCoord[i-1];
		      TheGraph->ListStatus[TheVert-1]=1;
		      IsAllDone=0;
		    }
		  else
		    {
		      MinusOne=-1;
		      ScalarMultVector(&(TheGraph->ListCoord[TheVert-1]), &(TheGraph->prov1), MinusOne);
		      AdditionVector(&(TheGraph->SingleCoord), &(TheGraph->prov1), &(TheGraph->TheDiff));
		      dist=VectorScalarProduct(&(TheGraph->TheDiff), &(TheGraph->TheDiff));
		      if (dist > *MaximalDiscrepancy)
			*MaximalDiscrepancy=dist;
		    }
		  PrevAdj=TheIAdj;
		}
	    }
	}
      if (IsAllDone == 1)
	break;
    }
}





void PartialMultVector(double **TheRadiusOld, int **SetPlus, double **TheResult, double MultPlus, int TheSize)
{
  int iVert;
  for (iVert=1; iVert<=TheSize; iVert++)
    if ((*SetPlus)[iVert-1] == 1)
      (*TheResult)[iVert-1]=(*TheRadiusOld)[iVert-1]*MultPlus;
    else
      (*TheResult)[iVert-1]=(*TheRadiusOld)[iVert-1];
}








void ScalarMult(PlanGraph *TheGraph, double **InputVect, double **OutputVect, double TheMult)
{
  int iVert;
  for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
    {
      //      fprintf(stderr, "Before i=%d  val=%f\n", iVert, (*InputVect)[iVert-1]);
      (*OutputVect)[iVert-1]=(*InputVect)[iVert-1]*TheMult;
      //      fprintf(stderr, " After i=%d  val=%f\n", iVert, (*InputVect)[iVert-1]);
    }
  //  fprintf(stderr, "TheMult=%f\n", TheMult);
}



double MaxMultiplier(double **TheRadius, int TheSize)
{
  double pi, TheMin, fact;
  int iVert;
  pi=3.141592653589792;
  TheMin=(pi/2)/(*TheRadius)[0];
  for (iVert=2; iVert<=TheSize; iVert++)
    {
      fact=(pi/2)/(*TheRadius)[iVert-1];
      if (fact < TheMin)
	TheMin=fact;
    }
  return TheMin;
}


double NormalizationValue(PlanGraph *TheGraph, double **TheRadius)
{
  int iVert;
  double res;
  DefectEvaluation(TheGraph, TheRadius, &(TheGraph->TheDefectBis));
  res=0;
  for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
    {
      res=res+TheGraph->TheDefectBis[iVert-1];
      //      fprintf(stderr, "defect=%f\n", TheGraph->TheDefectBis[iVert-1]);
    }
  return res;
}



void SimpleRenormalization(PlanGraph *TheGraph, double **InputVect, double **OutputVect)
{
  double a, b, c, ValA, ValB, ValC;
  a=0.00001;
  b=MaxMultiplier(InputVect, TheGraph->nbVert);
  ScalarMult(TheGraph, InputVect, &(TheGraph->VectA), a);
  ValA=NormalizationValue(TheGraph, &(TheGraph->VectA));
  ScalarMult(TheGraph, InputVect, &(TheGraph->VectB), b);
  ValB=NormalizationValue(TheGraph, &(TheGraph->VectB));
  //  fprintf(stderr, "ValA=%f, ValB=%f\n", ValA, ValB);
  while(1)
    {
      c=(a+b)/2;
      ScalarMult(TheGraph, InputVect, &(TheGraph->VectC), c);
      ValC=NormalizationValue(TheGraph, &(TheGraph->VectC));
      //      fprintf(stderr, "a=%f  b=%f\n", a, b);
      //      fprintf(stderr, "ValA=%f ValB=%f\n", ValA, ValB);
      if (fabs(a-b) < 0.000000001)
	break;
      if (ValC == 0)
	break;
      if (ValC<0)
	{
	  if (ValB >0)
	    {
	      a=c;
	      ValA=ValC;
	    }
	  else
	    {
	      b=c;
	      ValB=ValC;
	    }
	}
      else
	{
	  if (ValB < 0)
	    {
	      a=c;
	      ValA=ValC;
	    }
	  else
	    {
	      b=c;
	      ValB=ValC;
	    }
	}
    }
  ScalarMult(TheGraph, InputVect, OutputVect, c);
}

void MoharInitial(PlanGraph *TheGraph, double **OutputVect)
{
  int iVert;
  for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
    {
      TheGraph->MyAttempt[iVert-1]=1;
    }
  //  fprintf(stderr, "Before call to SimpleRenormalization\n");
  SimpleRenormalization(TheGraph, &(TheGraph->MyAttempt), OutputVect);
  //  fprintf(stderr, " After call to SimpleRenormalization\n");
}






void StepMoharAlgorithm(PlanGraph *TheGraph, double **TheRadiusOld, double **OutputVect)
{
  double a;
  double FuncA, FuncB;
  unsigned int iVert;
  int i, NbIteration;
  double TheDiff, TheFact, TheReduce;
  DefectEvaluation(TheGraph, TheRadiusOld, &(TheGraph->TheDefect));
  for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
    {
      if (TheGraph->TheDefect[iVert-1] >= 0)
	TheGraph->SetPlus[iVert-1]=0;
      else
	TheGraph->SetPlus[iVert-1]=1;
    }
  FuncB=SquareFunctional(TheGraph, TheRadiusOld);
  //  fprintf(stderr, " Func Initial=%f\n", FuncB);

  a=MaxMultiplier(TheRadiusOld, TheGraph->nbVert)-0.0001;
  PartialMultVector(TheRadiusOld, &(TheGraph->SetPlus), &(TheGraph->TheRadiusA_Step1), a, TheGraph->nbVert);
  SimpleRenormalization(TheGraph, &(TheGraph->TheRadiusA_Step1), &(TheGraph->TheRadiusA_Step2));
  FuncA=SquareFunctional(TheGraph, &(TheGraph->TheRadiusA_Step2));

  TheReduce=0.66;
  TheDiff=0.5;
  NbIteration=20;
  while(1)
    {
      //      fprintf(stderr, "We are passing here 1\n");
      TheFact=1-TheDiff;
      for (i=1; i<=NbIteration; i++)
	{
	  //	  fprintf(stderr, "i=%d\n", i);
	  if (FuncA < FuncB)
	    {
	      //	      fprintf(stderr, "a=%f \n", a);
	      for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
		(*OutputVect)[iVert-1]=TheGraph->TheRadiusA_Step2[iVert-1];
	      return;
	    }
	  a=1+(a-1)*0.70;
	  PartialMultVector(TheRadiusOld, &(TheGraph->SetPlus), &(TheGraph->TheRadiusA_Step1), a, TheGraph->nbVert);
	  //	  fprintf(stderr, "Before SimpleRenormalization\n");
	  SimpleRenormalization(TheGraph, &(TheGraph->TheRadiusA_Step1), &(TheGraph->TheRadiusA_Step2));
	  //	  fprintf(stderr, "After SimpleRenormalization\n");
	  FuncA=SquareFunctional(TheGraph, &(TheGraph->TheRadiusA_Step2));
	}
      TheDiff=TheDiff*TheReduce;
      NbIteration=NbIteration*2;
    }
}






void MoharAlgorithm(PlanGraph *TheGraph, double tolerance)
{
  int maxite, TheStep, iVert;
  double Norm1, Norm2, diff;
  MoharInitial(TheGraph, &(TheGraph->MyAttempt));
  maxite=20;
  TheStep=0;
  while(1)
    {
      TheStep=TheStep+1;
      StepMoharAlgorithm(TheGraph, &(TheGraph->MyAttempt), &(TheGraph->TheResult));
      Norm1=SquareFunctional(TheGraph, &(TheGraph->MyAttempt));
      Norm2=SquareFunctional(TheGraph, &(TheGraph->TheResult));
      fprintf(stderr, "%d : before %f after %f \n", TheStep, Norm1, Norm2);
      if (Norm2 < tolerance)
	{
	  fprintf(stderr, "success\n");
	  for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
	    TheGraph->RadiusVector[iVert-1]=TheGraph->TheResult[iVert-1];
	  return;
	}
      diff=fabs(Norm2-Norm1);
      for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
	TheGraph->MyAttempt[iVert-1]=TheGraph->TheResult[iVert-1];
    }
}



void PrintCircleOfOneBipartition(PlanGraph *TheGraph, int BipartVal, FILE *OUTPUT)
{
  int idx;
  unsigned int iVert;
  idx=0;
  for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
    {
      if (TheGraph->BipartitionVector[iVert-1] == BipartVal)
	{
	  idx=idx+1;
	  fprintf(OUTPUT, "%d,%f,%f,%f,%f\n", idx, TheGraph->ListCoord[iVert-1][0], TheGraph->ListCoord[iVert-1][1], TheGraph->ListCoord[iVert-1][2], TheGraph->RadiusVector[iVert-1]);
	}
    }
}


void PrintPositionOfOnePart(PlanGraph *TheGraph, int BipartVal, FILE *OUTPUT)
{
  int iVert;
  for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
    if (TheGraph->BipartitionVector[iVert-1] == BipartVal)
      fprintf(OUTPUT, "%d : %f,%f,%f\n", iVert, TheGraph->ListCoord[iVert-1][0], TheGraph->ListCoord[iVert-1][1], TheGraph->ListCoord[iVert-1][2]);
}



void PrintRadiusOfOnePart(PlanGraph *TheGraph, int BipartVal, FILE *OUTPUT)
{
  int iVert;
  for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
    if (TheGraph->BipartitionVector[iVert-1] == BipartVal)
      fprintf(OUTPUT, "%d rad %f\n", iVert, TheGraph->RadiusVector[iVert-1]);
}







void PrintEdges(PlanGraph *TheGraph, FILE *OUTPUT)
{
  int iAdj;
  int iVert, TheAdj;
  for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
    for (iAdj=1; iAdj<=TheGraph->ListNbAdj[iVert-1]; iAdj++)
      {
	TheAdj=TheGraph->ListAdjacencies[iVert-1][iAdj-1];
	if (TheAdj > iVert)
	  fprintf(OUTPUT, "%d %d\n", iVert, TheAdj);
      }
}




void PrintEdgesOfOnePart(PlanGraph *TheGraph, int BipartVal, FILE *OUTPUT)
{
  int iAdj, jAdj;
  int jVert;
  int iVert, lVert;
  int SoughtJAdj, NextJAdj;
  int TheJDeg;
  SoughtJAdj=-446;
  for (iVert=1; iVert<=TheGraph->nbVert; iVert++)
    if (TheGraph->BipartitionVector[iVert-1] == BipartVal)
      {
	for (iAdj=1; iAdj<=TheGraph->ListNbAdj[iVert-1]; iAdj++)
	  {
	    jVert=TheGraph->ListAdjacencies[iVert-1][iAdj-1];
	    TheJDeg=TheGraph->ListNbAdj[jVert-1];
	    for (jAdj=1; jAdj<=TheJDeg; jAdj++)
	      if (TheGraph->ListAdjacencies[jVert-1][jAdj-1] == iVert)
		SoughtJAdj=jAdj;
	    NextJAdj=NextIdx(TheJDeg, SoughtJAdj);
	    lVert=TheGraph->ListAdjacencies[jVert-1][NextJAdj-1];
	    if (lVert > iVert)
	      fprintf(OUTPUT, "%d %d\n", iVert, lVert);
	  }
      }
}





int main(int argc, char *argv[])
{
  FILE *OUTPUT;
  PlanGraph *TheGraph;
  double MaximalDiscrepancy;
  double tolerance;
  int iarg;
  int vert1, vert2, radius1, radius2, edge1, edge2, edge12;
  vert1=0;
  vert2=0;
  radius1=0;
  radius2=0;
  edge1=0;
  edge2=0;
  edge12=0;
  if ((TheGraph = (PlanGraph*)malloc(sizeof(PlanGraph))) == 0)
    exit (EXIT_FAILURE);
  if (argc <= 2)
    {
      fprintf(stderr, "Number of arguments is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "MoharSphere TheInput TheOutput [vert1] [vert2]\n");
      fprintf(stderr, "   [radius1] [radius2] [edge1] [edge2] [edge12]\n");
      fprintf(stderr, "TheInput: the graph being considered in format:\n");
      fprintf(stderr, "nbV\n");
      fprintf(stderr, "iPart1 nbAdj1 Adj[1] .... Adj[nbAdj1]\n");
      fprintf(stderr, "...\n");
      fprintf(stderr, "...\n");
      fprintf(stderr, "TheOutput: output of computed informations\n");
      fprintf(stderr, "  [vert1]: triggers output of xyz coordinates of vertices\n");
      fprintf(stderr, "           in part 1 of the bipartition\n");
      fprintf(stderr, "  [vert2]: triggers output of xyz coordinates of vertices\n");
      fprintf(stderr, "           in part 2 of the bipartition\n");
      fprintf(stderr, "[radius1]: triggers output of radius of vertices\n");
      fprintf(stderr, "           in part 1 of the bipartition\n");
      fprintf(stderr, "[radius2]: triggers output of radius of vertices\n");
      fprintf(stderr, "           in part 2 of the bipartition\n");
      fprintf(stderr, "  [edge1]: triggers output of edges between vertices\n");
      fprintf(stderr, "           of part 1\n");
      fprintf(stderr, "  [edge2]: triggers output of edges between vertices\n");
      fprintf(stderr, "           of part 1\n");
      fprintf(stderr, " [edge12]: triggers output of edges between vertices\n");
      fprintf(stderr, "           of part 1 and 2\n");
      return -1;
    }
  for (iarg=4; iarg<=argc; iarg++)
    {
      if (strcmp(argv[iarg-1],(char *)"vert1")==0) 
	vert1=1;
      else if (strcmp(argv[iarg-1],(char *)"vert2")==0) 
	vert2=1;
      else if (strcmp(argv[iarg-1],(char *)"radius1")==0) 
	radius1=1;
      else if (strcmp(argv[iarg-1],(char *)"radius2")==0) 
	radius2=1;
      else if (strcmp(argv[iarg-1],(char *)"edge1")==0) 
	edge1=1;
      else if (strcmp(argv[iarg-1],(char *)"edge2")==0) 
	edge2=1;
      else if (strcmp(argv[iarg-1],(char *)"edge12")==0) 
	edge12=1;
      else
	{
	  fprintf(stderr, " Error unrecognized option\n");
	  return -1;
	}
    }
  ReadGraph(argv[1], TheGraph);
  tolerance=0.00001;

  MoharAlgorithm(TheGraph, tolerance);
  ComputeCoordinates(TheGraph, &MaximalDiscrepancy);
  fprintf(stderr, "MaximalDiscrepancy = %f \n", MaximalDiscrepancy);
  OUTPUT=fopen(argv[2], "w");
  if (radius1)
    PrintRadiusOfOnePart(TheGraph, 1, OUTPUT);
  if (radius2)
    PrintRadiusOfOnePart(TheGraph, 2, OUTPUT);
  if (vert1)
    PrintPositionOfOnePart(TheGraph, 1, OUTPUT);
  if (vert2)
    PrintPositionOfOnePart(TheGraph, 2, OUTPUT);
  if (edge1)
    PrintEdgesOfOnePart(TheGraph, 1, OUTPUT);
  if (edge2)
    PrintEdgesOfOnePart(TheGraph, 2, OUTPUT);
  if (edge12)
    PrintEdges(TheGraph, OUTPUT);
  fclose(OUTPUT);
  return 0;
}
