#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <string.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <math.h>

using namespace std;


void ReadGraph(char* TheFile, vector<vector<unsigned int> >& TheGraph, vector<int> & BipartitionVector)
{
  ifstream* pFileStream = 0;
  long nbVert, iVert, i, nbAdj;
  vector<unsigned int> ListAdj;
  int ThePart;
  int TheVal;
  pFileStream=new ifstream(TheFile);
  (*pFileStream) >> nbVert;
  for (iVert=1; iVert<=nbVert; iVert++)
    {
      (*pFileStream) >> ThePart;
      BipartitionVector.push_back(ThePart);
      (*pFileStream) >> nbAdj;
      for (i=1; i<=nbAdj; i++)
	{
	  (*pFileStream) >> TheVal;
	  ListAdj.push_back(TheVal);
	}
      TheGraph.push_back(ListAdj);
      ListAdj.clear();
    }
  delete pFileStream;
}




double Alpha(double& Ru, double& Rv)
{
  return atan(tan(Ru)/sin(Rv));
}


void difAlpha(double& Ru, double& Rv, double& difRu, double& difRv)
{
  double W, s1, s2;
  W=tan(Ru)/sin(Rv);
  s1=cos(Ru);
  s2=sin(Rv);
  difRu=1/(s2*s1*s1*(1+W*W));
  difRv=-tan(Ru)*cos(Rv)/(s2*s2*(1+W*W));
}



vector<double> DefectEvaluation(vector<vector<unsigned int> >& TheGraph, vector<double>& TheRadius)
{
  double defect, pi, angle;
  long uVert, vVert, TheSize;
  vector<unsigned int> ListAdj;
  vector<unsigned int>::iterator iter;
  vector<double> TheDefect;
  TheSize=TheGraph.size();
  pi=3.141592653589792;
  for (vVert=1; vVert<=TheSize; vVert++)
    {
      defect=pi;
      ListAdj=TheGraph[vVert-1];
      iter=ListAdj.begin();
      while(iter != ListAdj.end())
	{
	  uVert=*iter;
	  angle=Alpha(TheRadius[uVert-1], TheRadius[vVert-1]);
	  //	  fprintf(stderr, "angle=%f\n", angle);
	  defect=defect-angle;
	  iter++;
	}
      TheDefect.push_back(defect);
      //      cout << "defect=" << defect << "\n";
    }
  return TheDefect;
}


double SquareFunctional(vector<vector<unsigned int> >& TheGraph, vector<double>& TheRadius)
{
  double sqr, defect;
  unsigned int iVert;
  vector<double> TheDefect;
  sqr=0;
  TheDefect=DefectEvaluation(TheGraph, TheRadius);
  for (iVert=1; iVert<=TheGraph.size(); iVert++)
    {
      defect=TheDefect[iVert-1];
      sqr=sqr+defect*defect;
    }
  return sqr;
}



vector<double> ScalarMultVector(vector<double> & InputVect, double & val)
{
  vector<double> w;
  unsigned int i;
  for (i=1; i<=InputVect.size(); i++)
    {
      w.push_back(val*InputVect[i-1]);
    }
  return w;
}

vector<double> AdditionVector(vector<double> & V1, vector<double> & V2)
{
  vector<double> w;
  unsigned int i;
  for (i=1; i<=V1.size(); i++)
    {
      w.push_back(V1[i-1]+V2[i-1]);
    }
  return w;
}

vector<double> VectorProduct(vector<double> & a, vector<double> & b)
{
  vector<double> c;
  double c1, c2, c3;
  c1=a[1]*b[2]-a[2]*b[1];
  c2=a[2]*b[0]-a[0]*b[2];
  c3=a[0]*b[1]-a[1]*b[0];
  c.push_back(c1);
  c.push_back(c2);
  c.push_back(c3);
  return c;
}


vector<double> RenormalizeVector(vector<double> & TheVect)
{
  double TheNorm;
  int iCol;
  TheNorm=0;
  vector<double> RetVect;
  for (iCol=1; iCol<=3; iCol++)
    {
      TheNorm=TheNorm+TheVect[iCol-1]*TheVect[iCol-1];
    }
  for (iCol=1; iCol<=3; iCol++)
    {
      RetVect.push_back(TheVect[iCol-1]/sqrt(TheNorm));
    }
  return RetVect;
}




double VectorScalarProduct(vector<double> & V1, vector<double> & V2)
{
  double scal;
  unsigned int i;
  scal=0;
  for (i=1; i<=V1.size(); i++)
    {
      scal=scal+V1[i-1]*V2[i-1];
    }
  return scal;
}

double Determinant3times3matrix(vector<vector<double> > & TheMat)
{
  vector<double> Line1, Line2, Line3;
  double TheVal;
  TheVal=0;
  Line1=TheMat[0];
  Line2=TheMat[1];
  Line3=TheMat[2];
  TheVal=TheVal+TheMat[0][0]*TheMat[1][1]*TheMat[2][2];
  TheVal=TheVal+TheMat[1][0]*TheMat[2][1]*TheMat[0][2];
  TheVal=TheVal+TheMat[2][0]*TheMat[0][1]*TheMat[1][2];

  TheVal=TheVal-TheMat[2][0]*TheMat[1][1]*TheMat[0][2];
  TheVal=TheVal-TheMat[0][0]*TheMat[2][1]*TheMat[1][2];
  TheVal=TheVal-TheMat[1][0]*TheMat[0][1]*TheMat[2][2];
  return TheVal;
}



void DifferentialSquareFunctional(vector<vector<unsigned int> >& TheGraph, vector<double>& TheRadius, vector<double> & Differential)
{
  double difRu, difRv, zero;
  long uVert, vVert, TheSize, iVert;
  vector<unsigned int> ListAdj;
  vector<unsigned int>::iterator iter;
  vector<double> TheDefect;
  Differential.clear();
  TheSize=TheRadius.size();
  for (iVert=1; iVert<=TheSize; iVert++)
    {
      zero=0;
      Differential.push_back(zero);
    }
  TheDefect=DefectEvaluation(TheGraph, TheRadius);
  for (vVert=1; vVert<=TheSize; vVert++)
    {
      ListAdj=TheGraph[vVert-1];
      iter=ListAdj.begin();
      while(iter != ListAdj.end())
	{
	  uVert=*iter;
	  difAlpha(TheRadius[uVert-1], TheRadius[vVert-1], difRu, difRv);
	  Differential[uVert-1]=Differential[uVert-1]-TheDefect[vVert-1]*difRu;
	  Differential[vVert-1]=Differential[vVert-1]-TheDefect[vVert-1]*difRv;
	  iter++;
	}
    }
}










vector<double> SteepestGradient(vector<vector<unsigned int> >& TheGraph, vector<double> & TheRadiusOld, int & maxite)
{
  double deltaLow, deltaHigh, deltaMiddle, NormLow, NormMiddle, NormHigh, deltaWork, NormWork, Norm, NewNorm;
  int nbite, nbmove, i;
  unsigned int iVert;
  int nbSubDiv;
  vector<double> Differential, Vect, NewVect, VectLow, VectMiddle, VectHigh, VectWork;
  DifferentialSquareFunctional(TheGraph, TheRadiusOld, Differential);
  deltaLow=0;
  deltaHigh=100;

  nbite=0;
  while(nbite<maxite)
    {
      Vect=TheRadiusOld;
      nbmove=0;
      while(true)
	{
	  Norm=SquareFunctional(TheGraph, Vect);
	  NewVect.clear();
	  for (iVert=1; iVert<=TheRadiusOld.size(); iVert++)
	    {
	      NewVect.push_back(Vect[iVert-1]-deltaHigh*Differential[iVert-1]);
	    }
	  NewNorm=SquareFunctional(TheGraph, NewVect);
	  if (NewNorm>Norm)
	    {
	      break;
	    }
	  Vect=NewVect;
	  nbmove=nbmove+1;
	}
      if (nbmove>0)
	{
	  break;
	}
      deltaHigh=deltaHigh/2;
      nbite=nbite+1;
    }
  deltaHigh=deltaHigh*2;
  nbite=0;
  nbSubDiv=100;
  while(nbite<maxite)
    {
      nbmove=0;
      VectLow.clear();
      for (iVert=1; iVert<=TheRadiusOld.size(); iVert++)
	{
	  VectLow.push_back(TheRadiusOld[iVert-1]-deltaLow*Differential[iVert-1]);
	}
      NormLow=SquareFunctional(TheGraph, VectLow);

      deltaMiddle=deltaLow;
      VectMiddle=VectLow;
      NormMiddle=NormLow;
      for (i=1; i<=nbSubDiv; i++)
	{
	  deltaWork=deltaLow+(deltaHigh-deltaLow)*(i/nbSubDiv);
	  VectWork.clear();
	  for (iVert=1; iVert<=TheRadiusOld.size(); iVert++)
	    {
	      VectWork.push_back(TheRadiusOld[iVert-1]-deltaWork*Differential[iVert-1]);
	    }
	  NormWork=SquareFunctional(TheGraph, VectWork);
	  if (NormWork < NormMiddle)
	    {
	      deltaMiddle=deltaWork;
	      VectMiddle=VectWork;
	      NormMiddle=NormWork;
	    }
	}
      VectHigh.clear();
      for (iVert=1; iVert<=TheRadiusOld.size(); iVert++)
	{
	  VectHigh.push_back(TheRadiusOld[iVert-1]-deltaHigh*Differential[iVert-1]);
	}
      NormHigh=SquareFunctional(TheGraph, VectHigh);
      cerr << "deltaLow=" << deltaLow << "   NormLow=" << NormLow << "\n";
      cerr << "deltaMiddle=" << deltaMiddle << "   NormMiddle=" << NormMiddle << "\n";
      cerr << "deltaHigh=" << deltaHigh << "   NormHigh=" << NormHigh << "\n\n";

      if (NormMiddle < NormLow && NormLow < NormHigh)
	{
	  deltaHigh=deltaMiddle;
	}
      if (NormMiddle < NormHigh && NormHigh < NormLow)
	{
	  deltaLow=deltaMiddle;
	}
      if (NormLow < NormMiddle && NormMiddle < NormHigh)
	{
	  deltaHigh=deltaMiddle;
	}
      if (NormHigh < NormMiddle && NormMiddle < NormLow)
	{
	  deltaLow=deltaMiddle;
	}
      if (NormHigh < NormMiddle && NormLow < NormMiddle)
	{
	  break;
	}
      nbite=nbite+1;
    }
  return VectMiddle;
}






double FuncDistance(double & Ru, double & Rv)
{
  return acos(cos(Ru)*cos(Rv));
}




int NextIdx(int & TheDeg, int & idx)
{
  if (idx < TheDeg)
    {
      return idx+1;
    }
  return 1;
}



void PrintVector(vector<double> & TheVect)
{
  unsigned int i;
  cerr << "(";
  for (i=1; i<=TheVect.size(); i++)
    {
      cerr << TheVect[i-1];
      if (i < TheVect.size())
	{
	  cerr << ",";
	}
    }
  cerr << ")\n";
}




void VectorsU2U3(vector<double> & u, vector<double> & v, vector<double> & u2, vector<double> & u3)
{
  vector<double> w, RSL;
  double scal, norm;
  int i;
  scal=VectorScalarProduct(u, v);
  w.clear();
  for (i=1; i<=3; i++)
    {
      w.push_back(v[i-1]-scal*u[i-1]);
    }
  norm=VectorScalarProduct(w, w);
  u2.clear();
  for (i=1; i<=3; i++)
    {
      u2.push_back(w[i-1]/sqrt(norm));
    }
  u3.clear();
  RSL=VectorProduct(u, u2);
  for (i=1; i<=3; i++)
    {
      u3.push_back(RSL[i-1]);
    }
}


/*
vector<int> Bipartition(vector<vector<unsigned int> > & TheGraph)
{
  vector<int> TheBipart, ListDone, ListStatus;
  vector<unsigned int> ListAdj;
  unsigned int iAdj;
  unsigned int iVert, i;
  int BipartVal, IsAllDone, jVert;
  for (i=1; i<=TheGraph.size(); i++)
    {
      TheBipart.push_back(0);
      ListDone.push_back(0);
      ListStatus.push_back(0);
    }
  TheBipart[0]=1;
  ListStatus[0]=1;
  while(true)
    {
      IsAllDone=1;
      for (iVert=1; iVert<=TheGraph.size(); iVert++)
	{
	  if (ListStatus[iVert-1] == 1 && ListDone[iVert-1] == 0)
	    {
	      IsAllDone=0;
	      ListDone[iVert-1]=1;
	      ListAdj=TheGraph[iVert-1];
	      for (iAdj=1; iAdj<=ListAdj.size(); iAdj++)
		{
		  jVert=ListAdj[iAdj-1];
		  BipartVal=3-TheBipart[iVert-1];
		  if (ListStatus[jVert-1] == 0)
		    {
		      TheBipart[jVert-1]=BipartVal;
		      ListStatus[jVert-1]=1;
		    }
		  else
		    {
		      if (BipartVal != TheBipart[jVert-1])
			{
			  cerr << "We have a bipartition error, PANIC !!!!\n";
			  cerr << "BipartVal=" << BipartVal << " Orig=" << TheBipart[jVert-1] << "\n";
			}
		    }
		}

	      
	    }
	}
      if (IsAllDone == 1)
	{
	  break;
	}
    }
  return TheBipart;
}
*/



void ComputeCoordinates(vector<vector<unsigned int> > & TheGraph, vector<double> & RadiusVector, vector<vector<double> > & ListCoord, double & MaximalDiscrepancy)
{
  vector<unsigned int> ListAdj;
  vector<int> ListStatus, ListDone;
  vector<double> SingleCoord, u2, u3, DirectionVect, TheDiff;
  vector<double> prov1, prov2;
  double angle, TheDist, MinusOne, val1, val2, dist;
  int jVert, iAdj, PrevAdj, TheIAdj, TheDeg, i;
  int PrevVert, TheVert;
  int IsAllDone;
  unsigned int iVert;
  PrevAdj=-1;
  for (i=1; i<=3; i++)
    SingleCoord.push_back(0);
  for (iVert=1; iVert<=TheGraph.size(); iVert++)
    {
      ListCoord.push_back(SingleCoord);
      ListStatus.push_back(0);
      ListDone.push_back(0);
    }
  ListAdj=TheGraph[0];
  iVert=1;
  jVert=ListAdj[0];
  TheDist=FuncDistance(RadiusVector[iVert-1], RadiusVector[jVert-1]);
  SingleCoord[0]=1;
  ListCoord[iVert-1]=SingleCoord;
  SingleCoord[0]=cos(TheDist);
  SingleCoord[1]=sin(TheDist);
  ListCoord[jVert-1]=SingleCoord;
  ListStatus[iVert-1]=1;
  ListStatus[jVert-1]=1;
  MaximalDiscrepancy=0;
  while(true)
    {
      IsAllDone=1;
      for (iVert=1; iVert<=TheGraph.size(); iVert++)
	{
	  if (ListStatus[iVert-1] == 1 && ListDone[iVert-1] == 0)
	    {
	      ListDone[iVert-1]=1;
	      ListAdj=TheGraph[iVert-1];
	      TheDeg=ListAdj.size();
	      for (iAdj=1; iAdj<=TheDeg; iAdj++)
		{
		  jVert=ListAdj[iAdj-1];
		  if (ListStatus[jVert-1] == 1)
		    PrevAdj=iAdj;
		}
	      VectorsU2U3(ListCoord[iVert-1], ListCoord[ListAdj[PrevAdj-1]-1], u2, u3);
	      angle=0;
	      for (iAdj=1; iAdj<=TheDeg; iAdj++)
		{
		  TheIAdj=NextIdx(TheDeg, PrevAdj);
		  PrevVert=ListAdj[PrevAdj-1];
		  TheVert=ListAdj[TheIAdj-1];
		  angle=angle+Alpha(RadiusVector[PrevVert-1], RadiusVector[iVert-1])+Alpha(RadiusVector[TheVert-1], RadiusVector[iVert-1]);
		  val1=cos(angle); prov1=ScalarMultVector(u2, val1);
		  val2=sin(angle); prov2=ScalarMultVector(u3, val2);
		  DirectionVect=AdditionVector(prov1, prov2);
		  TheDist=FuncDistance(RadiusVector[iVert-1], RadiusVector[TheVert-1]);
		  val1=cos(TheDist); prov1=ScalarMultVector(ListCoord[iVert-1], val1);
		  val2=sin(TheDist); prov2=ScalarMultVector(DirectionVect, val2);
		  SingleCoord=AdditionVector(prov1, prov2);
		  if (ListStatus[TheVert-1] == 0)
		    {
		      ListCoord[TheVert-1]=SingleCoord;
		      ListStatus[TheVert-1]=1;
		      IsAllDone=0;
		    }
		  else
		    {
		      MinusOne=-1;
		      prov1=ScalarMultVector(ListCoord[TheVert-1], MinusOne);
		      TheDiff=AdditionVector(SingleCoord, prov1);
		      
		      dist=VectorScalarProduct(TheDiff, TheDiff);
		      if (dist > MaximalDiscrepancy)
			MaximalDiscrepancy=dist;
		    }
		  PrevAdj=TheIAdj;
		}
	    }
	}
      if (IsAllDone == 1)
	break;
    }
}






vector<double> ConstantMultVector(vector<double> & TheRadiusOld, vector<int> & SetPlus, double& MultPlus, double & MultMinus)
{
  vector<double> TheResult;
  unsigned int iVert;
  TheResult.clear();
  for (iVert=1; iVert<=TheRadiusOld.size(); iVert++)
    {
      if (SetPlus[iVert-1] == 1)
	TheResult.push_back(TheRadiusOld[iVert-1]*MultPlus);
      else
	TheResult.push_back(TheRadiusOld[iVert-1]*MultMinus);
    }
  return TheResult;
}

vector<double> PartialMultVector(vector<double> & TheRadiusOld, vector<int> & SetPlus, double& MultPlus)
{
  vector<double> TheResult;
  unsigned int iVert;
  TheResult.clear();
  for (iVert=1; iVert<=TheRadiusOld.size(); iVert++)
    {
      if (SetPlus[iVert-1] == 1)
	TheResult.push_back(TheRadiusOld[iVert-1]*MultPlus);
      else
	TheResult.push_back(TheRadiusOld[iVert-1]);
    }
  return TheResult;
}








double NormalizationValue(vector<vector<unsigned int> >& TheGraph, vector<double> & TheRadius)
{
  vector<double> TheDefect;
  unsigned int iVert;
  double res;
  TheDefect=DefectEvaluation(TheGraph, TheRadius);
  res=0;
  for (iVert=1; iVert<=TheDefect.size(); iVert++)
    res=res+TheDefect[iVert-1];
  return res;
}


vector<double> ConstantVector(size_t nbV, double & val)
{
  vector<double> eVect;
  unsigned int iVert;
  for (iVert=1; iVert<=nbV; iVert++)
    {
      eVect.push_back(val);
    }
  return eVect;
}


vector<double> ScalarMult(vector<double> & InputVect, double & TheMult)
{
  vector<double> OutputVect;
  unsigned int iVert;
  for (iVert=1; iVert<=InputVect.size(); iVert++)
    OutputVect.push_back(InputVect[iVert-1]*TheMult);
  return OutputVect;
}



double MaxMultiplier(vector<double>& TheRadius)
{
  double pi, TheMin, fact;
  unsigned int iVert;
  pi=3.141592653589792;
  TheMin=(pi/2)/TheRadius[0];
  for (iVert=2; iVert<=TheRadius.size(); iVert++)
    {
      fact=(pi/2)/TheRadius[iVert-1];
      if (fact < TheMin)
	TheMin=fact;
    }
  return TheMin;
}


vector<double> SimpleRenormalization(vector<vector<unsigned int> >& TheGraph, vector<double> & InputVect)
{
  double a, b, c, ValA, ValB, ValC;
  vector<double> VectA, VectB, VectC;
  a=0.00001;
  b=MaxMultiplier(InputVect);
  //
  VectB=ScalarMult(InputVect, b);
  ValB=NormalizationValue(TheGraph, VectB);
  //
  VectA=ScalarMult(InputVect, a);
  ValA=NormalizationValue(TheGraph, VectA);
  //
  while(true)
    {
      c=(a+b)/2;
      VectC=ScalarMult(InputVect, c);
      ValC=NormalizationValue(TheGraph, VectC);
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
  return ScalarMult(InputVect, c);
}

vector<double> MoharInitial(vector<vector<unsigned int> >& TheGraph)
{
  vector<double> InputVect;
  double alpha;
  alpha=1;
  InputVect=ConstantVector(TheGraph.size(), alpha);
  return SimpleRenormalization(TheGraph, InputVect);
}






vector<double> StepMoharAlgorithm(vector<vector<unsigned int> >& TheGraph, vector<double> & TheRadiusOld)
{
  vector<double> TheDefect;
  vector<int> SetPlus, SetMinus;
  vector<double> TheRadiusA_Step1, TheRadiusA_Step2;
  double a;
  double FuncA, FuncB;
  unsigned int iVert;
  int i, NbIteration;
  double TheDiff, TheReduce;
  TheDefect=DefectEvaluation(TheGraph, TheRadiusOld);
  for (iVert=1; iVert<=TheRadiusOld.size(); iVert++)
    {
      if (TheDefect[iVert-1] >= 0)
	{
	  SetPlus.push_back(0);
	}
      else
	{
	  SetPlus.push_back(1);
	}
    }
  FuncB=SquareFunctional(TheGraph, TheRadiusOld);
  cerr << " Func Initial=" << FuncB << "\n";

  a=MaxMultiplier(TheRadiusOld)-0.0001;
  TheRadiusA_Step1=PartialMultVector(TheRadiusOld, SetPlus, a);
  TheRadiusA_Step2=SimpleRenormalization(TheGraph, TheRadiusA_Step1);
  FuncA=SquareFunctional(TheGraph, TheRadiusA_Step2);

  TheReduce=0.66;
  TheDiff=0.5;
  NbIteration=20;
  while(true)
    {
      for (i=1; i<=NbIteration; i++)
	{
	  if (FuncA < FuncB)
	    {
	      cerr << "a=" << a << "\n";
	      return TheRadiusA_Step2;
	    }
	  a=1+(a-1)*0.80;
	  TheRadiusA_Step1=PartialMultVector(TheRadiusOld, SetPlus, a);
	  TheRadiusA_Step2=SimpleRenormalization(TheGraph, TheRadiusA_Step1);
	  FuncA=SquareFunctional(TheGraph, TheRadiusA_Step2);
	}
      TheDiff=TheDiff*TheReduce;
      NbIteration=NbIteration*2;
    }
}






vector<double> MoharAlgorithm(vector<vector<unsigned int> > & TheGraph, double & tolerance)
{
  vector<double> MyAttempt, TheResult;
  int TheStep;
  double Norm1, Norm2;
  MyAttempt=MoharInitial(TheGraph);
  TheStep=0;
  while(true)
    {
      TheStep=TheStep+1;
      TheResult=StepMoharAlgorithm(TheGraph, MyAttempt);
      Norm1=SquareFunctional(TheGraph, MyAttempt);
      Norm2=SquareFunctional(TheGraph, TheResult);
      cerr << TheStep << " : before " << Norm1 << " after " << Norm2 << "\n";
      if (Norm2 < tolerance)
	{
	  cerr << "success\n";
	  break;
	}
      /*
      double diff=fabs(Norm2-Norm1);
      if (diff < 0.000000001)
	{
	  TheResult=SteepestGradient(TheGraph, MyAttempt, maxite);
	  Norm2=SquareFunctional(TheGraph, TheResult);
	  diff=fabs(Norm2-Norm1);
	  if (diff < 0.000000001)
	    {
	      cerr << "insufficient progress\n";
	      break;
	    }
	}
      */
      MyAttempt=TheResult;
    }
  return TheResult;
}



void PrintCircleOfOneBipartition(vector<int> & TheBipartition, vector<vector<double> > & ListCoord, vector<double> & RadiusVector, int & BipartVal, ofstream & oStream)
{
  int idx;
  unsigned int iVert;
  vector<double> TheVect;
  idx=0;
  for (iVert=1; iVert<=TheBipartition.size(); iVert++)
    if (TheBipartition[iVert-1] == BipartVal)
      {
	idx=idx+1;
	TheVect=ListCoord[iVert-1];
	oStream << idx << "," << TheVect[0] << "," << TheVect[1] << "," << TheVect[2] << "," << RadiusVector[iVert-1] << "\n";
      }
}


void PrintPositionOfOnePart(vector<int> & TheBipartition, vector<vector<double> > & ListCoord, int BipartVal, ofstream & oStream)
{
  unsigned int iVert;
  vector<double> TheVect;
  for (iVert=1; iVert<=TheBipartition.size(); iVert++)
    if (TheBipartition[iVert-1] == BipartVal)
      {
	TheVect=ListCoord[iVert-1];
	oStream << iVert << " : " << TheVect[0] << "," << TheVect[1] << "," << TheVect[2] << "\n";
      }
}

void PrintGanglInputFile(vector<vector<unsigned int> > & TheGraph, vector<int> & TheBipartition, std::vector<std::vector<double>> const& ListCoord, std::ostream & os)
{
  int idxVal=1;
  bool IsFirstV=true;
  int nbVert=TheBipartition.size();
  int nbVert1=0;
  for (int iVert=0; iVert<nbVert; iVert++)
    if (TheBipartition[iVert] == idxVal)
      nbVert1++;
  std::vector<int> MapDirect(nbVert,-1);
  std::vector<int> MapRev(nbVert1, -1);
  int idx1=0;
  for (int iVert=0; iVert<nbVert; iVert++) {
    if (TheBipartition[iVert] == idxVal) {
      MapDirect[iVert]=idx1;
      MapRev[idx1]=iVert;
      idx1++;
    }
  }
  os << "points=[";
  for (int iVert=0; iVert<nbVert; iVert++) {
    if (TheBipartition[iVert] == idxVal) {
      if (IsFirstV == false)
	os << ",";
      IsFirstV=false;
      os << "[";
      for (int i=0; i<3; i++) {
	if (i>0)
	  os << ",";
	os << ListCoord[iVert][i];
      }
      os << "]";
    }
  }
  os << "];\n";
  bool IsFirstE=true;
  os << "edges=[";
  for (unsigned int iVert=0; iVert<=size_t(nbVert); iVert++) {
    if (TheBipartition[iVert] == idxVal) {
      std::vector<unsigned int> ListIAdj=TheGraph[iVert];
      for (int iAdj=0; iAdj<int(ListIAdj.size()); iAdj++)
	{
	  int jVert=ListIAdj[iAdj];
	  std::vector<unsigned int> ListJAdj=TheGraph[jVert-1];
	  int TheJDeg=ListJAdj.size();
	  int SoughtJAdj=-1;
	  for (int jAdj=1; jAdj<=TheJDeg; jAdj++)
	    if (ListJAdj[jAdj-1] == iVert+1)
	      SoughtJAdj=jAdj;
	  std::cerr << "SoughtJAdj=" << SoughtJAdj << "\n";
	  unsigned int NextJAdj=NextIdx(TheJDeg, SoughtJAdj);
	  std::cerr << "NextJAdj=" << NextJAdj << "\n";
	  unsigned int lVert=ListJAdj[NextJAdj-1]-1;
	  if (lVert > iVert) {
	    int iVertMap=MapDirect[iVert];
	    int lVertMap=MapDirect[lVert];
	    if (IsFirstE == false) {
	      os << ",";
	    }
	    IsFirstE=false;
	    os << "[" << iVertMap << "," << lVertMap << "]";
	  }
	}
    }
  }
  os << "];\n";
  

  
}





void PrintRadiusOfOnePart(vector<int> & TheBipartition, vector<double> & RadiusVector, int BipartVal, ofstream & oStream)
{
  unsigned int iVert;
  vector<double> TheVect;
  for (iVert=1; iVert<=TheBipartition.size(); iVert++)
    if (TheBipartition[iVert-1] == BipartVal)
      oStream << iVert << " rad " << RadiusVector[iVert-1] << "\n";
}





set<int> VectorIntToSetInt(vector<int> & TheVect)
{
  set<int> TheReply;
  unsigned int iVert;
  for (iVert=1; iVert<=TheVect.size(); iVert++)
    {
      TheReply.insert(TheVect[iVert-1]);
    }
  return TheReply;
}


void PrintEdges(vector<vector<unsigned int> > & TheGraph, ofstream & oStream)
{
  unsigned int iAdj;
  unsigned int iVert, TheAdj;
  vector<unsigned int> ListIAdj;
  for (iVert=1; iVert<=TheGraph.size(); iVert++)
    {
      ListIAdj=TheGraph[iVert-1];
      for (iAdj=1; iAdj<=ListIAdj.size(); iAdj++)
	{
	  TheAdj=ListIAdj[iAdj-1];
	  if (TheAdj > iVert)
	    oStream << iVert << " " << TheAdj << "\n";
	}
    }
}




void PrintEdgesOfOnePart(vector<vector<unsigned int> > & TheGraph, vector<int> & TheBipartition, int BipartVal, ofstream & oStream)
{
  vector<double> iVect, lVect;
  unsigned int iAdj, jAdj;
  int jVert;
  unsigned int iVert, lVert;
  vector<unsigned int> ListIAdj, ListJAdj;
  int SoughtJAdj, NextJAdj;
  int TheJDeg;
  SoughtJAdj=-1;
  for (iVert=1; iVert<=TheBipartition.size(); iVert++)
    {
      if (TheBipartition[iVert-1] == BipartVal)
	{
	  ListIAdj=TheGraph[iVert-1];
	  for (iAdj=1; iAdj<=ListIAdj.size(); iAdj++)
	    {
	      jVert=ListIAdj[iAdj-1];
	      ListJAdj=TheGraph[jVert-1];
	      TheJDeg=ListJAdj.size();
	      for (jAdj=1; jAdj<=ListJAdj.size(); jAdj++)
		if (ListJAdj[jAdj-1] == iVert)
		  SoughtJAdj=jAdj;
	      NextJAdj=NextIdx(TheJDeg, SoughtJAdj);
	      lVert=ListJAdj[NextJAdj-1];
	      if (lVert > iVert)
		oStream << iVert << " " << lVert << "\n";
	    }
	}
    }
}





int main(int argc, char *argv[])
{
  vector<vector<unsigned int> > TheGraph;
  vector<double> RadiusVector;
  vector<vector<double> > ListCoord;
  double MaximalDiscrepancy;
  double tolerance;
  vector<int> TheBipartition;
  ofstream oStream;
  int iarg;
  int vert1, vert2, radius1, radius2, edge1, edge2, edge12;
  int gangl;
  vert1=0;
  vert2=0;
  radius1=0;
  radius2=0;
  edge1=0;
  edge2=0;
  edge12=0;
  gangl=0;
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
      else if (strcmp(argv[iarg-1],(char *)"gangl")==0) 
	gangl=1;
      else
	{
	  cerr << " Error unrecognized option\n";
	  return -1;
	}
    }

  ReadGraph(argv[1], TheGraph, TheBipartition);
  tolerance=0.00001;
  
  RadiusVector=MoharAlgorithm(TheGraph, tolerance);
  ComputeCoordinates(TheGraph, RadiusVector, ListCoord, MaximalDiscrepancy);
  cerr << "MaximalDiscrepancy = " << MaximalDiscrepancy << "\n";
  // TheBipartition=Bipartition(TheGraph);



  oStream.open(argv[2]);
  if (radius1)
    PrintRadiusOfOnePart(TheBipartition, RadiusVector, 1, oStream);
  if (radius2)
    PrintRadiusOfOnePart(TheBipartition, RadiusVector, 2, oStream);
  if (vert1)
    PrintPositionOfOnePart(TheBipartition, ListCoord, 1, oStream);
  if (vert2)
    PrintPositionOfOnePart(TheBipartition, ListCoord, 2, oStream);
  if (edge1)
    PrintEdgesOfOnePart(TheGraph, TheBipartition, 1, oStream);
  if (edge2)
    PrintEdgesOfOnePart(TheGraph, TheBipartition, 2, oStream);
  if (gangl)
    PrintGanglInputFile(TheGraph, TheBipartition, ListCoord, oStream);
  if (edge12)
    PrintEdges(TheGraph, oStream);
  oStream.close();
  std::cerr << "gangl=" << gangl << "\n";
  std::cerr << "radius1=" << radius1 << " radius2=" << radius2 << "\n";
  std::cerr << "vert1=" << vert1 << " vert2=" << vert2 << "\n";
  std::cerr << "edge1=" << edge1 << " edge2=" << edge2 << "\n";
  std::cerr << "edge12=" << edge12 << "\n";
  std::cerr << "gangl=" << gangl << "\n";
  return 0;
}
