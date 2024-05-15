#include "TestModel.h"



void PrintPopNew(OPTIONS *Opt, int Itter, POP *Pop, double *Vect)
{
	int Index;


	for(Index=0;Index<Opt->PopSize;Index++)
		Vect[Index] = CaclGenotype(Opt, Pop->Pop[Index]);

//		printf("%f\t", CaclGenotype(Opt, Pop->Pop[Index]));
//		printf("%f\t", Pop->Pop[Index]->Genome[0][0]+Pop->Pop[Index]->Genome[1][0]);

	printf("%d\t%f\t", Itter, gsl_stats_variance(Vect, 1, Opt->PopSize));


	printf("\n");
}

void TestModel(OPTIONS *Opt)
{
	int Itter;
	double *Vect;
	POP *CPop, *NPop;


	Vect = (double*)SMalloc(sizeof(double) * Opt->PopSize);
	CPop = CreatePop(Opt);
	NPop = CreatePop(Opt);

/*	for(Itter=0;Itter<Opt->NoGenerations;Itter++)
	{
		PrintPopNew(Opt, Itter, CPop, Vect);
		MutatePop(Opt, CPop);
	}
*/

	Opt->FitnessPar[1] = 4.0;

	for(Itter=0;Itter<Opt->NoGenerations;Itter++)
	{
		PrintPopNew(Opt, Itter, CPop, Vect);

//		PopFitness(Opt, Itter, CPop);

//		NewPop(Opt, CPop, NPop);
				
		MutatePop(Opt, CPop);

//		swap((void**)&CPop, (void**)&NPop);
	}

	free(Vect);
	exit(0);

}