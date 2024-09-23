/*
*  QTLSim 1.0
*
*  copyright 2023
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* GenomeSim is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GenLib.h"
#include "PassOptionFile.h"

void TestNoArg(char *Cmd, int argc, int Expect)
{
	if(argc != Expect)
	{
		printf("%s expecting %d args but found %d\n", Cmd, Expect, argc);
		exit(1);
	}
}

void TestInt(char *Cmd, char* Str)
{
	if(IsValidInt(Str) == FALSE)
	{
		printf("command %s: %s is not a valid int", Cmd, Str);
		exit(1);
	}
}

void TestDouble(char *Cmd, char* Str)
{
	if(IsValidDouble(Str) == FALSE)
	{
		printf("command %s: %s is not a valid float", Cmd, Str);
		exit(1);
	}
}
void Help(OPTIONS* Opt, int argc, char** argv)
{
	int Index;

	Index = 0;
	while(CMD_LIST[Index].Cmd != NULL)
	{
		printf("%s\n", CMD_LIST[Index].Cmd);
		Index++;
	}

}

void SetSampleNo(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestInt(argv[0], argv[1]);
	Opt->Sample = atoi(argv[1]);
}

void SetPopSize(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestInt(argv[0], argv[1]);
	Opt->PopSize = atoi(argv[1]);
}

void SetNoQTL(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestInt(argv[0], argv[1]);
	Opt->NoQTL = atoi(argv[1]);
}

void SetPloidy(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestInt(argv[0], argv[1]);
	Opt->Ploidy = atoi(argv[1]);

	if(!(Opt->Ploidy == 1 || Opt->Ploidy == 2))
	{
		printf("Ploidy must be 1 or 2.\n");
		exit(1);
	}
}

void SetDrift(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 1);

	if(Opt->Drift == TRUE)
		Opt->Drift = FALSE;
	else
		Opt->Drift = TRUE;
}

void SetEnvVar(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestDouble(argv[0], argv[1]);
	Opt->EnvVar = atof(argv[1]);
}


void SetMutScalar(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestDouble(argv[0], argv[1]);
	Opt->MuteVarScalar = atof(argv[1]);
}

void SetGenerations(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestInt(argv[0], argv[1]);
	Opt->NoGenerations = atoi(argv[1]);
}

void SetChangeFitnessSD(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 3);
	TestInt(argv[0], argv[1]);
	TestDouble(argv[0], argv[2]);

	Opt->ChangeFitnessGen = atoi(argv[1]);
	Opt->ChangeFitnessSD = atof(argv[2]);
}

void SetMoveFitnessSD(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 3);
	TestInt(argv[0], argv[1]);
	TestDouble(argv[0], argv[2]);

	Opt->MoveFitnessSDGen = atoi(argv[1]);
	Opt->MoveFitnessSD = atof(argv[2]);
}


void SetSeed(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	(void)sscanf(argv[1], "%lu", &Opt->Seed);
}

void SetPheotypeVar(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestDouble(argv[0], argv[1]);
	Opt->PheotypeVar  = atof(argv[1]);
}

void SetHeritabitlity(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestDouble(argv[0], argv[1]);
	Opt->Heritabitlity = atof(argv[1]);
}

void SetStandardisePhenotype(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 1);

	if(Opt->StandardisePhenotype == TRUE)
		Opt->StandardisePhenotype = FALSE;
	else
		Opt->StandardisePhenotype = TRUE;
}


void SetOptimumMoveType(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);

	if(strcicmp(argv[1], "fitness") == 0)
	{
		Opt->MoveType = FITNESS;
		return;
	}

	if(strcicmp(argv[1], "phenotype") == 0)
	{
		Opt->MoveType = PHENOTYPE;
		return;
	}

	if(strcicmp(argv[1], "continuous") == 0)
	{
		Opt->MoveType = CONTINUOUS;
		return;
	}
}

void SetAssortativeMating(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestDouble(argv[0], argv[1]);
	
	Opt->AssortativeMating = atof(argv[1]);

	if(Opt->AssortativeMating < 0 )
	{
		printf("AssortativeMating must be >0\n");
		exit(1);
	}
}

void SetInitFitnessSD(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestDouble(argv[0], argv[1]);
	
	Opt->InitFitnessSD = atof(argv[1]);

	if(Opt->InitFitnessSD < 0 )
	{
		printf("InitFitnessSD must be >0\n");
		exit(1);
	}
}

optfun GetOptFun(char* Str)
{
	int Index;

	Index = 0;
	while(CMD_LIST[Index].Cmd != NULL)
	{
		if(strcicmp(Str, CMD_LIST[Index].Cmd) == 0)
			return CMD_LIST[Index].Fun;

		Index++;
	}

	printf("Cmd not valid %s\n", Str);
	exit(1);
	return NULL;
}

int Comment(int argc, char **argv)
{
	if(argc == 0)
		return TRUE;

	if(argv[0][0] == '#')
		return TRUE;

	return FALSE;
}

void SetMutationRatePerGamete(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestDouble(argv[0], argv[1]);
	
	Opt->MutationRatePerGamete = atof(argv[1]);

	if(Opt->MutationRatePerGamete < 0 )
	{
		printf("MutationRatePerGamete must be >0\n");
		exit(1);
	}
}


void SetRecombinationPoissonMean(OPTIONS* Opt, int argc, char** argv)
{
	TestNoArg(argv[0], argc, 2);
	TestDouble(argv[0], argv[1]);
	
	Opt->RecombPoissonMean = atof(argv[1]);

	if(Opt->RecombPoissonMean < 0 )
	{
		printf("RecombPoissonMean must be >0\n");
		exit(1);
	}
}


void PassOptFile(OPTIONS* Opt, char* FName)
{
	TEXTFILE *TF;
	int Index, Tokes;
	char **Passed;
	optfun funpt;

	TF = LoadTextFile(FName, FALSE);

	Passed = (char**)SMalloc(sizeof(char*) * TF->MaxLine);
	
	for(Index=0;Index<TF->NoOfLines;Index++)
	{
		ReplaceChar('\r', '\0', TF->Data[Index]);
		ReplaceChar('\n', '\0', TF->Data[Index]);

		Tokes = MakeArgv(TF->Data[Index], Passed, TF->MaxLine);

		if(Comment(Tokes, Passed) == FALSE)
		{
			funpt = GetOptFun(Passed[0]);
			funpt(Opt, Tokes, Passed);
		}
	}

	FreeTextFile(TF);
	free(Passed);
}

void PrintHelp()
{
	int Index;

	Index=0;
	while(CMD_LIST[Index].Cmd != NULL)
		Index++;

}