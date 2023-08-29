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
#include "DataWindow.h"

DATA_WINDOW* CreateDataWindow(int Size)
{
	DATA_WINDOW* DW;

	DW = (DATA_WINDOW*)SMalloc(sizeof(DATA_WINDOW));
	
	DW->Size = Size;
	DW->Vect = (double*)SMalloc(sizeof(double) * DW->Size);
	DW->No = 0;

	return DW;
}

void FreeDataWindow(DATA_WINDOW* DataWindow)
{
	free(DataWindow->Vect);
	free(DataWindow);
}

void DataWindowAdd(DATA_WINDOW* DataWindow, double Val)
{
	int Pos;

	Pos = DataWindow->No % DataWindow->Size;

	DataWindow->Vect[Pos] = Val;

	DataWindow->No++;
}


double DataWindowAve(DATA_WINDOW* DataWindow)
{
	double Ave;
	int Index, Max;

	
	Max = DataWindow->Size;
	if(DataWindow->No < DataWindow->Size)
		Max = DataWindow->No;

	Ave = 0;
	for(Index=0;Index<Max;Index++)
		Ave += DataWindow->Vect[Index];

	return Ave / Max;
}