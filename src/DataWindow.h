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

#ifndef DATA_WINDOW_H
#define DATA_WINDOW_H

typedef struct
{
	double *Vect;
	int Size;
	int No;
} DATA_WINDOW;

DATA_WINDOW*	CreateDataWindow(int Size);
void			FreeDataWindow(DATA_WINDOW* DataWindow);
void			DataWindowAdd(DATA_WINDOW* DataWindow, double Val);
double			DataWindowAve(DATA_WINDOW* DataWindow);


#endif
