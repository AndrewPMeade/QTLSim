import sys
import pandas as pd
import numpy as np
from scipy import stats
import re

MIN_ROW = 0
KEY_ROW = ["P Survived", "Heritability", "w"]

START_MOVE = 4999

StartOfEnd = 5499


def LoadData(FName):

	Header = None
	Data = []
	with open(FName, "r") as FIn:
		for Line in FIn:
			Line = Line.rstrip()

			if Header != None:
				Data.append(Line.split("\t"))

			if re.match("^Generation\t", Line):
				Header = Line.split("\t")

			if re.match("^Pop all died.", Line):
				print(FName, "Population died.", sep='\t')
				sys.exit(1)

	df = pd.DataFrame(Data, columns=Header)
	df = df.apply(pd.to_numeric)
	df = df[df.Generation >= MIN_ROW]
	
	return df

def FindEndPoint(Start, Target, List):

	for i in range(Start, len(List)):
		if List[i] < Target:
			return i

def GetReg(x,y):

	if len(x) == 1:
		return -1, -1, -1

	Reg = stats.linregress(x, y=y)

	return Reg.slope, Reg.intercept, Reg.rvalue

def GetEndMean(Data, Header):
	DList = list(Data[Header])
	return np.mean(DList[StartOfEnd:])

def GetSlopeMean(Data):
	return GetEndMean(Data, "Surviving Slope"), GetEndMean(Data, "Surviving Standardised Slope")

def PrintHeader():
	Header = ["Input File", "GenotypeVar end", "GenotypeVar 95% end", "GenotypeVar start", "GenotypeVar end mean"]
	Header.extend(["GenotypeVar 10 slope"])

	Header.extend(["Phenotype 10 slope"])
	Header.extend(["Phenotype End slope"])

	Header.extend(["Phenotype Ratio 10 slope"])
	Header.extend(["Phenotype Ratio End slope"])

	Header.extend(["W End Mean"])

	Header.extend(["Surviving Slope End Mean", "Surviving Standardised Slope End Mean"])

	Header.extend(["Surviving Population End Mean"])
	
	print(*Header, sep='\t')


def ProcData(FName):

	Data = LoadData(FName)

	PrintHeader()

	GenotypeVar = list(Data["Genotype Var"])
	
	
	GenotypeVarEnd = GenotypeVar[StartOfEnd:]
	
	Mean = np.mean(GenotypeVarEnd)
	EndPoint = FindEndPoint(START_MOVE, Mean, GenotypeVar)

	Range = GenotypeVar[START_MOVE] - Mean
	EndPoint95 = FindEndPoint(START_MOVE, GenotypeVar[START_MOVE]-(Range * 0.95), GenotypeVar)

	
	print(FName, sep='\t', end='\t')

	print(EndPoint, EndPoint95, GenotypeVar[START_MOVE], Mean, sep='\t', end='\t')

	Slope, _, _ = GetReg(Data["Generation"][START_MOVE:START_MOVE+10], GenotypeVar[START_MOVE:START_MOVE+10])
	print(Slope,  sep='\t', end='\t')

	
	Phenotype = list(Data["Phenotype Mean"])
	PhenotypeDecay = Phenotype[START_MOVE:EndPoint+1]
	PhenotypeEnd = Phenotype[StartOfEnd:]
	PhenotypeDecay10 = GetReg(Data["Generation"][START_MOVE:START_MOVE+10], Phenotype[START_MOVE:START_MOVE+10])
	PhenotypeRunning = GetReg(Data["Generation"][StartOfEnd:], PhenotypeEnd)
	
	print(PhenotypeDecay10[0], PhenotypeRunning[0], sep='\t', end='\t')

	PhenotypeRatio = list(Data["Phenotype Ratio"])
	RatioDecay = PhenotypeRatio[START_MOVE:EndPoint+1]
	RatioEnd = PhenotypeRatio[StartOfEnd:]
	RatioDecay10 = GetReg(Data["Generation"][START_MOVE:START_MOVE+10], PhenotypeRatio[START_MOVE:START_MOVE+10])
	RatioRunning = GetReg(Data["Generation"][StartOfEnd:], RatioEnd)
	print(RatioDecay10[0], RatioRunning[0], sep='\t', end='\t')

	WEnd = list(Data["w"])[StartOfEnd:]
	print(np.mean(WEnd), end='\t')

	print(*GetSlopeMean(Data), sep='\t', end='\t')
	
	SurvivedEnd = list(Data["P Survived"])[StartOfEnd:]

	print(np.mean(SurvivedEnd), sep='\t', end='\t')

	print()

ProcData(sys.argv[1])