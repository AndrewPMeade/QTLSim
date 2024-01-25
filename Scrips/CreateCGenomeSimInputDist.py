import sys
import itertools
import numpy

NO_SIM = 15000

def OuputModel(No, Ploidy, MoveFitness, ChangeFitness, Heritability, PhenotypeVar, MutationScalar, MoveFitSD, Seed):

	with open(f"Input-{No:05d}.txt", "w") as FOut:
		print("Sample 1", file=FOut)
		print("PopSize 1000", file=FOut)
		print("NoQTL 1", file=FOut)
		print("Generations 6500", file=FOut)

		print(f"Ploidy {Ploidy}", file=FOut)

		if MoveFitness == True:
			print(f"MoveFitnessSD 5000 {MoveFitSD}", file=FOut)
		
		print(f"ChangeFitnessSD 5000 {ChangeFitness}", file=FOut)	
		
		print(f"PhenotypeVar {PhenotypeVar}", file=FOut)	

		print(f"MutationScalar {MutationScalar}", file=FOut)

		print(f"Heritability {Heritability}", file=FOut)	


		print(f"Seed {Seed}", file=FOut)
		
def TestWeibeull(Size, Shape, Scale):

	y = numpy.random.weibull(Shape, size=Size)
	y = Scale * y
	print(*y, sep='\n')

	sys.exit(0)


def CutList(X, Size, min_val, max_val):
	f_list = []

	for x in X:
		if x>=min_val and x<=max_val:
			f_list.append(x)
	
	if len(f_list) < Size:
		print("could not gnerate enough points")
		sys.exit(1)

	return f_list[:Size]

def GenExp(Alpha, Size, min_val, max_val):
	y = numpy.random.exponential(Alpha, size=Size*5)
		
	return CutList(y, Size, min_val, max_val)

def GenWeibeull(Shape, Scale, Size, min_val, max_val):

	y = numpy.random.weibull(Shape, size=Size*2)
	y = Scale * y
	
	return CutList(y, Size, min_val, max_val)

def GenGamma(Shape, Scale, Size, min_val, max_val):
	y = numpy.random.gamma(Shape, scale=Scale, size=Size*2)
	return CutList(y, Size, min_val, max_val)

Ploidy = 2
MoveFitness = True

Seeds = numpy.random.randint(0, numpy.iinfo(numpy.int32).max, size=NO_SIM)

#mu = -6.13; sigma = 1.21
#MutationScalars = numpy.random.lognormal(mu, sigma, NO_SIM)
#MutationScalars = numpy.random.gamma(1.1, scale=0.005, size=NO_SIM)
MutationScalars = GenGamma(1.1, 0.005, NO_SIM, 0.0001, 100000)


ChangeFitness = numpy.random.uniform(0.06, 0.8, NO_SIM)
MoveFitnessSD = numpy.random.uniform(1, 2, NO_SIM)
Heritability = GenWeibeull(1.6, 0.4, NO_SIM, 0.01, 0.99)

#PhenotypeVar = numpy.random.uniform(0.0001, 0.09, NO_SIM)
PhenotypeVar = GenExp(0.018, NO_SIM, 0.01, 0.09)



Header = ["SimNo", "Ploidy", "MoveFitness", "ChangeFitnessSD", "Heritability", "PhenotypeVar", "MutationScalar", "MoveFitnessSD", "Seed"]
print(*Header, sep='\t', end='\n')

for No, (PhenotypeVar_P, H, CFit, MutS, MoveFitSD, Seed) in enumerate(zip(PhenotypeVar, Heritability, ChangeFitness, MutationScalars, MoveFitnessSD, Seeds)):
	Model = [No, Ploidy, MoveFitness, CFit, H, PhenotypeVar_P, MutS, MoveFitSD, Seed]
	OuputModel(*Model)
	print(*Model, sep='\t')
