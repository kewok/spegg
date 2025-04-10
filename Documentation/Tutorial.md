# sPEGG Tutorial

Version 0.5 

Were you able to get everything from the Quickstart to work? Great :thumbsup: Once you have everything in place, this tutorial will walk you through what was going on behind the curtains there. Another goal of this tutorial is to provide you with a template that you can start modifying to simulate your model.

## Object-Oriented Programming (OOP) Kick-off:

Initially, some familiarity with object-oriented programming would help ([this](http://sepwww.stanford.edu/sep/jon/family/jos/oop/oop1.htm) is an excellent and very accessible introduction; it's written to introduce students to java. What you would learn from this is how the author made useful real-life analogies to multiple OOP concepts, which are supposed to be hard to grasp but become easier with these analogies. But don't worry - what it says is almost entirely applicable to C++ as well).

Ironically, the official author group of Java, **_[Oracle](https://www.oracle.com/])_** also created a **[tutorial web page for OOP in Java](https://docs.oracle.com/javase/tutorial/index.html)**. Please refer to it if you desire to delve even deeper into the OOP world with the Java language! The topics that may be most relevant to this project are **_[Learning the Java Language](https://docs.oracle.com/javase/tutorial/java/index.html)_** and **_[Collections](https://docs.oracle.com/javase/tutorial/collections/index.html)_**.

However, for `C++`, there is no official tutorial page from its author, but [this tutorial page](https://www.tutorialspoint.com/cplusplus/cpp_object_oriented.htm) also achieves the likely purpose to the Oracle's Documentation.

## My first **sPEGG** project

**sPEGG** makes extensive use of the **[thrust](https://github.com/thrust/thrust)** library to facilitate parallelization on the GPU. As such, we encourage users who are eager to try out **sPEGG** to become at least somewhat acquainted with how [thrust works](https://nvidia.github.io/cccl/thrust/). 

The rest of this guide assumes you are looking to create a GPU simulation using forward-time population genetics for your eco-evolutionary model. The basic design is more or less identical for the CPU version; just be sure to checkout [the CPU version of the code](https://github.com/kewok/spegg/tree/CPU_version) at the start, use ```host_vector``` throughout instead of ```device_vector```, and modify the random number generator as needed along the way. See the examples in [the CPU version](https://github.com/kewok/spegg/tree/CPU_version/Examples) to get a feel for how to do that.

Some of the tasks described here in the tutorial can be streamlined using a program like [NSight Eclipse](https://docs.nvidia.com/cuda/nsightee-plugins-install-guide/index.html). However, in what follows, I will assume all you are using is a basic text editor, a file browser, and the [prerequiste tools](/README.md#prereqs).

To get started, first obtain a copy of the source-code for sPEGG's GPU version into your working directory. If you have [git](http://git-scm.com/), go to your terminal and clone the git repository as we did in the Quickstart:

```sh
$ git clone --recursive https://github.com/kewok/spegg
```

Alternatively, you could just download the [zip folder](https://github.com/kewok/spegg/archive/master.zip) and unpack it to your working directory.

If you want, rename the folder as needs be (e.g., "My_sPEGG_project") and navigate into the folder.

The whole tutorial can be skimmed to give you a broad overview of the logic of **sPEGG** projects. If you are new to C++, it might take a few hours or separate sessions to get through it with your own example up and running. If you haven't yet been confident with C++, visit the section [Object-Oriented Programming (OOP) Kick-off](#object-oriented-programming-(oop)-kick-off).

## Code organization in a sPEGG project

The hierarchical organization of the directories should be familiar to those who have worked on other large coding projects. There is a header folder that stores the .h files and a src folder that stores the CUDA/C (.cu) files. Finally, code related to generating the executable and the various input files will be stored in your main directory. 

In the main codebase, within each of these folders, the code is organized thematically. Code relating to the ecological variables (e.g., abiotic factors such as nutrient availablility, climactic input etc... as well as biological variables like prey or predator densities of species that are not tracked on an individual-level) are stored in the folder ```environ```. Code for conducting basic numerical operations (calculating summary statistics, random number generation, etc...) is stored in the folder ```math```. Miscellaneous code that mostly extends thrust functionality and does some basic vector manipulation (e.g., removing duplicates, flattening 2D arrays, etc...) as well as basic diagnostics such as a timer and code for checking GPU memory availability are stored in the folder ```util```. Finally, code relating to biological processes is stored in the folder ```species.``` Within species, code pertaining to reproduction events are stored in ```add_kids```, that pertaining to selection and updating phenotypes in ```update``` and (if migration between demes is modeled) that pertaining to  migration in a folder ```movement```. Generic code that applies to all species are stored as loose files within the ```species``` folder.

In general, we've found that when setting up a **sPEGG** project, creating a folder to store your specie's code separately helps. The remaining code base can be completely agnostic about your additions. This enforces modularity in your code and allows you to work with most of **sPEGG**'s functionality much as you might a C++ library. 

### Creating an object for my species

Now you're ready to begin building your simulation! For this tutorial, we'll just simulate a stochastic change in allele frequencies at two neutral loci in two demes (i.e., simulate drift). We'll assume the demes are unlinked by migration. Let's pretend these are penguin populations living on two distant islands unlinked by migration, and that the neutral loci determine the phenotypic value for crown color, which we assume governs neither survival nor reproduction.

The first thing you'll want to know is that **sPEGG** organizes individuals in a species using the [inds](/header/species/inds.h) object. [inds](/header/species/inds.h) is very basic - it merely:

>    1. stores data on individuals, 
	
>    2. let's you export this data to a csv file, 
	
>    3. sorts individuals by deme, and 
	
>    4. keeps track of who's alive and who's dead. 

Because **sPEGG** aims to also accommodate deterministic models, we actually will use a derived variant of ```inds``` called ```inds_stochastic``` instead. The class ```inds_stochastic``` is identical to its parent ```inds```, except that it has associated with it a random number generator for the GPU, placeholders for parameters determining the baseline survival rates and fecundities, as well as unspecified reproduction and updating methods associated with it.

### The header file

Obviously, by itself the functionality of ```inds_stochastic``` won't get you very far in your simulation. So first, you will want to create a header file (.h) that lets you add a few more things for the species you want to simulate (in our case, penguins). Let's store this header file in a folder we'll create called ```header``` inside your ```My_sPEGG_project``` folder. Once the folder is created, navigate into ```header``` and create a file called Penguins.h there (i.e., so you have the following path: ```My_sPEGG_project/header/Penguins.h```). The contents of this header file will look something like:

```c++
#ifndef PENGUINS_H
#define PENGUINS_H

#include <species/inds_stochastic.h>

#include "Penguin_Parents.h"
#include <species/add_kids/neonates_class.h>

class Penguins : public inds_stochastic
	{
	public:
		Penguins(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val);
		void addKids();

		using inds_stochastic::update;
		void update(inds_stochastic **species);

	protected:
		void initialize_demes();
		void setPhenotype(int index, int n);
		void assignSex(int index, int n);
	};
#endif
```
Let's walk through some salient features that form the pre-requisites for our class. First, because the ```Penguins``` class is based on ```inds_stochastic```, we need to import the header file that defines the ```inds_stochastic``` class (```<species/inds_stochastic.h>```, which is stored in the main **sPEGG** directory). Second, when we simulate reproduction we need to be able to keep track of who the parents and offspring are. To do this, we will use a parents and neonates class, and so we import their header files. For this example, the parents class would be customized for penguins as we will describe below. By contrast, in this example, for the neonates class we will just use the [default base class for this](/header/species/add_kids/neonates_class.h) provided by **sPEGG**. The class definition for neonates can be found inside the ```species``` subfolder within the header folder in your main **sPEGG** directory.

We need to declare the methods used by our ```Penguins``` class. The first method is the constructor that creates our class. The second method ```addKids()``` is where we will specify the routines used in reproduction. The third method ```update(inds_stochastic **species)``` will handle how our phenotypes change between reproduction events. We also have to unhide the ```update``` function from ```ind_stochastic```. This is because the ```update``` method is actually defined in a couple of different ways inside the ```inds_stochastic``` class to allow some flexibility for the user. We'll get into this in greater detail below.

Let's look at the protected components. These are data and members that can only be called from within the ```Penguins``` class. The methods deal with initializing our simulation, and we'll also want to incorporate the newborns into our data structures by assigning them phenotypes using setPhenotype(). Since penguins are a sexually reproducing species, it's important to keep track of who the males and females in the population are. We'll therefore also create a method for assigning sex to newborns.

### The .cu file

Now let's actually write the code to handle the business logic for Penguins. Go back to the **My_sPEGG_project** directory and create and navigate to a folder called ```src/```. There, we'll make a file called Penguins.cu. In what follows, we will piece together different parts described above in this .cu file step-by-step.

The first thing you want to do is include some header files.
<a name="header">

```c++
#include "Penguins.h"

#include <species/add_kids/genotype_phenotype_map.h>
#include <species/update/updatebehavior.h>

#include <thrust/copy.h>
#include <thrust/fill.h>
```
</a>

The first header is for the class definition. The second header specifies an object that will handle translating the genotypic values of individuals into their phenotypic values. The third header is for an external class that models phenotypic change.  Each of the remaining headers specify a few basic operations (like assigning a fixed value to all elements of your data vectors) that are to be handled on the GPU via thrust.

#### The constructor

The constructor is where you will normally handle initialization. Our constructor will perform the following operations.

1. Allocate memory on the GPU and assign all individuals a sequential ID, set everyone's status as alive, and read in the parameter values (using the functionality from ```inds```).

2. Seed and initialize your random number generator (using the functionality from ```inds_stochastic```).

3. Initialize the age of all individuals to zero.

4. Assign half the individuals to deme #0, and the other half to deme #1. This will be handled by a [helper method](#helper_methods) initialize_demes() that will determine the number of individuals in each deme.

5. Assign allelic values at the maternally inherited loci (represented by the array of vectors fgenotype, where the *i*th element of the array stores the allelic values of individuals at the *i*th locus inherited from their female parent) and the paternally inherited loci (represented by the array of vectors mgenotype storing allelic values inherited from the male parent) according to a gaussian distribution (thus we are assuming a continuum-of-alleles approximation) with mean=0 and standard deviation=0.1.

6. Based on the allelic values, determine everyone's phenotype. We'll use the ```setPhenotype()``` function [described below](#gen_phen_map), which is used to determine offspring phenotypes at birth, for initialization.

7. Assign the sex of individuals. For the initialization step, we'll simply assume odd numbered individuals to be male, and even numbered individuals to be female.

The good news is that we don't have to write code to do all of this. Because **sPEGG** follows an object-oriented design, step (1) is handled by the base class ```inds``` and step (2) is handled by Penguin's base class (not to be confused with THE ```parents_class```) ```inds_stochastic```. Notice that these are steps that have to be carried out by most species we'll model, and how you go about doing that is basically the same. By contrast, the remaining steps might warrant customization (e.g., you might want the age distribution to resemble the distribution observed in the field at a given point in time), so they still have to be specified from within the ```Penguins``` class.

Here is the code:

```c++
Penguins::Penguins(int size_val, int maxsize_val, int seed_val, int ndemes, int species_ID_val) : inds_stochastic(size_val, maxsize_val, seed_val, ndemes, species_ID_val)
	{
	// Assume everyone starts at age = 0.
	thrust::fill(age.begin(), age.begin() + size, 0);

	initialize_demes();
	
	// Specify the genetics by assuming allelic values are gaussian-distributed
	for (int i=0; i < nloci; i++)
		{		
		draw_gaussian(size, 0, 0.1, fgenotype[i], gen);
		draw_gaussian(size, 0, 0.1, mgenotype[i], gen);
		}

	//Set phenotype
	setPhenotype(0, size);

	// To start, assign odd numbered individuals to be male, even numbered individuals to be female
	thrust::device_vector<int> twos(size);
	thrust::fill(twos.begin(), twos.end(), 2);
	thrust::transform(id.begin(), id.begin() + size, twos.begin(), sex.begin(), thrust::modulus<int>());
	}
```
Note in the last line of code that we pass ```id.begin() + size``` instead of ```id.end()``` to thrust's ```transform()``` function. The reason for this is that the vector ```id``` is of length ```maxsize```. However, the actual number of individuals on whose calculations we are operating at any given point in time is ```Penguins->size```, not ```maxsize```. Therefore, if we were to pass ```id.end()``` instead here, this would be longer than the vector ```twos``` and cause a memory error. Just be sure to keep this distinction in mind as you put together your simulation code.

You're now ready to start performing calculations on the individuals of our class ```Penguins```.

### Mating

In this tutorial, mating is simulated using the method ```addKids()```. Simulating mating consists of a few basic steps:

1. [Identify parents and their reproductive contributions](#id_potential_parents)
2. [Simulate the transmission of genetic information from parents to offspring](#simulate_gene_transfer)
3. [Incorporate the newborns into our species data object](#integrate_newborns)

Below, we describe each of these steps.
<a name="id_potential_parents">
</a>

####  Determine who the potential parents are
We will create a class to manage potential parents and perform basic operations on them (such as determining their reproductive potential). As hinted at [above](#header) we will use a derived version of the [Parents class](/header/species/add_kids/parents_class.h) to accomplish this.

The header file for our parents class will look something like this:

```c++
#ifndef PENGUIN_PARENTS_CLASS_H
#define PENGUIN_PARENTS_CLASS_H

#include <species/add_kids/parents_class.h>

class Penguin_Parents : public Parents
	{
	class Penguins *species;
	public:
		Penguin_Parents(Penguins *species);
		
	protected:
		int FECUNDITY_PHENOTYPE_INDEX;
		void determine_probability_individual_becomes_female_parent();
		void determine_probability_individual_becomes_male_parent();
	};

#endif
```
The ```Penguin_Parents``` class is a special case of a base ```Parents``` class. The base ```Parents``` class handles most of the generic functionality of determining who can be parents, and how likely it is that an individual parent will sire a particular offspring. We do, however, have to specify two penguin-specific functions: ```determine_probability_individual_becomes_female_parent()``` and ```determine_probability_individual_becomes_male_parent()```. These functions will assume that whatever the fecundity score of a given individual is will reflect it's probability of reproducing in the next generation.

We can now implement the functionalities of the parents class. Inside your ```src/``` directory in My_sPEGG_project, create a file called ```Penguin_Parents.cu```. The constructor to the ```Penguins_Parents``` class allows us to link the classes ```Penguin_Parents``` and ```Penguins``` and looks something like:
```c++
Penguin_Parents::Penguin_Parents(Penguins *species) : Parents(species)
	{
	FECUNDITY_PHENOTYPE_INDEX = demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"];
	}
```
Notice here we have to use some terminology specific to ```Penguins```. In particular, the variable ```FECUNDITY_PHENOTYPE_INDEX``` here is a species-specific variable that is specified uniquely for penguins. In **sPEGG**, you can access such species-specific variables with reference to a [demeParemeters](#demeParameters) object (more on this below).

For the remaining methods inside ```Penguin_Parents```, we will make use of the [```transform``` operation in thrust](#https://thrust.github.io/doc/group__transformations.html#gacca2dd17ae9de2f7bcbb8da6d6f6fce4). For females, the determination of their probability of becoming a parent will be determined by multiplying a vector ```will_reproduceF``` (calculated by the base ```Parents``` class) that lists whether a female is reproductively capable, by her anticipated fecundity. Here is the relevant function definition:

```c++
void Penguin_Parents::determine_probability_individual_becomes_female_parent()
	{
	thrust::multiplies<float> op;
	thrust::transform(will_reproduceF.begin(), will_reproduceF.begin() + size, phenotype[FECUNDITY_PHENOTYPE_INDEX].begin(),  probability_individual_becomes_female_parent.begin(), op);
	}
```
We can similarly define an exactly analogous function for ```determine_probability_individual_becomes_male_parent()``` as:


```c++
void Penguin_Parents::determine_probability_individual_becomes_male_parent()
	{
	thrust::multiplies<float> op;
	thrust::transform(will_reproduceM.begin(), will_reproduceM.begin() + size, phenotype[FECUNDITY_PHENOTYPE_INDEX].begin(),  probability_individual_becomes_male_parent.begin(), op);
	}
```

Once we have defined our Penguin_Parents class, simply declare an instance of it inside ```Penguins```, let's call it ```exampleParents```, and initialize it:

```c++
Penguin_Parents *exampleParents;
exampleParents = new Penguin_Parents(this);
```
The argument ```this``` refers to an instance of our Penguins class. 

Our first step will be to determine what the reproductive output should be for each individual. By default, this function identifies females, looks for a parameter in the parameter inputs specifying the per-capita fecundity for females, and sums the maximum reproductive output possible for all individuals. On occasion, we may have to modify this function to distinguish reproductively viable individuals from individuals that are not. For instance, in real-world penguin populations, females must often have a certain amount of energy reserves to be able to reproduce. Thus, behind the scenes, the method ```setup_parents()``` will have to on the phenotypes (e.g., current energy reserves) of individuals to determine this value. Although we won't pursue this approach here, for species with complicated rules governing the determination of parents, these rules can be specified further in the derived, species-specific ```Parents``` class.

In most situations, we cannot allow all individuals to reproduce as many offspring as they are able to, because this will lead to exponential population growth which is neither biologically realistic nor computationally tractable (even on a GPU!). Therefore, to check the per-capita fecundity, the function ```setup_parents()``` will rescale parental reproductive output to no longer exceed the maximum population size in each deme. The sum of these scaled contributions will be stored in a variable called ```Potential_Number_of_Kids``` inside the Parents class. If there are no spaces left for new individuals, the value of ```Potential_Number_of_Kids``` is set to 0. These operations will all be called from within the ```addKids()``` method. 
[//]: # (Generally, you will want to set the maximum population size well above the equilibrium number of individuals in your simulation so as to not artificially introduce this sort of density-dependent regulation into your model.)

<a name="simulate_gene_transfer">

#### Simulate the transmission of genetic information from parents to offspring
</a>
Provided there will be a reproduction event in this time step, we now want to create a class that manages the newborns for us. To do this, we will use the EggsNeonates class imported in the [header](#header) much as we used the Parents class to help us in the previous step. We'll initialize this as follows:

```c++
EggsNeonates *exampleNeonates;
exampleNeonates = new EggsNeonates (this, exampleParents->kids_per_mom);
```
As before, the argument ```this``` refers to our Penguins class. The initialization involves a new argument ```kids_per_mom``` from our ```Parents``` class. This is an integer vector which, for each reproductive female in our species, stores the number of offspring she will birth in this time step. It is calculated during the setup_parents() function described above.

Once the class for newborns is initialized, the transmission of genetic information is handled by its method ```inherit_genotypes```:

```c++
exampleNeonates->inherit_genotypes(exampleParents->probability_individual_becomes_female_parent,  exampleParents->probability_individual_becomes_male_parent);
```
These two variables, ```probability_individual_becomes_female_parent``` and ```probability_individual_becomes_male_parent``` are members of the ```Parents``` class which weigh the reproductive contribution to the next generation of females and males, respectively. The function ```inherit_genotypes()``` simulates the mating/matching of parents, mutation and recombination for each offspring. The default behavior is to model sexual reproduction of a diploid species; versions of exampleNeonates that handle asexual reproduction, where we would only be interested in the reproductive contribution of female parents (since everyone would be female) are also supported, and future expansions will include support for other genetic systems such as haplodiploidy and horizontal gene transfer.

<a name="integrate_newborns">

#### Integrate the newborns into our species
</a>

Once the offspring have inherited the genotypes, we will determine their phenotypes at birth using the `setPhenotype(exampleNeonates -> previous_pop_size, exampleNeonates->Total_Number_of_Neonates)` function from our initialization routine. The differences between this call to `setPhenotypes()` and the previous call during initialization are twofold. Instead of the index individual being zero, the function will begin starting with the number of individuals in the population prior to the reproduction event (which will be stored in `exampleNeonates -> previous_pop_size` during the construction of exampleNeonates). Second, it will evaluate the phenotypes for the total number of new born offspring, rather than the total number of individuals in the simulation. This is performed using the [setPhenotype function](#gen_phen_map) whose discussion we defer for now. We will also assign the sex of the newborns using the [helper function](#helper_methods) `assignSex()`.

Finally, we perform a few book-keeping operations. To make sure all individuals in a given deme are contiguous in GPU memory, we will call a function `sortByDeme()` from inds, and recalculate the population sizes in each deme using `demeCalculations()`, also from the inds base class. The parents and neonates classes are then cleared from memory. 

Putting all the components described above together, the function addKids() which simulates mating and  reproduction then looks like:

```c++
void Penguins::addKids()
	{   
	demeCalculations();

	Penguin_Parents *exampleParents;
	exampleParents = new Penguin_Parents(this);
	exampleParents->setup_parents();

	if (exampleParents->Potential_Number_of_Kids > 0)
		{
		EggsNeonates *exampleNeonates;
		exampleNeonates = new EggsNeonates (this, exampleParents->kids_per_mom);
		exampleNeonates->inherit_genotypes(exampleParents->probability_individual_becomes_female_parent,  exampleParents->probability_individual_becomes_male_parent);

		setPhenotype(exampleNeonates->previous_pop_size, exampleNeonates->Total_Number_of_Neonates);
		assignSex(exampleNeonates->previous_pop_size, exampleNeonates->Total_Number_of_Neonates);
		sortByDeme();
		delete exampleNeonates;     
		}
	demeCalculations();
	delete exampleParents;
	}
```

<a name="gen_phen_map">

### The genotype-phenotype map
</a>

One function that is needed both during initialization and in simulating reproduction is a representation of the genotype-phenotype map. This segment of the code can actually be somewhat involved, so we have deferred discussing it until here. The basic idea is that for each phenotype, we create an object that will use a genotype-phenotype map to translate the genotypic values of individuals into their phenotypic values. 

Why would we do things this way? When multiple traits are involved (as they often are in models of eco-evolutionary dynamics), there are a lot of different potential genotype-phenotype maps. As a result, it might be helpful to create separate classes to manage these maps, rather than spell out the logic of each genotype-phenotype map for all traits inside ```Penguins```. This could facilitate modularity in our code, and also allow us to potentially reuse a genotype-phenotype map for another species or change some of the rules for updating phenotypes without having to change the code inside the Penguins class. This tries to apply the so-called [strategy pattern](https://web.archive.org/web/20160625125844/http://r3dux.org/2011/07/an-example-implementation-of-the-strategy-design-pattern-in-c/) to our simulation. This isn't always necessary, but for illustrative purposes we show how this kind of approach could be implemented in case you're interested in trying it. We'll discuss a little about the scope associated with pre-written **sPEGG** code.

We represent a genotype-phenotype map using the class GenotypePhenotypeMap. The class definition can be found in [/header/species/add_kids/genotype_phenotype_map.h](/header/species/add_kids/genotype_phenotype_map.h) The create_genotype_phenotype_map() function will initialize the particular genotype phenotype-map handler we want to create for phenotype *i*. We will then perform the actual calculations required to set the individual's phenotypes based on their genetic values. This class is defined as follows:

```c++
class GenotypePhenotypeMap
	{
	public:
		static GenotypePhenotypeMap *create_genotype_phenotype_map(inds *species, int phenotype_index, int index_case, int num_kids);

		 GenotypePhenotypeMap(inds *species, int phenotype_index, int index_case, int num_kids)
						{
						this->phenotype_index = phenotype_index;
						this->Parameters = species->demeParameters->GeneticArchitecture->phen_gen_map_parm[phenotype_index];
						this->index_case = index_case;
						this->num_kids = num_kids;
						}
	
		virtual void calculate_phenotype(inds *species)=0;

	protected:
		GenotypePhenotypeMapParameters *Parameters;

		int phenotype_index;
		int index_case;
		int num_kids;

		curandGenerator_t gen;		
	};
```
The variables involved in this class can be described as follows. The class GenotypePhenotypeMapParameters is responsible for storing the parameters that go into determining the genetic architecture. For instance, suppose you define a genotype-phenotype map for an additive trait as ```F(x) = A + B * x_1 + C * x_2```, where the vector *x* consist of the allelic values at each locus and F(x) describes the genetic (breeding) value of your phenotype. Associated with this map are coefficients A, B and C. Then the class GenotypePhenotypeMapParameters will simply be a class that stores the values for [A,B,C]. It will also specify what loci go into the genotype-phenotype map. This is discussed a little further in the section in this tutorial on [input files](#demesettings). 

The other variables are fairly straightforward. ```phenotype_index``` is the index of the current phenotype whose genotype-phenotype map we are modeling. ```index_case``` identifies the index of the first individual that needs to have its phenotype calculated. ```num_kids``` specifies the number of individuals whose phenotypes at birth (or upon initialization) we have to determine. ```gen``` is simply a random number generator that runs on the GPU.

The ```GenotypePhenotypeMap``` class is responsible for instantiating a specific genotype-phenotype map and evaluating that map using a function called ```calculate_phenotype(inds)```. To call a genotype-phenotype map class from Penguins, we simply define a method:
```c++
void Penguins::setPhenotype(int index, int num_inds_to_calculate)
	{	
	for (int i=0; i < nphen; i++)
		{	
		GenotypePhenotypeMap *genotype_phenotype_map;
		genotype_phenotype_map = genotype_phenotype_map->create_genotype_phenotype_map(this, i, index, num_inds_to_calculate);
		genotype_phenotype_map->calculate_phenotype(this);
		delete genotype_phenotype_map;
		}
	}
```

The for loop will go through each of ```nphen``` phenotypes, doing the calculations for determining what the phenotypic values of individuals should be at that phenotype. The class ```genotype_phenotype_map``` will take care of delegating calculations related to the genotype-phenotype map.

This is all simple enough, but we need to go into some detail about how this is actually implemented in order to use this feature effectively.

<a name="gen_phen_map_header"></a>
The method ```create_genotype_phenotype_map``` can be created following a template in ```src/cuda/species/add_kids/genotype_phenotype_map_TEMPLATE.cu```. For a very simple simulation where we only model, say, the phenotypes directly involved in fitness (fecundity and mortality) and the neutral phenotype (here, crown color), we can define this method inside ```src``` of ```My_sPEGG_project``` as:

```c++
GenotypePhenotypeMap *GenotypePhenotypeMap::create_genotype_phenotype_map(inds *species, int phenotype_index, int index_case, int num_kids)
	{
	if (phenotype_index == species->demeParameters->species_specific_values["FECUNDITY_PHENOTYPE_INDEX"]) 
	   {       
	   return new fecundity_phenotype(species, phenotype_index, index_case, num_kids);
	   }

	if (phenotype_index == species->demeParameters->species_specific_values["MORTALITY_PHENOTYPE_INDEX"])
	   {        
	   return new mortality_phenotype(species, phenotype_index, index_case, num_kids);
	   }

   if (phenotype_index == species->demeParameters->species_specific_values["CROWN_COLOR_INDEX"])
	   {
		   return new crown_color_phenotype(species, phenotype_index, index_case, num_kids);
	   }
	}
```
What this does is create instances of genotype-phenotype map objects, which we define below, that are in charge of managing the translation of genetic data into phenotypic values in our model. 

However, notice again that in order to call ```create_genotype_phenotype_map```, we have to use some terminology specific to ```Penguins```. For instance, the variable ```FECUNDITY_PHENOTYPE_INDEX``` is not defined for all ```inds```, but only for ```inds_stochastic```; more restrictive still is the variable ```CROWN_COLOR_INDEX``` that only makes sense for penguins. Similarly, the generic ```GenotypePhenotypeMap``` class will not know what the classes ```fecundity_phenotype``` and ```mortality_phenotype``` are. This is because the genetic architecture underlying these traits could vary by species. In order for ```GenotypePhenotypeMap``` to be able to make sense of all this, you need to pass the header file where you actually declare, e.g., what exactly ```mortality_phenotype``` and ```fecundity_phenotype``` are. Because each species you model may very well have a different set of idiosyncratic phenotypes (wing length, prey preference, flowering time, etc...), as well as different ways for determining their baseline mortality and fecundity schedules, this header file will likely be specific to your species. Therefore, the top of the genotype_phenotype_map.cu file should have a set of instructions to include the header file such as:

```c++
#include "Penguins_genotype_phenotype_maps.h"
```

Below, we describe what the contents of this header file can look like for our penguins example.

#### Example genotype-phenotype maps.

In the ```create_genotype_phenotype_map()``` function above, we constructed three classes: ```fecundity_phenotype```, ```mortality_phenotype``` and ```crown_color_phenotype``` which actually handle our genotype-phenotype map calculations. To define these, we'll add a header file called ```Penguins_genotype_phenotype_maps.h``` to the ```header``` folder in the ```My_sPEGG_project``` directory.

We'll start by importing headers into Penguins_genotype_phenotype_maps.h. Because each of these genotype-phenotype maps are examples of the more generic GenotypePhenotypeMap class (i.e., they are derived classes), we'll need to include the header file for the GenotypePhenotypeMap parent class definition. 

Now we're ready to define our genotype-phenotype maps. Lets start by defining the fecundity_phenotype class:

```c++
#include <species/add_kids/genotype_phenotype_map.h>

class fecundity_phenotype : public GenotypePhenotypeMap
	{
	public:
		fecundity_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{};
	
		void calculate_phenotype(inds *species);
	};
```

The class inherits from GenotypePhenotypeMap, and implements the ```calculate_phenotype()``` method according to the genotype-phenotype map for the fecundity phenotype. The arguments from the constructor to GenotypePhenotypeMap will be used to querry the parameters involved in determining the genotype-phenotype map and set up other variables of the GenotypePhenotypeMap, including the index that represents the first offspring in the arrays of the Inds class and the total number of offspring. 

<a name="functor_genphenmap">

Now we are finally able to specify our genotype-phenotype map. To do this, we will use what is called a functor. Functors are explained further in the [thrust introduction](https://nvidia.github.io/cccl/thrust/api_docs/algorithms.html). For now, you can think of them as customizable functions `F()` that map vectors like ```x=[x1,x2,x3,...,xn]``` and  ```y=[y1,y2,y3,...,yn]``` to a third vector ```F(x,y) = [F(x1,y1), F(x2,y2), F(x3,y3),...,F(xn,yn)]```. In thrust, functors get specified as structures. 

</a>

Here's what our functor mapping the individual's genotypes to their phenotypes might look like:
```c++
struct fecundity_calculator
	{
	float *map_constant; 
	float *map_coefficient_0;
	float *map_coefficient_1;

	fecundity_calculator(float* map_cons, float* map_coef0, float* map_coef1) : map_constant(map_cons), map_coefficient_0(map_coef0), map_coefficient_1(map_coef1)
	{};

	/* 
	Elements in the tuple.
	---------------------
	0: individual's deme
	1: individual's maternally inherited allelic value at their locus 0
	2: individual's paternally inherited allelic value at their locus 0
	3: individual's maternally inherited allelic value at their locus 1
	4: individual's paternally inherited allelic value at their locus 1
	5: individual's fecundity phenotypic value
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int ind_deme = thrust::get<0>(t);
		thrust::get<5>(t) = map_constant[ind_deme] + map_coefficient_0[ind_deme] * (thrust::get<1>(t) + thrust::get<2>(t) )/2 + map_coefficient_1[ind_deme] * (thrust::get<3>(t) + thrust::get<4>(t) )/2;
		}   
	};
``` 
Let's walk through what is going on in this code above. The functor will access the parameter value at the individual's deme and apply an additive linear model, based on the underlying genetic values, to determine the phenotype.

<a name="tuples_logic">
</a>

In **sPEGG**, as much as possible we try to have functors operating on more than two variables to operate on what are called tuples. For our purposes right now, you can think of tuples as thrust objects which basically glue together elements from different vectors. For example, suppose you represent each individual's properties with three vectors: **X**, **Y** and **Z**. An individual *i*'s properties can then be specified as a point in this three-dimensional space (X<sub>i</sub>, Y<sub>i</sub>, Z<sub>i</sub>). A tuple allows you to create an ordered structure of these points, which can be envisioned as a list looking something like: [(X<sub>1</sub>, Y<sub>1</sub>, Z<sub>1</sub>), (X<sub>2</sub>, Y<sub>2</sub>, Z<sub>2</sub>), ..., (X<sub>n</sub>, Y<sub>n</sub>, Z<sub>n</sub>)]. As you can imagine, this is very convenient for performing operations on multiple variables that go into describing an individual. For example, suppose we have three loci, whose allelic values we store in vectors **X**, **Y** and **Z**. If we want to evaluate the sum of all three allelic values for individual *i*, we would like a function that can perform operations on the triplet (X<sub>i</sub>, Y<sub>i</sub>, Z<sub>i</sub>). Here, we are going to glue together six vectors: the vector storing the demes of the individuals, two vectors storing the paternally inherited and maternally inherited allelic values of individuals at locus 0, two vectors storing the same data for locus 1, and the vector storing the phenotypic value of each individual at the fecundity phenotype. The elements of this last vector are calculated inside the ```operator()``` function. If, instead of looking at just the additive effects of different loci, we were interested in an epistatic effect, we might use something like:

```c++
void operator()(tuple t) {
		int ind_deme = thrust::get<0>(t);
		thrust::get<5>(t) = map_constant[ind_deme] + map_coefficient_0[ind_deme] * (thrust::get<1>(t) + thrust::get<2>(t)) * (thrust::get<3>(t) + thrust::get<4>(t) );
		}	
```
instead.

We highlight that the thrust functor can also work with data that get passed into the functor as arguments instead of as tuples. In our example functor ```fecundity_calculator()```, these are represented as floating point arrays. This is not ideal, but is often necessary in the current version of Thrust. Indeed, currently Thrust only allows one to have up to 10 vectors "glued" together using tuples. The limitation is problematic. Suppose our gentoype-phenotype map needed to reference 12 genes rather than two. If we're interested in both maternally and paternally inherited allelic values, this would entail at least 24 vectors storing the individual's genotypes (one vector for each gene), and we would be unable to use tuples. In contrast to tuples, you can pass as many vectors as you need to the functor itself as arguments, and reference the elements of the vectors separately. There is some (potentially nontrivial) performance penalty, and there are also more elaborate solutions like having tuples of tuples one could try in the interim. [Thrust will eventually handle more than 10 vectors in a tuple](https://github.com/thrust/thrust/issues/524). Until then, there isn't a single obvious way around this problem.

In any event, we'll call this functor from the method ```calculate_phenotype()``` defined for our ```fecundity_phenotype``` class. We can store the code inside a .cu file, ```Penguins_genotype_phenotype_maps.cu```, inside ```src```. Inside calculate_phenotype(), the functor fecundity_calculator (which remember is a structure) will then get instantiated. We will then "glue" the appropriate vectors together using two thrust operations: ```make_tuple``` and ```make_zip_iterator```. They are described further in the thrust documentation [here](#https://thrust.github.io/doc/group__tuple.html#ga48cf9f1740f033f22386c8c0310c0510) and [here](https://thrust.github.io/doc/group__fancyiterator.html#ga338d4f994660c4dc89e2bb3cf0a6a60f). For now, just know that these are calls to thrust that are necessary to build your tuples. 

Here is the code for calling your functor:
```c++
void fecundity_phenotype::calculate_phenotype(inds *species)
	{
	// wrap the parameters from vectors to arrays:
	float *constants = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT")[0]);
	float *coefficient_0 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);
	float *coefficient_1 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF1")[0]);

	//Instantiate functor
	fecundity_calculator fecundity_calculator_functor(constants, coefficient_0, coefficient_1);

	//Perform genotype-phenotype map operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case, species->fgenotype[0].begin() + index_case, species->mgenotype[0].begin() + index_case, species->fgenotype[1].begin() + index_case, species->mgenotype[1].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case)),
			thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case + num_kids, species->fgenotype[0].begin() + index_case + num_kids, species->mgenotype[0].begin() + index_case + num_kids, species->fgenotype[1].begin() + index_case + num_kids, species->mgenotype[1].begin() + index_case + num_kids, species->phenotype[phenotype_index].begin() + index_case + num_kids)),
			fecundity_calculator_functor);
	}
```
The first three instructions using the raw_pointer_cast() operations are necessary to reference thrust vectors using pointers to arrays that can be passed as arguments distinct from the tuple to our functor. Notice how pointers to these arrays are referenced when we actually instantiate the functor. We then invoke the for_each() function from thrust to perform our genotype-phenotype map calculations on the GPU.

The advantages of having used the strategy pattern become more apparent when we then go on to define the genotype-phenotype map for our next phenotype - a phenotype setting an individual's probability of dying in each time step. Here is the full code for the class definition that could be added to ```Penguins_genotype_phenotype_maps.h```:

```c++
class mortality_phenotype : public GenotypePhenotypeMap
	{
	public:
		mortality_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{};

		void calculate_phenotype(inds *species);
	};


struct mortality_calculator
	{
	float *map_constant; 
	float *map_coefficient_0;

	mortality_calculator(float* map_cons, float* map_coef0) : map_constant(map_cons), map_coefficient_0(map_coef0)
	{};

	/* 
	Elements in the tuple.
	---------------------
	0: individual's deme
	1: individual's maternally inherited allelic value at their locus 0
	2: individual's paternally inherited allelic value at their locus 0
	3: individual's mortality phenotypic value
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int ind_deme = thrust::get<0>(t);
		thrust::get<3>(t) = map_constant[ind_deme] + map_coefficient_0[ind_deme] * (thrust::get<1>(t) + thrust::get<2>(t) )/2 ;
		}   
	};
```	

and the code for performing the calculate_phenotype() inside ```src/Penguins_genotype_phenotype_maps.cu``` can be:

```c++
void mortality_phenotype::calculate_phenotype(inds *species )
	{
	// wrap the parameters from vectors to arrays:
	float *constants = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT")[0]);
	float *coefficient_0 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);

	//Set up functor
	mortality_calculator mortality_calculator_functor(constants, coefficient_0);

	//Perform mortality operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case, species->fgenotype[0].begin() + index_case, species->mgenotype[0].begin() + index_case, species->phenotype[phenotype_index].begin() + index_case)),
					 thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case + num_kids, species->fgenotype[0].begin() + index_case + num_kids, species->mgenotype[0].begin() + index_case + num_kids, species->phenotype[phenotype_index].begin() + index_case + num_kids)),
					 mortality_calculator_functor);
	}
```
Notice that the logic of mortality_phenotype is identical to our fecundity_phenotype except rather than using both loci 0 and 1 to determine the phenotype, we only use one (locus 0). We could have forced further abstractions (e.g., a derived class of ```GenotypePhenotypeMap``` for additively-determined traits) but at this point there is probably little need to in this examle. The advantage of the strategy pattern is that we don't have to modify anything in our code for Penguins to add any other phenotypes. That responsibility gets delegated to the GenotypePhenotypeMap class. This keeps your code self-contained, modular (so you can re-use it for other species or add more phenotypes later, for instance), and makes it easier to test it as you go along.

The final genotype-phenotype map will calculate crown color at birth. The corresponding entries in the header file will be ```Penguins_genotype_phenotype_maps.h```

```c++
class crown_color_phenotype : public GenotypePhenotypeMap
	{
	public:
		int CROWN_COLOR_INDEX;

		crown_color_phenotype(inds *species, int phenotype_index, int index_case, int num_kids) : GenotypePhenotypeMap(species, phenotype_index, index_case, num_kids)
			{
			CROWN_COLOR_INDEX = phenotype_index;
			};

		void calculate_phenotype(inds *species);
	};

struct crown_color_calculator
	{
	float *map_constant; 
	float *map_coefficient_0;

	crown_color_calculator(float* map_cons, float* map_coef0) : map_constant(map_cons), map_coefficient_0(map_coef0)
	{};

	/* 
	Elements in the tuple.
	---------------------
	0: individual's deme
	1: individual's maternally inherited allelic value at their locus 0
	2: individual's paternally inherited allelic value at their locus 0
	3: individual's crown color phenotypic value
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int ind_deme = thrust::get<0>(t);
		thrust::get<3>(t) = map_constant[ind_deme] + map_coefficient_0[ind_deme] * (thrust::get<1>(t) + thrust::get<2>(t) )/2 ;
		}   
	};
```
and in the ```Penguins_genotype_phenotype_maps.cu``` file

```c++
void crown_color_phenotype::calculate_phenotype(inds *species )
	{
	// wrap the parameters from vectors to arrays:
	float *constants = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_CONSTANT")[0]);
	float *coefficient_0 = raw_pointer_cast(&Parameters->get_vector_ptr("GENPHEN_MAP_COEF0")[0]);

	//Set up functor
	crown_color_calculator crown_color_functor(constants, coefficient_0);

	//Perform genotype-phenotype map operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case, species->fgenotype[2].begin() + index_case, species->mgenotype[2].begin() + index_case, species->phenotype[CROWN_COLOR_INDEX].begin() + index_case)),
					 thrust::make_zip_iterator(thrust::make_tuple(species->deme.begin() + index_case + num_kids, species->fgenotype[2].begin() + index_case + num_kids, species->mgenotype[2].begin() + index_case + num_kids, species->phenotype[CROWN_COLOR_INDEX].begin() + index_case + num_kids)),
					 crown_color_functor);
	}
```
will implement the calculations for this genotype-phenotype map using the third locus.

<a name="updating">

### Updating phenotypes
</a>

The basic idea here is to describe how the phenotypes of individuals change between one time step of the simulation to the next. This could happen, for instance, as a result of ontogenic development (e.g., increase in body mass representing growth), or as a result of interactions with the environment (e.g., reduced viability due to parasitism or prolongued drought). As with the case of the genoytpe-phenotype map, again there are many ways we might want to represent how the phenotypes of individuals get updated. Thus, we'll try to re-use the [strategy pattern](https://web.archive.org/web/20160625125844/http://r3dux.org/2011/07/an-example-implementation-of-the-strategy-design-pattern-in-c/) to assign the updating behavior its own class. As with the genotype-phenotype map, this isn't strictly necessary, but we will include it here for illustrative purposes for how this could work for the updating component of our simulation.

As was the case for the genotype-phenotype map, we'll have a base class to manage calling the appropriate updating behavior for our species. We'll call this class UpdateBehavior, and it is defined in ```/header/species/update/updatebehavior.h``` as:

```c++
#ifndef UPDATE_BEHAVIOR_H
#define UPDATE_BEHAVIOR_H

#include <species/inds_stochastic.h>
#include <species/update/survivorship_kernel_functors.h>
#include <environ/environment.h>
#include <thrust/sequence.h>

//virtual interface; 

class UpdateBehavior
	{
	public:
		static UpdateBehavior *create_updateBehavior(inds_stochastic **species, environment *habitat, int species_ID);

		/* the actual updating */
		virtual void update()=0;
		virtual ~UpdateBehavior() {};

		/* Generic functionality for simulating mortality by changing the vital state variable according to a bernoulli RV */
		void determine_mortality(inds_stochastic *species);
	};

#endif
```
<a name="update_behavior_header"></a>

The `create_updateBehavior()` function will then be defined by us in a file we'll call ```UpdateBehavior.cu``` in the ```src/``` of the ```My_sPEGG_project``` folder; a template is provided in ```src/cuda/species/update/UpdateBehavior_TEMPLATE.cu``` of the main **sPEGG** repository.

We'll call this function create_updateBehavior inside our local ```UpdateBehavior.cu``` as:

```c++
#include "update_Penguins.h"

UpdateBehavior * UpdateBehavior::create_updateBehavior(inds_stochastic **species, environment *habitat, int species_ID)
	{
	if (species_ID==0)
	   return new update_Penguins(species[0]);
	}
```
Notice, once again that as in the case of [the genotype-phenotype map above](#gen_phen_map_header), we instantiate a class ```update_Penguins``` that is specific to our tutorial species. Therefore, we need to include an appropriate header file (```update_Penguins.h```) that define this class, so that we can apply the strategy pattern.

To call ```UpdateBehavior``` from the ```Penguins``` class, simply define the ```update()``` member function of ```Penguins``` as:
<a name="update_method">

```c++
void Penguins::update(inds_stochastic **species)
	{
	UpdateBehavior *updatebehavior;
	updatebehavior = updatebehavior->create_updateBehavior(species, NULL, 0);
	updatebehavior->update();
	}
```
</a>

Since we aren't modeling the environmental variables in this tutorial, this is passed as `NULL` when `create_updateBehavior()` is invoked in the `Penguins` class. If we were explicitly modeling the environmental variables and wanted to describe how they affect our phenotypes of interest, we would have to redefine ```void Penguins::update(inds_stochastic **species)``` to instead be ```void Penguins::update(inds_stochastic **species, environment *habitat)```. ```inds_stochastic``` supports the potential use of both ```update``` methods. However, when a new class is created based on ```inds_stochastic```, the update functions are actually hidden from the new class. This is why we had to invoke ```using inds_stochastic::update;``` when we were defining our ```Penguins``` class.  The only remaining variable we worry about here is the species ID, which, because we only have one species, is set to zero. Note that this code can be easily expanded to specify the updates for other species, and the code for ```update()``` in Penguins remains unchanged except for a flag telling create_updateBehavior what species to work with as the focal species.

The first line of the method ```update()``` creates an instance of the UpdateBehavior class. You specify for this class which particular UpdateBehavior you want to simulate. The function ```create_updateBehavior``` will then create an instance of this specification. Finally, the ```update()``` method is called. 

In our example, there are two potential phenotypes that can get updated at each time step. First, we can simulate phenotypic plasticity whereby the crowns "grey" as the penguins age; in our model, we will simulate this greying by allowing the numerical value characterizing crown color to decline to zero. There is one more "phenotype" that can get updated during each time step - the individual's vital status (i.e., whether they are alive = 1 or dead = 0). The ```update_Penguins``` class, therefore, can be defined in a header file ```update_Penguins.h``` in our local ```header``` folder as something like this:

```c++
#include <species/update/updatebehavior.h>
#include <util/amplify.h>

#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>

class update_Penguins : public UpdateBehavior
	{
	class inds_stochastic *species;
	// Constructor
	public:
		update_Penguins(inds_stochastic *species) 
			{
			this->species = species;
	
			// Copy the constants 
			this->size = species->size;
			this->Number_of_Demes = species->Num_Demes;
			}
	
		void update();

	protected:
		void update_color();

		// variable names
		int MORTALITY_PHENOTYPE_INDEX;
		int CROWN_COLOR_INDEX;
		int size;
		int Number_of_Demes;
		void survive();
	};
```

Notice that the class `Penguins.h` is not defined, and instead a place holder for an object of the ```inds_stochastic``` class is specified instead. We present this to illustrate an alternative approach to implementing the strategy pattern like we have here. One potential advantage of the approach we illustrate is that it can reduce the dependence of update_Penguins on the actual penguins class, and only requires knowledge of the inds_stochastic class. This might be helpful if, for instance, we want to use the functionality associated with changing the vital state and the crown color in a simulation of a different bird species altogether.

The ```update()``` method, which gets called by [Penguins](#update_method), serves as the interface to this class that is responsible for updating the penguins. We define its functionality in a .cu file ```update_Penguins.cu``` inside our ```src``` folder. It, too, has a pretty straightforward structure:

```c++
#include "update_Penguins.h"

void update_Penguins::update()
	{
	determine_mortality(species);

	// Change the crown color
	update_color();

	//Increment age for all survivors.
	thrust::transform_if(species->age.begin(), species->age.begin() + size, species->status.begin(), species->age.begin(), unary_plus<int>(1), thrust::identity<int>());
	}
```

You can conceivably add other functions to the ```update()``` method, e.g., something like ```grow()``` that simulates somatic growth. 

The actual updating of the vital status takes place inside the method ```determine_mortality()```. Because the operations described by this method apply quite generally, this function is included as part of the **sPEGG** [code base](/src/cuda/species/update/determine_mortality.cu) as a member of the class UpdateBehavior, rather than being specific to Penguins. The ```determine_mortality()``` method is also potentially complex. We describe it in some depth below. You can skip this for now if you just need a broad overview of things at this stage. 

To implement the ```determine_mortality()``` operation, first we'll describe a functor simulate_mortality() that actually changes the vital states for a specific individual. Analogously to the example of the ```calculate_phenotype()``` method of [the GenotypePhenotypeMap class](#functor_genphenmap), our functor will map an individual's current vital state at time step t-1 to the vital state at time step t (i.e., our functor F() maps a vector of vital states ```x(t-1)=[x1(t-1),x2(t-1),x3(t-1),...,xn(t-1)]``` to a vector ```x(t) = [F(x1), F(x2), F(x3),...,F(xn)]```). Here's how our functor actually does this: it takes as its argument a vector of *n* uniformly distributed random variables ```[U1, U2, ..., Un]```. For each individual in i=[1,...,n] if their corresponding random number is less than their probability of survivorship, their vital state is set to 1 (i.e., they survive). Otherwise, the vital state is set to 0 (i.e., they die).

Here is what the functor that actually modifies the vital state variable looks like.

```c++
// Functor that calculates an individual's fate.
struct simulate_mortality
	{
	/* 
	* Functor related to survivorship: to be invoked in the update_mySpecies class's determine_mortality() function, specified in \loc cuda/species/update/survivorship_kernel.cu.
	*/

	float *uniform_rv;

	simulate_mortality(float* uniform_random) : uniform_rv(uniform_random)
	{};

	/* 
		Elements in the tuple.


		----------------------
		0: individual index
		1: individual's vital state
		2: probability of survivorship

	*/ 
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int vital_state = thrust::get<1>(t);

		if (vital_state==1) /* If the individual is alive */
			{
			float unif = uniform_rv[thrust::get<0>(t)];
			float surv_val = thrust::get<2>(t);
			
			// Actually kill them
			if (unif < surv_val)
				vital_state = 1;
			else
				vital_state = 0;
			thrust::get<1>(t) = vital_state;
			}
		}	
	};
```
Notice that we've used [tuples](#tuples_logic) again and access them through the ```thrust::get``` operation. Our tuple consists of three vectors: the vector storing the indices of the individuals (i.e,. [0,1,2,..,n]; these are not the same as the individual's IDs), the vector storing their vital states, and the vector storing their probability of survivorship (which is equal for everyone). Technically, we could have also included the vector of uniform random variates in this list, but we chose to leave that as an argument to illustrate how you can also pass the vectors as arguments to your functors.

Next, we'll actually invoke this functor in the function "determine_mortality()" described above. This method will first create the indices of individuals. Next, it will draw a uniform random number for each individual. The functor will then get instantiated and called, using the "glued" vectors that contain the indices, the vital states, and the probabilities of survivorship. 

<a name="survive_kernel">
</a>

```c++
void UpdateBehavior::determine_mortality(inds_stochastic *species)
	{
	//Specify the individuals indices
	thrust::device_vector<int> individuals(species->size);
	thrust::sequence(individuals.begin(), individuals.begin() + species->size, 0);

	// Draw the random numbers
	thrust::device_vector<float> rand(species->size);
	// wrap the vector
	float *rand_ptr = raw_pointer_cast(&rand[0]);
	curandGenerateUniform(species->gen, rand_ptr, species->size);

	//Set up mortality.
	simulate_mortality mortality_functor(rand_ptr);

	//Perform mortality operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(individuals.begin(), species->status.begin(), species->phenotype[species->MORTALITY_PHENOTYPE_INDEX].begin())),
			 thrust::make_zip_iterator(thrust::make_tuple(individuals.begin() + species->size, species->status.begin() + species->size, species->phenotype[species->MORTALITY_PHENOTYPE_INDEX].begin() + species->size)),
			 mortality_functor);
	}
```

Admittedly, this ```determine_mortaility()``` function is a lot more complicated, particularly when we consider also the associated functor, than what an equivalent pure C++ version which could look like:

```c++
void UpdateBehavior::determine_mortaility()
	{	
	for (int i=0; i < species->size; i++)
		   {
		   u = draw_uniform(0,1);
		   if (u < species->phenotype[species->MORTALITY_PHENOTYPE_INDEX][i])
					  species->status[i] = 1;
		   else
					  species->status[i] = 0;
		   }
	}
```

It turns out for this particular problem, the performance gains for moving to the GPU are modest. But as we argue in our [accompanying manuscript](/Documentation/1603.09255v1.pdf), once one works with something beyond relatively trivial models, the speed-up provided by modeling updating of the phenotype using CUDA can be considerable.

<a name="update_color">
</a>

In addition to simulating the mortality of individuals, our ```update()``` kernel includes a function for updating the color of the crown for each individual to simulate greying as they age. Our approach to implementing this functionality is very similar to our approach to implementing mortality. First, we will define a method to simulate this greying inside the class ```update_Penguins```. Then, we will specify a functor to offload the actual calcuations across individuals on the GPU. To do this, we begin by specifying a functor, ```crown_color_updater``` that will actually change the color of the crown for each individual. Because this is a functor, we will again declare it as a struct inside the header file associated with the updating operations. We note, however, that this operation is modeled as being specific to simulated penguins (unlike updating mortality, which we assumed was generic to any inds_stochastic object). We therefore include the functor ```crown_color_updater``` in the same ```update_Penguins.h``` header where the object ```update_Penguins``` is defined. Here is the functor:

```c++
struct crown_color_updater
	{
	float *crown_phenotypes;
	float *target_crown_color;

	crown_color_updater(float* crown_phens, float* targ_crown_color) : crown_phenotypes(crown_phens), target_crown_color(targ_crown_color)
	{};

	/* 
	Elements in the tuple.
	---------------------
	0: individual's index
	1: individual's deme
	2: the rate at which crown color decays for the individual's deme.
	*/ 

	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) 
		{
		int ind_index = thrust::get<0>(t);
		int ind_deme = thrust::get<1>(t);
		float delta_crown_color = thrust::get<2>(t);

		float old_crown_color = crown_phenotypes[ind_index];
		crown_phenotypes[ind_index] = old_crown_color - delta_crown_color*(old_crown_color - target_crown_color[ind_deme]);
		}   
	};
```
The functor takes as its arguments two arrays: an array describing the current crown color for each individual, and an array describing the target crown color to which each individual's crown color will decay. This latter target color is assumed fixed for the entire deme. We emphasize that we could have made both of these part of our tuple. However, we include them as vector arguments to illustrate how one would handle situations where more than 10 vectors needs to be fed into a functor. Moreover, because the target crown color does not vary across individuals within a given deme, the length of the array storing the target crown color is different than the length of the crown_phenotypes array. Therefore, forcing both arrays into a single tuple would require additional operations (see below for the rate of crown color decay) that we can avoid by feeding them as array parameters for our functor rather than tuples. 

The workhorse ```operator()``` of our functor takes only the bare minimum information for each individual: its index, its deme, and the rate at which the crown color decays. It then calculates the new crown color for this individual based on a simple arithemtic operation. The rate governing the change in crown color is also deme-, rather than individual, specific. As such, it is based on an array whose length differs from the array storing the individual indices and their demes. We will illustrate below how we make sure these arrays have the same length before we glue them together into a tuple.

To invoke the functor ```crown_color_updater```, we use the method ```update_color()```, defined in ```update_Penguins.h``` (because this is a penguin-specific operation). The function's behavior has a broad outline very similar to the method ```determine_mortality()```. Here is the function:

```c++
void update_Penguins::update_color()
	{
	int CROWN_COLOR_INDEX = species->demeParameters->species_specific_values["CROWN_COLOR_INDEX"];

	// Step 1. Cast the thrust vectors into arrays 
	float *crown_colors = raw_pointer_cast(&species->phenotype[CROWN_COLOR_INDEX][0]);
	float *target_crown_colors = raw_pointer_cast(&species->demeParameters->get_vector_ptr("TARGET_CROWN_COLOR")[0]);
	
	// Step 2. Set up functor
	crown_color_updater crown_color_functor(crown_colors, target_crown_colors);

	// Step 3. Create a vector which, for each individual, stores the crown color decay rate associated with their deme.
			// Step 3a. Create a one-time vector, demewise_color_decays that copies the deme-specific color decay rates.
	thrust::device_vector<float> demewise_color_decays(species->Num_Demes);
	thrust::copy(species->demeParameters->get_vector_ptr("CROWN_COLOR_DECAY"), species->demeParameters->get_vector_ptr("CROWN_COLOR_DECAY") + species->Num_Demes, demewise_color_decays.begin());

			// Step 3b. Create a vector, color_decays of length size (i.e., the number of living individuals) which stores, for each individual, the crown color decay rate associated with their deme.
	thrust::device_vector<float> color_decays(size);
	amplify_float(demewise_color_decays, species->deme_sizes, color_decays);

	// Step 4. specify the indices of individuals
	thrust::device_vector<int> individuals(size);
	thrust::sequence(individuals.begin(), individuals.begin() + size, 0);

	// Step 5. Perform greying operation with for_each.
	thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(individuals.begin(), species->deme.begin(), color_decays.begin())),
					 thrust::make_zip_iterator(thrust::make_tuple(individuals.begin() + size, species->deme.begin() + size, color_decays.begin() + size)),
					 crown_color_functor);
	}
```
Casting the thrust vectors storing the individual crown colors to an array pointer ```crown_colors```, as we do in step (1), should be familiar at this point. The ```demeParameters``` object contains a vector that stores, for each deme, the target crown color for that deme. Thus, this line simply casts this vector into an array to be used when we instantiate our functor in Step 2. 

Recall that in our functor definition above, we include the crown color decay rate as a tuple element. However, the crown color decay rate is associated with a deme, rather than an individual. As such, the vector storing the crown color decay rates in demeParameters is not going to necessarily have the same length as the array ```crown_colors``` storing the crown colors for each individual. This is problematic, as our functor requires calculating, for each individual, their new color based on their specific crown color decay rate. Thus, in Step 3 we use a generic method ```amplify_float``` that comes with **sPEGG** to create a new vector, ```color_decays``` of length equal to the number of individuals. The *i*th element of this vector corresponds to the crown color decay rate of the *i*th individual. 

A worked example could help illustrate what has to happen to create this vector. Suppose we are simulating a total of four individuals with indices [0,1,2,3]. Suppose Individuals [0,1,2] belong to deme 0, and individual [3] belongs to deme 1. Then the corresponding vector of demes (stored in the vector ```species->deme```) is [0,0,0,1]. Let's pretend we only have two demes, and the color decay rate associated with deme 0 is 0.5, while the rate associated with deme 1 is 0.25. To assign every individual the decay rate associated with their deme-wide, we need a vector that looks like ```color_decays=```[0.5,0.5,0.5,0.25] in order to correctly carry out the operation specified by the functor ```crown_color_updater```. The function amplify_float in step 3b creates the vector ```color_decays``` from the vectors [0.5,0.25] and [3,1]; that is, it creates a vector that repeats the first element 3 times, and the second element once. 

Step 4 then simply creates a vector storing the lists of individuals indices, and the functor is actually evaluated in step 5.

Once the survivors have been determined and ```update_color()``` has been called, the ```update()``` [kernel](#update_method) ends by incrementing the age of survivors by one. 

<a name="helper_methods">

### Helper methods
</a>

Hopefully, by now you've gotten a bit of a sense for how thrust will work with the data elements in **sPEGG**. In particular, you have probably gotten a good feel for the centrality of an individual's deme affiliation to the simulation. Earlier, we glossed over the helper function ```initialize_demes()``` when the ```Penguins``` object is initialized. Discussion of this function was relegated to here so the code above can be made somewhat more concise without getting too much into the book-keeping details.

We'll write a short function to place individuals into demes at the initialization stage. In more realistic models, you might want to allocate the starting deme sizes to have different numbers of individuals, perhaps reflective of spatial variation in starting population sizes. Here, we will assume that each deme has the same carrying capacities and the same starting number of individuals.
```c++
void Penguins::initialize_demes(int ndemes)
	{
	thrust::fill(max_deme_sizes.begin(), max_deme_sizes.begin() + Num_Demes, maxsize/Num_Demes);

	for (int i=0; i < Num_Demes; i++)
		{
		float startsize = (float) (size/Num_Demes);
		int temp1 = (int) startsize * i;
		int temp2 = (int) startsize * (i+1);
		thrust::fill(deme.begin() + temp1, deme.begin() + temp2, i);
		}
	}
```
A final helper function, assignSex, will randomly assign sex to new borns. We will create a temporary vector, neonates_sex, that will temporarily store a series of bernoulli random variates. We will then call a built-in **sPEGG** function, draw_bernoulli, that takes 4 arguments: the number of random variates to be simulated, the probability of success, the vector that will store the random variates, and a random number generator. Finally, we will copy the elements of the temporary vector to the corresponding elements of the vector that stores each individual's sex.

```c++
void Penguins::assignSex(int index, int num_inds_to_calculate)
	{
	thrust::device_vector<int> neonates_sex(num_inds_to_calculate);
	draw_bernoulli(num_inds_to_calculate, 0.5, neonates_sex, gen);
	thrust::copy(neonates_sex.begin(), neonates_sex.begin() + num_inds_to_calculate, sex.begin() + index);
	}
```

The resulting Penguins.cu file that integrates everything we've covered so far can be downloaded [here](Penguins.cu).

---------------------------------

# Putting it all together

## The simulation class

Now that you've defined your Penguins class and its associated methods, you'll want to actually program the simulation with it. To do this, we'll create a ```Simulation``` object we'll call ```Penguin_Drift_Simulator``` that will run our simulation for us. This can be done using the pre-defined Simulation class that comes with **sPEGG**. This file can be found by navigating to the ```header``` directory of the main **sPEGG** code base.

We will only need to use a handful of headers to define our simulation class. We'll obviously need the ```Simulation_Class``` definition and we'll also want the ```Penguins``` class definition. Finally, in our simulation we would be interested in tracking what goes on both at the individual level (e.g., individual genotype values), as well as at the deme-level (e.g., mean phenotypic value). For this we will rely on the ```Statistics``` class, also provided by **sPEGG**, that will allow us querry our penguin class to get summary statistics during the simulation. The statistics class does this by interfacing with ```inds``` objects to summarize the individual-level data into deme-level data. Finally, we could use an array of inds objects to store potentially different penguin species. Here is what our header file ```Penguin_Drift_Simulator.h``` inside our local ```header``` directory could look like:

```C++
#ifndef PENGUIN_SIMULATOR_H
#define PENGUIN_SIMULATOR_H

#include <Simulation_Class.h>
#include "Penguins.h"
#include <math/statistics_class.h>

class Penguin_Drift_Simulator : public Simulation
	{
	public:
		Penguin_Drift_Simulator();
		~Penguin_Drift_Simulator();
		void run();
	private:
		inds_stochastic **array;
		Statistics *stats_penguins;
		void initialize_classes();

		int nspecies;
	};
#endif
```

Some comments on the class methods. Both the ```Penguins``` class and the statistics class ```stats_penguins``` will need initializing (in the case of the penguins class, we need to specify initial conditions), which is what the method ```initialize_classes()``` will accomplish. Our simulator will run when we tell it to via the function ```run()```, and the accompanying classes will be freed from GPU memory when the destructor to ```Penguin_Drift_Simulator()``` is called. Let's look at these methods a little more closely by defining them within a corresponding ```.cpp``` file ```Penguin_Drift_Simulator.cpp``` inside the ```src``` folder of ```My_sPEGG_project```. Notice here we can use a ```.cpp``` file extension here because all the cuda material is actually used in other classes, and so we can go back to using basic C++.

The method ```initialize_classes()``` will create an array of ```inds_stochastic``` objects, which might or might not represent penguins. But let's assume the first element of this array is indeed a ```Penguins``` object representing individuals in ```Penguins```. We'll just call the constructor for ```Penguins```, as well as the constructor for the ```Statistics``` class. Its constructor takes a single argument specifying how many demes we are simulating. Finally, once initialization is done, we'll store our new penguins data in an external CSV file we'll call initial_data.csv. Here is the code:

```C++
void Penguin_Drift_Simulator::initialize_classes()
	{
	array = new inds_stochastic *[nspecies];
	int species_ID = 0;
	array[0] = new Penguins(initpop, maxpop, seed, demes, species_ID);

	stats_penguins = new Statistics(demes);

	array[0]-> exportCsv("initial_data.csv");
	}
```
 
Just as we had to initialize our ```Penguins``` class, we will also need to initialize the Penguin_Drift_Simulator() using a constructor. For our purposes, let's suppose we'll only model one species, and we'll seed our two demes with 50 individuals each. We'll also arbitrarily set a maximum number of individuals at a 1000 individuals per deme. The constructor therefore looks like:

```C++
Penguin_Drift_Simulator::Penguin_Drift_Simulator() : Simulation()
	{
	initpop =50*demes;
	maxpop = 1000*demes;

	nspecies = 1;
	
	initialize_classes();
	}
```

The actual simulation can be described inside the ```run()``` method. This will iterate the simulation for ```nsteps``` generations for each species, during which reproduction, updating (i.e., modifying individual's phenotypes) and removal of the dead will take place. When the simulation is done, the mean genotypic values at, say, locus 2, which codes for crown color at birth, for each deme will be calculated and printed. Finally, we'll store the final state of our simulation for the first species of penguins (with all the information about living individuals) in a CSV file we'll call final_data.csv. 

```C++
void Penguin_Drift_Simulator::run()
	{
	for (int t = 0 ; t < nsteps ; t++) 
		{
		for (int i = 0 ; i < nspecies ; i++) 
			{
			array[i]->addKids();
			array[i]->update(array);
			array[i]->removeDead();
			}
		}

	int genotype_index_of_interest = 2;

	stats_penguins->calculate_mean_genotypes_by_deme(array[0], genotype_index_of_interest);
	std::cout << "The average genotypes at locus 1 for the two demes are:" << std::endl;
	stats_penguins->print_mean_genotypes_by_deme();

	array[0]-> exportCsv("final_data.csv");
	}
```
Notice we're only interested in the summary statistic and CSV file for the single species of ```Penguins```. So we only refer to the pointer ```array[0]``` here.

Finally, once our simulation is done, we will want to clear our objects ```array``` and ```stats_penguins``` from the GPU memory. To do this, we'll define our destructor as follows:

```C++
Penguin_Drift_Simulator::~Penguin_Drift_Simulator()
	{
	/* cleanup */
	for (int i=0; i < nspecies; i++)
		{
		delete array[i];
		}

	delete[] array;
	delete stats_penguins;
	}
```


<a name="main_cpp_fixing">

## The main.cpp file

To actually run our simulation, we'll create a C++ file in the main directory (```My_sPEGG_project```). There, create a file called ```main.cpp```. In addition to the standard C++ headings, a header for the ```Penguin_Drift_Simulator``` class will be referenced. We will instantiate our simulation, run it, and free the memory on the GPU. Here is the code:

```C++
#include "Penguin_Drift_Simulator.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>

int main(void)
	{
	Penguin_Drift_Simulator *Penguin_model;
	Penguin_model = new Penguin_Drift_Simulator();
	Penguin_model->run();
	delete Penguin_model;
	}
```
</a>

In the main.cpp file, we'll have only a single function, ```int main(void)```. By initializing the ```Penguing_Drift_Simulator``` object, our ```main()``` function sets in motion the following sequence of activities:

>    1. read in the basic simulation settings, 
	
>    2. initialize our classes, 
	
>    3. save the individual-level data at the start of the simulation into a CSV spreadsheet

Then, when it calls the ```run()``` function, the following happen:
	
>    4. run the simulation,

>    5. calculate the average genotypes at locus 2 for the two demes at the end of the simulation,

>    6. save the individual-level data for the single penguin species at the end of the simulation into a CSV spreadsheet

The ```delete Penguin_model``` step ensures we:

>    7. free our objects from memory.

We're almost done. Now you can setup your configuration (text) files specifying how you want your simulation to run. We'll start with Simulation.conf, which controls how the arguments like seed, nloci etc... get set.

## Inputting parameters

Ideally, once you have a **sPEGG** simulation coded for the problem you are trying to model, you want to minimize having to modify C++ code to change things too much. Moreover, you want your model to be able to interface with data, which means that at some point, the focus on the analysis moves from model development to comparing model input (numerical parameter values measured in the laboratory or the field) to output (model-generated v. observed data). **sPEGG** tries to facilitate your efforts to minimize coding, and aims to allow you to focus on getting the inputs into the model with as little hassle as possible once the hard work of coding is finished. 

To that end, each model implemented in **sPEGG** has associated with it four input files:

  * a very simple Simulation.conf file that specifies basic model parameters like how many demes are simulated, how long the simulation is to run, etc..., 
  * a file that configures the species and their demes (deme_config.txt),
  * an optional file that configures the environment in which the species will be interacting (environment_config.txt) and
  * an optional file that specifies the migration rate between demes. 

Below, we briefly describe each of these input files, which you can create right now in the folder ```My_sPEGG_project```.

### Simulation.conf

This file is deliberately kept short and limited. It configures the simulation by specifying (i) the number of time steps (n_timesteps), which can be thought of as corresponding to the number of reproductive bouts, (ii) the number of demes (ndemes), (iii) the number of time steps simulated between reproductive events within a time step, if any, (intra_step_time_steps - e.g., months of a year for an annually reproducing species; this might also allow you to build a model that approximates continuous-time processes between reproductive events), (iv) the number of biotic variables (num_biotic_variables) in the environment, (v) the number of abiotic variables (num_abiotic_variables) in the environment, and (vi) the starting random number seed. For this tutorial, our Simulation.conf file looks something like:

```sh
n_timesteps = 2
ndemes = 2
intra_step_time_steps = 0
num_biotic_variables = 0
num_abiotic_variables = 0
random_seed = 314159
```

<a name="demesettings">

### deme_config.txt and environment_config.txt
</a>

Since we're simulating genetic drift and our loci are neutral markers, there might not seem to be a need for specifying differences among the demes. However, our code still makes use of several parameter inputs, such as the coefficients defining the allelic effects in the genotype-phenotype map and the decay rate for crown color. <a name="demeParameters"> This input gets stored in an object called ```demeParameters```</a> included in the ```inds``` base class. So although some of the following code might seem a bit redundant/unnecessary, it's in this tutorial because these same principles will apply when you want to use **sPEGG** to model a more conceptually interesting or realistic problem.

#### the deme_config.txt file

All parameters pertaining to your species are stored in this file. Input from the .txt file is handled by the [libconfig library](http://www.hyperrealm.com/libconfig/). Here is a pseudo-code skeleton strucuture of a deme_config.txt file

```sh
number_of_demes = n
species_data = ({inputs for species 0}, {inputs for species 1}, ..., {inputs for species n})
species_genetics = ({genetic inputs for species 0}, {genetic inputs for species 1}, ..., {genetic inputs for species n})
```
At this point, it would be good to review some important features of libconfig. libconfig distinguishes betwen ```lists```, ```arrays``` and ```groups```. Lists are demarcated by parentheses (), arrays are demarcated by brackets [], whilst groups are demarcated by curly brackets {}. All elements of arrays must be the same type (e.g., all integers, all characters, etc...). The main difference between the lists and arrays and groups is that in libconfig, groups contain one or more of what libconfig calls ```settings```. ```settings``` are any entities within the config file that have the form ```name = value```. For our purposes, ```value``` consists only of groups, lists and single numeric, boolean or character values (e.g., ```3.14``` or ```tree_height```). Thus, in our skeleton above, ```number_of_demes```, ```species_data``` and ```species_genetics``` are all settings. These points and others are explained very well in [libconfig's documentation](http://www.hyperrealm.com/libconfig/libconfig_manual.html#Configuration-Files). Below, we'll discuss each of these settings in turn.

The parameter number_of_demes is straightforward. It stores an integer - which should match the value in Simulation.conf (later versions might dispense of it because of this potential redundancy) - that tells us how many demes we are simulating.

Each of the ```{inputs for species *x*}``` within the list ```species_data``` specifies, for each species:

  * the number of non-parametric variables (which we will call ```species_specific_values```) that apply across all demes
  * the names of the non-parametric variables
  * the values of these these non-parametric variables.
  * the number of parameters governing the simulated dynamics of this species that could potentially differ across demes
  * the names of these parameters
  * for each deme, the value of those parameter values
  
The distinction between ```species_specific_values``` and parameters in **sPEGG** are that the former (1) will never vary across demes, and (2) help manage the code rather than modify the trajectory of the simulation. An example will probably help clarify the distinction.

In our example, we only have one species - our ```Penguins```. Therefore, our species_data setting should look something like ```species_data = ({inputs for demes of species 0})```. Here is how we specify the inputs for demes:

```sh
species_data = ({
				number_of_species_specific_values = 3
				species_specific_values_names = ["FECUNDITY_PHENOTYPE_INDEX","MORTALITY_PHENOTYPE_INDEX","CROWN_COLOR_INDEX"]

				FECUNDITY_PHENOTYPE_INDEX  =  0.000000
				MORTALITY_PHENOTYPE_INDEX  =  1.000000 
				CROWN_COLOR_INDEX = 2.00000

				number_of_parameters = 4
				parameter_names = ["M_reproductive_advantage","F_reproductive_advantage","TARGET_CROWN_COLOR","CROWN_COLOR_DECAY"]

				deme_specifications = (
										{#deme0
										M_reproductive_advantage = 1.0
										F_reproductive_advantage = 1.0
										TARGET_CROWN_COLOR = 0.808080
										CROWN_COLOR_DECAY = 0.1
										},
										{#deme1
										M_reproductive_advantage = 1.0
										F_reproductive_advantage = 1.0
										TARGET_CROWN_COLOR = 0.809120
										CROWN_COLOR_DECAY = 0.01
										}
									  )               
				}
			   )
```

In our example, we consider two parameters: (1) a parameter controlling the reproductive skew towards males with higher fecundity (i.e,. if *y[i]* is a male *i*'s fecundity, the parameter ```M_reproductive_advantage``` controls the relative probability that a male will sire an offspring according to *y[i]^(M_reproductive_advantage)*), and (2) a parameter controlling the reproductive skew towards females with higher fecundity. For our tutorial, we'll assume these values are identical across the two demes; however, you can imagine a situation where there is greater reproductive inequality among individuals in deme 0 than in deme 1 - e.g., if you want to explore a gradient of this value to see whether varying reproductive skew affects dynamics. Under these cases, your deme_specifications list may look instead like:


```sh
	   deme_specifications = (
										{#deme0
										M_reproductive_advantage = 1.0
										F_reproductive_advantage = 2.5
										TARGET_CROWN_COLOR = 0.808080
										CROWN_COLOR_DECAY = 0.1
										},
										{#deme1
										M_reproductive_advantage = 1.0
										F_reproductive_advantage = 1.0
										TARGET_CROWN_COLOR = 0.809120
										CROWN_COLOR_DECAY = 0.01
										}
									  ) 
```

We should add at this point that libconfig doesn't care about whitespace, and comments can be added by placing ```#``` right before the comment (thus, the phrases "deme0" and "deme1" above are ignored at run-time).

The ```species_genetics``` list is defined in much the same way. For each species you model, you will specify:

  * the number of loci modeled
  * (optionally) the names of the loci - e.g., COX1
  * the number of phenotypes modeled
  * the names of the phenotypes
  * the recombination rates among loci, which we are going to assume is constant across demes for our species
  * for each locus, parameters controlling its mutation rate and mutation magnitude (which could vary by deme)
  * for each phenotype, the genotype-phenotype parameter values by each deme.

Here's a basic outline, using pseudo-code
```sh
species_genetics = ({#species0
				number_of_loci = nL
				number_of_phenotypes = nP 
				phenotype_names = ["PHEN1","PHEN2"]

				recombination_rates = [recombination rates]

				locus_specifications = (
					{#locus0
					mutation_parameters = (
								{#deme0
								MUTATION_RATE = 0.0001
								MUTATION_MAGNITUDE = 0.1
								},
								{#deme1
								MUTATION_RATE = 0.0001
								MUTATION_MAGNITUDE = 0.2
								}
								)
					},
					{#locus1
					mutation_parameters = (
								{#deme0
								MUTATION_RATE = 0.0001
								MUTATION_MAGNITUDE = 0.3
								},
								{#deme1
								MUTATION_RATE = 0.0001
								MUTATION_MAGNITUDE = 0.4
								}
								)
					}, ...
					)
				phenotype_specifications = (
										{#phenotype0
										number_of_genotype_phenotype_map_parameters = nP0
										names_of_genotype_phenotype_map_parameters = [parameter names]
										genotype_phenotype_map_parameters = ({parameter values in deme 0},{parameter values in deme 1})
										},
									   {#phenotype1
										number_of_genotype_phenotype_map_parameters = nP1
										names_of_genotype_phenotype_map_parameters = [parameter names]
										genotype_phenotype_map_parameters = ({parameter values in deme 0},{parameter values in deme 1})
										}, ...
									  )               
				}
			   )
```
Most of the terms here are pretty self-explanatory. The basic logic is that only the numerical parameters governing the genotype-phenotype map, which are stored in a list called ```genotype_phenotype_map_parameters```, would vary by deme (for instance, if you want to model how different genotype-phenotype maps - i.e., genetic architectures - would affect your simulation's results, you'd simulate multiple demes to look at how varying the constituent parameters affect the evolutionary trajectory). The rest of the genetic architecture is assumed to be constant across all demes, as is the recombination rates among loci. Potentially the recombination rate among loci can be modified according to deme, and this is a candidate feature for future releases. Note also that you can render some loci unimportant in specific demes by setting the parameter corresponding to their allelic effects to zero.

Our resulting genetic specification (refer to [the section on the genetic architectures of our model's phenotypes](#gen_phen_map) for the genotype-phenotype map) looks something like this:

```sh

species_genetics = ({#species0
				number_of_loci = 3
			   
				number_of_phenotypes = 3 
		phenotype_names = [ "FECUNDITY_PHENOTYPE_INDEX", "MORTALITY_PHENOTYPE_INDEX", "CROWN_COLOR_INDEX" ]
		recombination_rates = [ 0.5, 0.5, 0.5 ]

		locus_specifications = (
					{#locus0
					mutation_parameters = (
								{#deme0
								MUTATION_RATE = 0.0001
								MUTATION_MAGNITUDE = 0.1
								},
								{#deme0
								MUTATION_RATE = 0.0001
								MUTATION_MAGNITUDE = 0.2
								}
								)
					},
					{#locus1
					mutation_parameters = (
								{#deme0
								MUTATION_RATE = 0.0001
								MUTATION_MAGNITUDE = 0.3
								},
								{#deme0
								MUTATION_RATE = 0.0001
								MUTATION_MAGNITUDE = 0.4
								}
								)
					},
					{#locus2
					mutation_parameters = (
								{#deme0
								MUTATION_RATE = 0.0001
								MUTATION_MAGNITUDE = 0.5
								},
								{#deme0
								MUTATION_RATE = 0.0001
								MUTATION_MAGNITUDE = 0.6
								}
								)
					}
					)

				phenotype_specifications = (
										{#phenotype0
					number_of_genotype_phenotype_map_parameters = 3
					names_of_genotype_phenotype_map_parameters = ["GENPHEN_MAP_CONSTANT", "GENPHEN_MAP_COEF0","GENPHEN_MAP_COEF1"]
					genotype_phenotype_map_parameters  = (
															{#deme0
										GENPHEN_MAP_CONSTANT = 10.0
															GENPHEN_MAP_COEF0 = 0.0
															GENPHEN_MAP_COEF1 = 0.0
															}, 
															{#deme1
										GENPHEN_MAP_CONSTANT = 10.0
															GENPHEN_MAP_COEF0 = 0.0
															GENPHEN_MAP_COEF1 = 0.0
															}
														  )
										},
										{#phenotype1
					number_of_genotype_phenotype_map_parameters = 2
					names_of_genotype_phenotype_map_parameters = ["GENPHEN_MAP_CONSTANT", "GENPHEN_MAP_COEF0"]
					genotype_phenotype_map_parameters  = (
															{#deme0
								GENPHEN_MAP_CONSTANT = 0.5
															GENPHEN_MAP_COEF0 = 0.0
															}, 
															{#deme1
								GENPHEN_MAP_CONSTANT = 0.5
															GENPHEN_MAP_COEF0 = 0.0
															}
														  )
										},
					{#phenotype2
					number_of_genotype_phenotype_map_parameters = 2
					names_of_genotype_phenotype_map_parameters = ["GENPHEN_MAP_CONSTANT", "GENPHEN_MAP_COEF0"]
					genotype_phenotype_map_parameters  = (
															{#deme0
								GENPHEN_MAP_CONSTANT = 0.25
															GENPHEN_MAP_COEF0 = 3.0
															}, 
															{#deme1
								GENPHEN_MAP_CONSTANT = 0.25
															GENPHEN_MAP_COEF0 = 3.0
															}
														  )
										}
									  )               
				}
			   )
```
Since we don't assume there are genetic differences among individuals in fitness or fecundity, their coefficients governing how they map onto phenotypes controlling fitness are zero, and our genotype-phenotype map looks something like (*F(x) = A*, where *A* is the value in the variable ```GENPHEN_MAP_CONSTANT```). 

#### the environment_config.txt file

The environment_config.txt file works much the same way. You'll specify the number of demes, but instead of lists containing species-specific parameters and genetic details, you'll have separate lists storing parameters related to abiotic and biotic variables. These lists will consist of groups of parameter values by deme. Our model doesn't have environmental variables that feedback to our evolutionary dynamics, but a minimal environment_config.txt would look something like:

```sh
number_of_demes = 2
abiotic_variable_names = ["dummy"]

abiotic_variable_initialization = 
	(	{#deme0
		 dummy=0.0
		},
		{#deme1
		 dummy=1.0
		}
	)
```

If there is sufficient interest, a later tutorial module might go into some detail about how the environmental components that aren't modeled on an individual level would be simulated. In the mean time, to get a feel for what the environment_config.txt for a model with resource dynamics (where the resource - e.g., prey - is not tracked at the individual level) might look like, see [this example](/Examples/PhysiologicalStructure_LifeHistoryEvolution/environment_config.txt).

### Using demes to represent replicates
It's worth reiterating at this stage that different demes can also represent different replicate simulations. In such cases, their settings as specified in these configuration files would be identical. You can also vary demes across a gradient of a parameter value, and have replicate demes for each level of the parameter's gradient. In **sPEGG**, "demes" thus represent a versatile concept that could differentially represent real subpopulations in a metapopulation, as well as a convenient way of organizing your simulated replicates and treatment units.

--------------------------

# Putting it all together
That pretty much covers the basics of what you need to build your own **sPEGG** simulation. A few final points on getting this to run on your computer.

<a name="makefile_fixing">
</a>

## The Makefile

In the directory <My_sPEGG_project> you will also create a file called ```Makefile```. This can be somewhat involved, but one thing you may have to change is the directory in which the main sPEGG codebase was saved and built. As you add source files to your project, you will want to keep adding compilation instructions for those files to your Makefile - i.e.:

```Makefile
(OBJDIR)/My_New_Source_File.o : ${LOCATION_OF_NEW_SOURCEFILE}/My_New_Source_File.cu
	nvcc -c $(CFLAGS) ${HEADERS}  ${LOCATION_OF_NEW_SOURCEFILE}/My_New_Source_File.cu -o $(OBJDIR)/My_New_Source_File.o
```

Also add instructions to append the resulting *.o file to the line beginning

```Makefile
a.out :
```

and the line containing the compilation instructions:

```sh
nvcc -O3 -lcurand -lrt -lcuda -lconfig++ 
```

by preceding the .o file with the tag $(OBJDIR)/ e.g., 

```sh
a.out : $(OBJDIR)/My_New_Source_File.o
nvcc -O3 -lcurand -lrt -lcuda -lconfig++ $(OBJDIR)/My_New_Source_File.o
```

A Makefile for this tutorial can be found [here](/Examples/Tutorial_Simulation/Makefile).

When your Makefile is complete, type ```make``` in the terminal, and you can run your simulation via:

```sh
./a.out
```

And you are now good to go! You can get the full source-code used in this tutorial [here](/Examples/Tutorial_Simulation).

Of course, how you organize the **sPEGG** project is entirely up to you, and not all of the solutions we identified would be ideal for your model. The approach here only describes the workflow I've found works for me as I tried to simulate a range of different problems :grin:.

## Exercises
1. One parameter that might be of interest may be to vary the demographic rates across demes. Let's assume that deme #0 permits a individuals to live longer (on average, one generation) than individuals in deme #1. For now, we'll keep it simple and won't model the mechanism behind this difference. How would you model this situation without changing your underlying C++/CUDA thrust code? (hint: look at the GENPHENO_MAP_CONSTANT in your deme_config.txt)

A more elaborate (non-stochastic) way to simulate of this might involve defining a functor describing death due to senescence as: 

```c++
struct death_from_old_age
	{
	float *max_age;
	death_from_old_age(float *maximum_age) : max_age(maximum_age)
	{};
	/* 
		Elements in the tuple.
		---------------------
		0: individual's index
		1: individual's vitality status
		2: probability of survivorship
		3: individual's deme
		4: individual's age
	*/ 
	
	template <typename tuple>
	__host__ __device__
	void operator()(tuple t) {
		if (thrust::get<1>(t)==1) /* If the individual is alive */
			{
			int deme = thrust::get<3>(t);
			float age = (float) thrust::get<4>(t); 
			if (age > max_age[deme])
				 {
				 thrust::get<2>(t) = 0;
				 }
			}
	};
```

Next, we could add a function to to our ```update_Penguins``` class that includes the following lines:

```c++
int *max_ages = raw_pointer_cast(species->demeParameters->get_vector_ptr("Maximum_Age"));
death_from_old_age death_from_old_age(rand_ptr, max_ages);
```
The function ```get_vector_ptr()``` in the demeParameters class returns a pointer to the float vector storing the maximum ages of individuals permissible in the different demes. This then gets passed as an argument to the death_from_old_age functor. The relevant segments in our corresponding new deme_config.txt file now looks like:
```c++
				number_of_parameters = 5
				parameter_names = ["M_reproductive_advantage", "F_reproductive_advantage","TARGET_CROWN_COLOR","CROWN_COLOR_DECAY", "Maximum_Age"]

				deme_specifications = (
										{#deme0
										M_reproductive_advantage = 1.0
										F_reproductive_advantage = 1.0
										TARGET_CROWN_COLOR = 0.808080
										CROWN_COLOR_DECAY = 0.1
										Maximum_Age = 5.0
										},
										{#deme1
										M_reproductive_advantage = 1.0
										F_reproductive_advantage = 1.0
										TARGET_CROWN_COLOR = 0.809120
										CROWN_COLOR_DECAY = 0.01
										Maximum_Age = 25.0
										}
									  ) 
```

However, we're not done yet - we still need to add a for_each() instruction. We'll leave this part as an exercise for the reader.

2. Study the code in [Examples/CoevolvingSpecies](/Examples/CoevolvingSpecies), particularly lines 15-32 of migrate_coevolvingSpecies.h, migrate_CoevolvingSpecies.cu, MigrationBehavior.cu and line 36 of coevolution_Simulator.cpp. Now see if you can implement a version of the tutorial WITH migration, to see what happens to the trajectory of alleles in the two populations. Hint: write a routine to change the individual's demes at random to simulate panmixia.


<a name="setup_issues">

## Some possible issues with getting sPEGG set up
</a>

Depending on your system's configuration, there may be some manual steps you would need to keep in mind when setting up **sPEGG**. Unfortunately, some of these issues are the result of a few of the pre-requisites that I've found don't always work "out of the box". On Ubuntu, using some of the default packages through 

```sh
$ sudo apt-get <package>
```
can cause problems as of late 2015/early 2016. Here are some alternative approaches to consider until the situation is improved.

### Libconfig issues
There is a debian package for libconfig, but if you manually download and install libconfig from its [webpage](http://www.hyperrealm.com/libconfig/) and run the configure and make scripts that come with the download as root, the executables etc... are by default saved to the ```/usr/local``` directories. You will need to update your library path, as described [below](#library_path) to point to the relevant directory for the .so file (often, ```/usr/local/lib``` is sufficient).

### CUDA issues
Installing CUDA and its compiler nvcc on Linux can be a pain, although this is rapidly improving on the [distros for which NVIDIA supports CUDA](http://docs.nvidia.com/cuda/cuda-getting-started-guide-for-linux/#system-requirements). Part of the problem is that the drivers that come with the distro to handle NVIDIA GPUs (e.g., Nouveau) might not be compatible with NVIDIA's CUDA. Here's a [tutorial](http://http://www.r-tutor.com/gpu-computing/cuda-installation/cuda7.5-ubuntu) to get you started for Ubuntu. Generally, analogous approaches to installing CUDA should work on all [supported distros](http://docs.nvidia.com/cuda/cuda-getting-started-guide-for-linux/#system-requirements), but I have occasionally had some issues with getting CUDA to work with (currently) unsupported distros, such as Linux Mint. 

### CPU version on non-CUDA enabled machines
The CPU version of **sPEGG** should work without a CUDA-capable device, but this has not been extensively tested. Interested readers are encouraged to see if they can get a [basic thrust example](https://github.com/thrust/thrust/wiki/Quick-Start-Guide) running on such systems first.

### Non-root issues
Unfortunately, as best we can gather the CUDA issues noted above will require root access to resolve. If you don't have root access to your computer and can't have the administrator install libconfig and libgsl, a version of the Makefile for linking to these libraries when you aren't root is currently in the works.

Finally, if you want to call nvcc as you go, you would need to install CUDA in your .bashrc file (or the equivalent for your distro). This assumes you are installing CUDA v.7.5.

<a name="library_path">

```sh
export CUDA_HOME=/usr/local/cuda-7.5 
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

PATH=${CUDA_HOME}/bin:${PATH} 
export PATH 
```
</a>

## License

GPLv3 © Kenichi Okamoto
