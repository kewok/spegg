# sPEGG


**sPEGG** is a simulator and framework for studying eco-evolutionary dynamics. Using forward-time population genetics, **sPEGG** allows researchers to study the interplay between ecological and evolutionary processes in multispecies metacommunities.

 **sPEGG** harnesses the massive parallelization potential of relatively economical General-Purpose Graphics Processing units (GPGPUs) installed on modern desktops. Depending on the complexity of the ecological and evolutionary components of your model, **sPEGG** can provide performance gains on your desktop comparable to gains obtainable on a small-to-medium sized cluster.

<center> Version 0.5 </center>

**sPEGG** continues to be under active development. Should you have any comments or if there are aspects of the usage that are unclear, please feel free to note any [issues](https://github.com/kewok/spegg/issues) you might have, to help make **sPEGG** more useful :smile:.

---


Currently, **sPEGG** has only been tested on Linux (debian-based distributions to be specific). If you encounter platform-specific commands that need to be added while porting **sPEGG** to other operating systems, please [create an issue](https://github.com/kewok/spegg/issues) and I will update the documentation. Or, if you can get it to work reliably, consider becoming a [contributor](https://github.com/kewok/spegg/contributors)!

## New updates for Colab users

If users for some reason are not able to afford a GPU backend, please don't just clone the repository first. Instead, first, download the [firstStepSpegg.ipynb](firstStepSpegg.ipynb) file. On Colab, create a new notebook, then choose `File (top-left corner) > Upload notebook > Upload > Browse (right in the center of the window) > location of the firstStepSpegg.ipynb file you just downloaded`. Then, the notebook should be correctly imported into your current Colab notebook. Set the `PATH_TO_YOUR_SPEGG_IN_DRIVE` variable specifying where in your GoogleDrive you want to store this repository. Then, press the play buttons from the top until you hit the end of the notebook. When these steps are done, you should find your spegg directory in your specified location in your GoogleDrive. Follow step 8 in the `firstStepSpegg.ipynb` when you come back and work with spegg codebase again.

<a name="prereqs">

## Pre-requisites
</a>

  * Linux OS: [Ubuntu-Linux](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/#ubuntu-installation)

  * The latest [NVIDIA CUDA](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/#system-requirements) for your distro; appropriate drivers and the required thrust library should come with this.
  
  * [gcc](https://gcc.gnu.org/)
  
  * [libconfig](https://hyperrealm.github.io/libconfig/)
  
  * [ConfigFile](http://ai.stanford.edu/~gal/Code/FindMotifs/ConfigFile.h)(included)
  
  * CUDA-enabled GPU (generally this will require [compute capability >= 2.0](http://en.wikipedia.org/wiki/CUDA#Supported_GPUs))
  
### In addition, the CPU version requires

  * [GNU Scientific Library](http://www.gnu.org/software/gsl/)
</a>  
Please see here for some [possible issues](Documentation/Tutorial.md#setup_issues) you may run into when setting up **sPEGG**.

Once you have everything in place, you can run your first **sPEGG** simulation!

***
##  Quick-start example: Simulating genetic drift </h4> </center>

We'll use the code for the simulation that is created from the [tutorial](Documentation/Tutorial.md) to simulate genetic drift in two demes of 10000 individuals each across 100000 generations on the GPU.

Open a terminal (ctrl+Alt+T in recent ubuntu versions) and navigate to your project directory.

```sh
$ cd </path/to/My_sPEGG_project>
```

Using [git](http://git-scm.com/), clone sPEGG's git repository
```sh
$ git clone --recursive https://github.com/kewok/spegg/
```

Enter the main code base for sPEGG 

```sh
$ cd spegg
```
in your terminal.

Now, simply type 

```sh
$ make
```
in the terminal. Depending on your system, this might take awhile.

Now, navigate into the directory:

```sh
$ cd Examples/Tutorial_Simulation
```
and type

```sh
$ make
```

in the terminal. If everything has been installed correctly, this will create an executable (a.out). You can run it as:

```sh
$ ./a.out
```

However, the genetic drift simulator you downloaded is only configured to run for a single generation. Open the file Simulation.conf using a text editor (e.g., gedit) and change this line:

```
n_timesteps = 1
```

to something like:

```
n_timesteps = 10000
```

Now launch  the **sPEGG** simulation by running this executable you created from the terminal via 

```sh
./a.out 
```

If all went well, you should have output that looks something like:

```
The average genotypes at locus 1 for the two demes are:
-0.596273 0.99806 
```

The exact numbers may very well differ, however depending on your hardware.

This just gives you a flavor for running a **sPEGG** simulation. If you are interested in simulating more biologically relevant models, have a look at the other examples described  ([here](Documentation/1603.09255v1.pdf)). To configure and run those models, you'll want to modify the deme_config.txt and environment_config.txt files as well.

For those of you who'd like to build your own **sPEGG** model, have a look at our more in-depth [Tutorial](Documentation/Tutorial.md).


## License
**sPEGG** is released under the terms of the GNU Public License v3.0. Please refer to the LICENSE file. © Kenichi Okamoto, 2016.
