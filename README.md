# sPEGG


**sPEGG** is a simulator and framework for studying eco-evolutionary dynamics. Using forward-time population genetics, **sPEGG** allows researchers to study the interplay between ecological and evolutionary processes in multispecies metacommunities.

 **sPEGG** harnesses the massive parallelization potential of relatively economical General-Purpose Graphics Processing units (GPGPUs) installed on modern desktops. Depending on the complexity of the ecological and evolutionary components of your model, **sPEGG** can provide performance gains on your desktop comparable to gains obtainable on a small-to-medium sized cluster.

<center> Version 0.5 </center>

**sPEGG** continues to be under active development. Should you have any comments or if there are aspects of the usage that are unclear, please feel free to note any [issues](https://github.com/kewok/spegg/issues) you might have, to help make **sPEGG** more useful :smile:.

---


Currently, **sPEGG** has only been tested on Linux (debian-based distributions to be specific). If you encounter platform-specific commands that need to be added while porting **sPEGG** to other operating systems, please [create an issue](https://github.com/kewok/spegg/issues) and I will update the documentation. Or, if you can get it to work reliably, consider becoming a [contributor](https://github.com/kewok/spegg/contributors)!

<a name="prereqs">
## Pre-requisites
</a>

  * Linux OS

  * The latest [NVIDIA CUDA](https://developer.nvidia.com/cuda-downloads) for your distro; appropriate drivers and the required thrust library should come with this.
  
  * [gcc](https://gcc.gnu.org/)
  
  * [libconfig](http://www.hyperrealm.com/libconfig/)
  
  * [ConfigFile](http://ai.stanford.edu/~gal/Code/FindMotifs/ConfigFile.h)(included)
  
  * CUDA-enabled GPU (generally this will require [compute capability >= 2.0](http://en.wikipedia.org/wiki/CUDA#Supported_GPUs))
  
### In addition, the CPU version requires

  * [GNU Scientific Library](http://www.gnu.org/software/gsl/)
</a>  
Please see the bottom of this page for some [possible issues](#setup_issues) you may run into when setting up **sPEGG**.

Once you have everything in place, you can run your first **sPEGG** simulation!

***
##  Quick-start example: Simulating genetic drift </h4> </center>

We'll use the code for the simulation that is created from the [tutorial](https://github.com/kewok/spegg/Examples/Tutorial) to simulate genetic drift in two demes of 10000 individuals each across 100000 generations on the GPU.

Open a terminal (ctrl+Alt+T in recent ubuntu versions) and navigate to your project directory.

```sh
$ cd </path/to/My_sPEGG_project>
```

Using [git](http://git-scm.com/), clone sPEGG's git repository
```sh
$ git clone https://github.com/kewok/spegg/
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

```sh
n_timesteps = 1
```

to something like:

```sh
n_timesteps = 10000
```

Now launch  the **sPEGG** simulation by running this executable you created from the terminal via 

```sh
$ ./a.out 
```

If all went well, you should have output that looks something like:

```sh
The average genotypes at locus 1 for the two demes are:
-0.596273 0.99806 
```

The exact numbers may very well differ, however depending on your hardware.

This just gives you a flavor for running a **sPEGG** simulation. If you are interested in simulating more biologically relevant models, have a look at the other examples described  ([here](arXiv.org)). To configure and run those models, you'll want to modify the deme_config.txt and environment_config.txt files as well.

For those of you who'd like to build your own **sPEGG** model, have a look at our more in-depth [Tutorial](https://github.com/kewok/spegg/Documentation/tutorial.md).

## License
**sPEGG** is released under the terms of the GNU Public License v3.0. Please refer to the LICENSE file. © Kenichi Okamoto, 2016.
