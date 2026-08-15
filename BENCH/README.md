Benchmarking scripts
====================

The given directory aimed at providing some basic benchmarking scripts to
compare CROCO on various tests configurations (physical options, parallelism, 
compilation option, compilation type of build...).

Concepts
--------

First you can give a look to the `config.jsonc` file. You will see the concepts
of :

- `cases` : The cases to run (BASIN, SOLITON...) and some case size changing
  the grid size (eg. BASIN-LARGE, BASIN-HUGE).
- `variants`: The parallelization mode (openmp, mpi, gpu) with specification on
  the ressources to use because CROCO needs to be rebuilt for every case du to
  some currently static definitions in the CROCO code on the use of those 
  ressources.
- `hosts`: You will find some variables to adapt the running to the host you
  are running on.

Note: you can force the host to be selected if needed by using `--host NAME`.

Dependencies
------------

As it will build CROCO for you, you will need to have the basic dependencies
available (see main README.md to know what is required).

For bench script you will need a python envrionment with 
some required packages listed in requirements.txt.

```sh
# if you already have an activated venv
pip install -r requirements.txt
```

Otherwise you first needs to create a `virtualenv` :

```sh
# only the first time
python3 -m venv venv

# for every terminal used to run the bench script
source venv/bin/activate

# only once
pip install -r requirement.txt
```

Running
-------

The simplest way to run is :

```sh
# simplest way
./bench-croco.py

# use your own config file
./bench-croco.py -c config-mine.jsonc
```

It will run the `cases` and `variants` combination defined in 
`default_selection` of the config file.

You can also choose by hand some of those :

```sh
# just a single
./bench-croco.py --case BASIN --variants sequential

# or many
./bench-croco.py --case BASIN,BASIN-LARGE --variants sequential,openmp-8,mpi-8
```

Build & results
---------------

The script will build CROCO for each `case`/`variant` in the workdir (`./rundir`
by default).

The results are stored in `./results` (by default) and a graph will be produced
for each case.

By default the benchmarking script will run CROCO 1 time for each case. 
You can change the number of times each run is done to plot
the error bars, for example 4 times by using :

```sh
./bench-croco.py --runs 4
```

The script will keep the results of each run by `{hostname}-{date}` but you can
also ask the script to add a given name if you are making different tries
`{title}-{hostname}-{date}` .

```sh
# while you are developping your optimizations (trash dir)
./bench-croco.py --title dev
# Produce clean final perf ready for paper
./bench-croco.py --title final-perf
```

Partial re-reun. If you make a new run on the same machine not running all the
cases the script will search the missing results in the previous result
directories in order to aggregate no the charts (with same `title`). This is
usefull in dev process to compare without re-running everything.

Configuration
-------------

The configuration file can be provided under `json` which also supports 
comments. It can also be splitted in subdirectory mode by merging several 
subfiles.

In the main `config.jsonc` file, the inclusion is described by the `imports`
section.

Whildcards
----------

You can also use whildcards for `case` and `variant` to select them :

```sh
./bench-croco.py --case 'BASIN*' --variants 'openmp-*'
```

Meta variants
-------------

You will notice in the `config.jsonc` file the definition of `meta_variants`
which allow you to define sets of variants which you might use often for your
run.

Then you can just use them easily:

```sh
./bench-croco.py --case BASIN --variants @cpu
# you can remove some of the elements by using '-' :
./bench-croco.py --case BASIN --variants @cpu,-openmp-4
```

Meta case
---------

If you want to simply run all available cases and all variants, you can just:

```sh
./bench-croco.py --case @all --variants @all --auto-skip
# you can remove some of the elements by using '-' :
./bench-croco.py --case @all,-BASIN-LARGE --variants @all --auto-skip
# you can also use wildcard in the expr for removing :
./bench-croco.py --case '@all,-BASIN-*' --variants '@all,-openmp-*' --auto-skip
```

Validation report
-----------------

If you want to run `bench-croco.py` for validation and get a summary report
on what is valid or not, you can simply :

```sh
./bench-croco.py --case @all --variants @all -r --report
```

You can select what is the reference variant to compare to (default is 
`sequential`) :

```sh
./bench-croco.py --case @all --variants @all -r --report --compare-to openmp-2
```

You will get an output like :

```plain
------------------------------- REPORT SUMMARY ---------------------------------
CASE-VARIANT                      Build     Run           Check   
SOLITON-sequential                None      True          True    
SOLITON-openmp-2                  None      True          True    
SOLITON-LARGE-sequential          None      True          True    
SOLITON-LARGE-openmp-2            None      True          True 
```

Plus the results stored in `status_report.json`.

Reference data
--------------

For some reason, the reference sequential run can cost a lot to produce if the
case is large. In order not to redo it all the time, you can use two approach :

The simpler one consist in making the run a first time to get the reference data
generated then to run without `sequential`:

```sh
# first time
./bench-croco.py --case @all --variants @all -r --report
# others
./bench-croco.py --case '@all,-*sequential' --variants @all -r --report
```

The other approach consist in producing a reference directory which can be
zipped and uploaded somewhere. This is possibly usefull for CI integration
to ensure non regression even on the `sequential` run.

```sh
# first time
./bench-croco.py --case @all --variants @all -r --build-ref ./croco-refs-2024-05
# others
./bench-croco.py --case '@all' --variants @all -r  --report --use-ref ./croco-refs-2024-05
```

Auto-skip
---------

When arriving on an host, `bench` might try to run some `variants` which 
requires more ressources than available (`cores` or `gpu`). You can ask to 
automatically skip those configs :

```sh
./bench-croco.py --case @all --variants @all -r --auto-skip
```

Performance
-----------

In order to perform correct performance measurement you can ask `bench-croco.py`
to repeat several times the run to get some error bars.

By default it makes only a single run per case/variant.

You can set it to 4 by :

```sh
./bench-croco.py --case BASIN --variants @all --runs 4
```

Variables
---------

You might notice in the configuration file the use of some variables to be
replaced at parsing time on the form of : `{tuning.gnu}`.

It avoids duplicating too many entries in the file.

You will mostly find them at two places : 

- On the cases to unpack king of template not to fully duplicate all cases for
  OpenMP and MPI.
- To inject some values from `host` into the `jobcomp-configure` call.
- To inject the mpi extra options in `command_prefix`.

**Note**:

- For the template base system you can recurse the template var (`{thread}`, `tasks`)
  to build the final var name (as it is done in the default config).
- For the rest you cannot currently recurse even it can be trivial to patch the code to enable it
  (just that I didn't had time to validate it really works by doing it).
  For this case, you have acces to two base keys : `tuning` & `case`
  corresponding to the corresponding par of the config file you are running
  (host & current case).

Configuration patching
----------------------

You will notice in the cases configuration files (`config.d/cases/*`) that there
is for most of them a patching part which explain how to modify the CROCO
configuration files to get the case setted up.

The patching is like this :

```json
"patches": {
  // Define a file to patch
  // Note : it can also contains an array of patch elements if needs to change
  //        several lines.
  "param.h": {
    // Patch mode
    "mode": "replace",
    // What to replace (expected to be present)
    "what":    "      parameter (LLm0=60,   MMm0=50,   N=10)",
    // By what to replace
    "by":      "      parameter (LLm0=200,  MMm0=200,  N=50)",
    // Description to be printed in the terminal when applied
    "descr":   "Set mesh size to (200 x 200 x 50)"
  },
},
```

Imports
-------

You will also notice at the head of the configuration file the `imports` section
which allow to split the config file in multiple sub files and merge them
at loading time.

This is usefull for exemple for the `hosts` and `cases` definitions as it 
permits a user to define its own hosts without requirement to patch the files 
tracked by git and to keep its definition for himself until he finally decide 
to commit it.

Authors
-------

- First implementation by Sébastien Valat (INRIA / LJK) - 2023 - 2024
- Adaptation for CROCO gitlab-ci by Rachid Benshilla (CNRS), Gildas Cambon (IRD) 
  and Solène Le Gac (IFREMER) - 2024 - 2025
