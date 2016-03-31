---
layout: page
title: Download
permalink: /download/
---

You can download the latest version of GeoPDEs in the following link: 

## [Download GeoPDEs 3.0.0]()

# File contents

* The GeoPDEs package.
* The NURBS package.
* A technical report, explaining the design of GeoPDEs.
* A README file, with detailed installation instructions.
* For MATLAB users, the mex files of the NURBS toolbox.

# Installation instructions (quick guide)

**Octave users**

Type at the Octave prompt  the three following commands:

1. pkg install nurbs-_version_.tar.gz
2. pkg install geopdes-_version_.tar.gz
3. pkg load geopdes

**Matlab users**

Uncompress and untar the files containing the NURBS package (nurbs-_version_.tar.gz) and the GeoPDEs package (geopdes-_version_.tar.gz). 

Add the two folders, with all their subfolders, to the Matlab path.

More detailed installation instructions are found in the README file.

# <a name="examples"></a>Getting started: running the examples

GeoPDEs contains a set of simple examples, that can be run with the commands:

* geopdes_base_examples
* geopdes_elasticity_examples
* geopdes_fluid_examples
* geopdes_maxwell_examples
 
These will guide you through a set of menus to choose some examples for solving different model problems.

# GeoPDEs in GitHub

The source code of GeoPDEs is now available in [GitHub](https://github.com/rafavzqz/geopdes). Older versions of GeoPDEs are also available in the repository.