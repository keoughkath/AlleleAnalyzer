### Requirements

All noted as "specific" only need to be downloaded for that specific use-case. Most users won't need those extra packages. 

* Operating systems: Unix, OSX
* Python 3.x, 2.7.x for CRISPOR portion

#### Python modules:

These should be installed using your favorite Python packages manager (i.e. conda/pip).

* Pandas (>= 0.21.0)
* Numpy (>= 1.13.3)

**set-cover specific**
* pulp (>= 1.6.8)

**preprocessing specific**
* pyfaidx (>= 0.5.1)

#### External binaries (preprocessing-specific)
* bcftools (>= 1.5)

### Download the software

Clone the latest version from Github: 

`git clone https://github.com/keoughkath/ExcisionFinder`

### Update ExcisionFinder

As ExcisionFinder gains new features and bugs are solved, changes will be made. To keep your installation current, cd into the installation directory and pull the latest version: 

`cd ExcisionFinder`

`git pull`