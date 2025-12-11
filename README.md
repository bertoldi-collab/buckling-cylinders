# buckling-cylinders

Repo with all code and data used for "Complex deformation in soft cylindrical structures via programmable sequential instabilities (2024)" (DOI 10.1002/adma.202406611). Main file `cylinder_obj.py` generates Abaqus models and input files. Tested on Abaqus/2019. All units are assumed to be in (mm, MPa). Authors of this code are: Helen Read, originally David Melancon, Antonio Elia Forte, and Ahmad Zareeeei. If you have questions, please contant me directly at `heread@g.harvard.edu` or my advisor (Katia Bertoldi) at `bertoldi@seas.harvard.edu`.

## How to use this
`cylinder_obj.py` wants to be called from a sub-folder one level down. Additionally, post-processing scripts expect a folder called `data_out` that is in the main repo to exist (from your sub-folder where your specific job creating scripts are `../data_out` should be a valid address). A python file in your job-creating sub-folder should begin with

```
import sys
sys.path.append('../')

from cylinder_obj import *
from geo_prop import *
```

Importing `cylinder_obj` imports the classes that create Abaqus models for these thin-walled cylinders. Importing `geo_prop` imports some named tuples that contain all geometric variables and are used to initialize instances of the cylinder classes. This might look like:

```
geo_props_use = geoProps(radius, height, cap_width, thickness, young_modulus)
#this is the simplified geoProps, assuming the thickness of the cylinder is constant

proj_name = 'this_is_a_project_name'
cylinder_example = full_shell(project = proj_name, simpProps = geo_props_use, imperfection = 0.002)
#we use simpProps because geo_props_use is type geoProps, not geoPropsFull
```

Now that the abstract cylinder object is created, we can call different methods to make various types of abaqus analyses. These should all return job names, and you can use `run_job()` with the returned job name to run the job. (Note: `run_job()` is not a method in a class, and is instead defined outside any class. Here are a list of the types of simulations that are currently written

- `make_nonlin_multi_buckle()`
-  `make_riks_model()`
-  `make_nonlin_model()`
-  `make_static_dyn_model()`
-  `make_linear_model()`
-  `run_linear_model()`
-  `make_force_buckling_model()`

An example call to make a job and run that job might look like this
```
cylinder_example.run_linear_model()
#this actually makes and runs a linear buckling analysis; we need the buckling modes to run any nonlinear analysis so that we can have an imperfection
#notably, this generates a .fil file that is needed to add an imperfection into the model

#this makes a dynamic implicit model that removes 30% of the volume of the cylinder
jname_nonlin = cylinder_example.make_nonlin_model(temp_mult = 0.3)
run_job(jname_nonlin, num_cores = 4)
```

Now, we might want to post-process our generated `.odb` file. Here's a list of post-processing scripts

- `post_process_multi_buckle()`
- `post_process_multi_pv()`
- `post_process_pv()`
- `post_process_multi_contraction_twist()`
- `post_process_contraction_twist()`
- `post_process_bend_angle()`
- `post_process_lin_centernodes()`
- `post_process_num_folds()`
- `post_process_centernodes()`

and an example of how we might use one or two
```
cylinder_example.post_process_pv()
#this saves a file with pressure-volume data to ../data_out/this_is_a_project_name_pcav_cvol.txt

cylinder_example.post_process_contraction_twist()
#creates similar txt file in data_out
```

Hope this helps! Functions that I actually used for this paper should be reasonably documented with what various arguments are doing.
