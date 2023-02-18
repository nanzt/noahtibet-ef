# Modified Noah LSM v3.4.1 associated with a manuscript submitted to Earth's Future(2022EF002652: Qinghai-Tibet Plateau permafrost at risk in the late 21st century)

## A short note

Key improvements in the modified Noah LSM codes include a modified thermal roughness scheme for sparse vegetation typical of the QTP, improved parameterization of thermal and hydraulic conductivities to account for coarse-grained QTP soils and ground-ice, and an expansion of the simulation depth to below 15 m that allows for vertical soil heterogeneity.

For more technical details, please refer to:

* Wu X, Nan Z*, Zhao S, Zhao L, Cheng G. Spatial modelling of permafrost distribution and properties on the Qinghai-Tibet Plateau. Permafrost and Periglacial Processes. 2018, 29(2): 86-99. doi:10.1002/ppp.1971. [https://onlinelibrary.wiley.com/doi/full/10.1002/ppp.1971](https://onlinelibrary.wiley.com/doi/full/10.1002/ppp.1971)

This release is the version used for the simulation in the work submitted to Earth's Future. We are continuourly working to improve the LSM more suitable for permafrost simulation on the QTP. We decided to disclose these codes in the hope of encouraging sharing by the permafrost science community.

It is based on the Unified Noah LSM 3.4.1, which was originally developed by scientists from Research Applications Laboratory, NCAR. More information as well as the source codes of Noah LSM 3.4.1 can be found at [its official site]( https://ral.ucar.edu/solutions/products/unified-noah-lsm). 

If you have made modifications based on this version or published a paper using this version, please share your good news with us (nanzt@njnu.edu.cn). Unfortunately, since our time is limited, we cannot provide any technical assistance. Please check our group site for further improvements to this version. 

Permalab team
([https://permalab.science](https://permalab.science))
May 7, 2022


## Contents

In this folder, three subfolders：“Modified Noah LSM v3.4.1 code”，“Model_forcing”，and “Model_output” are included, alongside with this "README.md" file.

The codes of the modified Noah LSM v3.4.1 can be found in folder `Modified Noah LSM v3.4.1 code`

Source code: all `*.F` files, and the main driver code is `simple_driver.F`. 

Parameters files：
* General parameters -- GENPARM.TBL	
* Soil Parameters	-- SOILPARM.TBL
* Vegetation Parameters -- VEGPARM.TBL
* Urban Parameters -- URBPARM.TBL
	
Txt file `forcing_example.txt` in `Model_forcing` is an example input for drivering the LSM.

Txt file `forcing_example.txt` in `Model_output` directory is the corresponding output of the modified Noah LSM.


## Notes on compiling the source codes

### Set up compilation environment

Noah LSM runs in a Linux environment.

If you are using Windows 10, you can either enable [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/about) or install a virtual machine using [Oracle VirtualBox](https://www.virtualbox.org).


### Compile the modified Noah LSM source codes

Edit user_build_options to choose compiler, compiler options, libraries, etc. 

Adjust compiler options and library paths as necessary from the default settings provided in user_build_options.

Invoke "make" to compile. 

Once the model successfully compiled, the excecutable file `driver.exe` is generated in the `Modified Noah LSM v3.4.1 code` directory.

### Run the model

Copy `driver.exe` and four `*.TBL` parameters files from the `Modified Noah LSM v3.4.1 code` directory to the `Model_forcing` directory.

Run:  "driver.exe forcing_example.txt"

Once the model is successfully run, a txt file `forcing_example.txt` is produced and saved in the `Model_output` directory.

