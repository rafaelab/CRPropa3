# Plugin template

To create your own CRPropa modules in C++ we recommend using plugins. Plugins are small programs that can be installed as a separate modules. 
Here, we provide a template to create such a plugin. 

## General structure of the template

All include should be in `include/myPlugin/`, so that the `myPlugin` folder can be moved to the install location with
minimum risk of overwriting any existing files in the install location.
Additionally a file `include/myPlugin.h` should be provided where all header files in `include/myPlugin/` are included so the user has an easy way of including the whole project.

The source files should be placed in `src/`.

All files related to the python package (so the swig header and eventual python tools) should be placed in `python/`.

Finally you should create tests that test if your code is working (technically and physically) to `test/`.

## Adjusting custom module name

When creating your own module you need to adjust the name of your module at multiple occasions:
- File and folder names, such as `include/myPlugin/` and `include/myPlugin.h`
- The project name in your CMakeLists.txt (you only need to change the project name on top, everything else should change automatically)
- `python/myPlugin/__init__.py`: The directory name and the content of the init-file have to be changed: `.myModule` to `.<MyModuleName>`
- `myPlugin.i`: at two positions the header file is listed. The lines have to be adapted accordingly. 

# Installation of a plugin

There are some ways to install your plugin.
- The first one would be to use an existing crpropa version and let it be found automatically by cmake.
- The second one would be to use the option to use a builtin crpropa, this option first pulls the current master repository and builds it alongside your plugin.

Regardless which of the options above you use, you will then need to configure your plugin and start building:

```sh
cd /path/to/your/pluginfolder
mkdir build && cd build

cmake .. -D CMAKE_INSTALL_PREFIX=/your/install/location -D USE_CRPROPA_BUILTIN=OFF
cmake --build .
cmake --install .

```

Set `USE_CRPROPA_BUILTIN=ON` if you want to install crpropa alongside your project.
Instruction how to install crpropa yourself can be found [here](https://crpropa.github.io/CRPropa3/pages/Installation.html).


## optional testing

You can simply test your plugin with:

```sh
ctest --output-on-failure --repeat until-pass:3
```