# STD_Library
Standard library using ASM/C for ASM/C/C++/Rust

## Development Structure:

> ### /src/, /inc/
> library source files

> ### /int/, /bin/
> .(s)o output directory
> 
> in bin/static a .h file equivalent .o file can be found for static linking

> ### /test/
> testing files and programs

## Preferred Clion Config

> ## build config
> #### system dependencies
> * make
> * gnu toolchain

> ## testing tool config
> "test" is a config that will first run the aforementioned build config after which it will compile and run the test program "/test/main.c"