# Code_Optimization
A program that aims to implement optimization on someone else's code (original_code.c) and see the improvement in runtime. <br/>
The original code was `original_code.c` and the new code after optimizations is `optimized_code.c`

## Code Optomizations ##
For optimize the code and analyzed him I used GPROF - profiling tool. <br/>
Some examples of the optimizations that i did:
* `Code motion`
* `Using registers` - reading and writing registers is much faster than reading / writing to memory.
* `Using macro expression` - write the math operation on the expression as macro instead of doing it during compilation time (that call to math libary...).
* 'Etc'

### The improvment in runtime is: <br/>
Original runtime code: <br/>
![picture alt](/Original_Time.png "Original Time") <br/>
Optimized run time code: <br/>
![picture alt](/Optimized_Time.png "Optimized time")
