# Platforms
> ## x86_64
> this is the target architecture for this project and thus the mathfunctions are x86_64 calling convention complaint.
> ### Discrepancies with x86_32:
>> * #### Floating points with the x86_32 calling convention
 >>  the x86_32 architecture does not support **XMM** registers.
 >>  This means that most of the code dealing with floating points wont work because the x86_64 calling convention passes its floating points in **XMM** registers.
 >>  A compiler flag (**"-msoft-float"**) can be added to force gcc to pass floats on the stack like x86_32 does this will not make up for the fact that we need **XMM** registers and it would make the code slower in general so this is **NOT** added.

> ## Resources
> ### SSE:op
>> * https://students.mimuw.edu.pl/~zbyszek/asm/en/instrukcje-sse.html
