### Math Functions
*[ ] f64 acos(f64 x)
*[ ] f64 asin(f64 x)
*[ ] f64 atan(f64 x)
*[ ] f64 atan2(f64 y, f64 x)
*[x] f64 cos(f64 x)
*[ ] f64 cosh(f64 x)
*[x] f64 sin(f64 x)
*[ ] f64 sinh(f64 x)
*[x] f64 tan(f64 x)
*[ ] f64 tanh(f64 x)
*[x] f64 exp(f64 x)
*[ ] f64 frexp(f64 x)
*[ ] f64 ldexp(f64 x)
*[x] f64 log(f64 x)
*[x] f64 log10(f64 x)
*[ ] f64 fsplit(f64 x, f64 *int)
*[x] f64 pow(f64 x, f64 y)
*[x] f64 sqrt(f64 x)
*[ ] f64 ceil(f64 x)
*[x] f64 fabs(f64 x)
*[ ] f64 floor(f64 x)
*[ ] f64 fmod(f64 x, f64 y)
*[x] f64 csc(f64 x)
*[x] f64 sec(f64 x)
*[x] f64 cot(f64 x)
* f32 and f64 versions of these functions are being implemented.
* f80 versions are not possible with x87 instructions as far is i know so an other method should be used for that later
* add newtons method to the log functions to improve precision

### Other Features
*[ ] add custom types and operations to gcc toolchain https://hsandid.github.io/posts/risc-v-custom-instruction/
*[ ] add vector and matrix math in asm using vectored registers (xmm, ymm, zmm, etc..) (and) or on the gpu